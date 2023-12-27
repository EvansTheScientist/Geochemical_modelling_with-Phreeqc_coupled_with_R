# Time-stamp: "Last modified 2022-02-21 17:24:37 delucia"
library(RedModRphree)
options(width = 116)


##' This function performs a combinatorial analysis, creating all
##' possible PHREEQC models with a number "len" of equilibrium phases
##' across a pool of primary and secondary minerals.
##'
##' Requires packages doParallel and foreach (for parallelization),
##' but they should already be installed as dependency of RedModRphree
##' @title Perform combinatorial analysis of a given pool of
##'     equilibrium minerals
##' @param initsol the base PHREEQC script
##' @param primary character vector containing minerals which are
##'     inserted in the models as primary minerals
##' @param secondary character vector containing minerals which are
##'     inserted in the models as secondary minerals (i.e., initial 0
##'     amount)
##' @param len number of phases in each model
##' @param procs integer, how many CPUs can we use? Defaults to 4
##' @return a list containing the parsed results (list of many blocks)
##' @author MDL
DoCombSim <- function(initsol, db, primary, secondary, len, procs=4L) {
    require(parallel)
    require(doParallel)
    require(foreach)
    
    ## create all combinations of primary and secondary phases of length len
    phases <- c(primary, secondary)
    combs <- combn(phases, len, FUN = NULL, simplify=FALSE)

    cat(":: Going to do ", length(combs), " simulations\n")

    ## create the phreeqc scripts
    .addPhaseComb <- function(x) {
        inp <- initsol
        for (phase in x) {
            inp <- AddProp(inp, name=phase, values=ifelse(phase %in% primary, "0.0 2", "0.0 0"), cat="pphases")
        }
        return(inp)
    }

    ## apply this function to all entries
    biginp <- lapply(combs, .addPhaseComb)

    ## workhorse function to run simulations
    .runPQC <- function(input) {
        phreeqc::phrSetOutputStringsOn(TRUE)
        phreeqc::phrRunString(input)
        tmpout <- phreeqc::phrGetOutputStrings()
        res <- RedModRphree::ReadOut(tmpout)[[1]]
        return(res)
    }

    if (procs > 1) {
        if (Sys.info()[["sysname"]]=="Windows") {
            ThisRunCluster <- parallel::makePSOCKcluster(procs)
        } else {
            ThisRunCluster <- parallel::makeForkCluster(procs)
        }

        doParallel::registerDoParallel(ThisRunCluster)
        cat(":: Registered default doParallel cluster with ", procs, "nodes\n")
        parallel::clusterCall(cl=ThisRunCluster, phreeqc::phrLoadDatabase, db)

        res <- foreach::foreach(i = seq_along(biginp)) %dopar% 
            .runPQC(biginp[[i]])
        cat("[ DONE ]\n")
        parallel::stopCluster(ThisRunCluster)

    } else {
        ## revert to sequential computation
        cat(":: Firing up PHREEQC onsingle CPU...")
        res <- lapply(biginp, .runPQC)
        cat("[ DONE ]\n")
        
    }
    
    return(res)
}

#############################################################

##' Computes a specific metric allowing for selection of the
##' components to be included
##'
##' @title Compute metric selecting components
##' @param data the matrix or data.frame containing all the results
##'     from the PHREEQC simulations. Its columns need to be named!
##' @param target the named vector with the target concentrations
##' @param FUN the name of the metric function
##' @param comp optional, a char vector with the names of the
##'     components. If unspecified, all components are selected
##' @param ... further parameter passed to FUN, such as "na.rm"
##' @return numeric vector with the computed metric
##' @author Marco
ComputeMetric <- function(data, target, FUN="rmse", comp=colnames(data), ...){
    ## find the metric function
    .Fun <- match.fun(FUN)

    ## retain only the columns given as argument
    tmp <- subset(data, select = comp)

    cvec <- target[comp]
    ## compute stuff using apply
    res <- apply(tmp, 1, function(x) .Fun(cvec, x, ...))
    return(res)
}


##### Filtering

Filter <- function(lin, delta=0.5) {
    excluded <- sapply(lin, function(x) any(abs(x$pphases$delta)>delta))
    out <- lin[which(!excluded)]
    return(out)
}

Filter2 <- function(lin, delta=0.5) {
    excluded <- sapply(lin, function(x) any(abs(x$pphases$delta)>delta))
    return(which(excluded))
}

FilterAll <- function(lin, delta=0.5) {
    retain <- sapply(lin, function(x) all(abs(x$pphases$delta)<delta))
    out <- which(!retain)
    return(out)
}

## Some metrics
## ## root mean square logarithmic error

## root mean square error
rrmse <- function(y_true, y_pred, na.rm=TRUE)
    sqrt(mean(((y_true - y_pred)/y_true)^2, na.rm = na.rm))

## mean absolute percent error
rmape <- function(y_true, y_pred, na.rm=TRUE)
    mean(abs((y_true- y_pred)/y_true), na.rm = na.rm) * 100



## from "Evaporated_Rainwater_Conc.txt"
ewat <- c("SOLUTION 1",
          "units mol/kgw",
          "temp 25",
          "pressure 1",
          "pH 3.749",
          "C(4)  1.147e-03",
          "Ca 1.504e-03",
          "Cl 1.321e-03",
          "K  1.559e-04",
          "Mg 1.081e-04",
          "Na 1.692e-04",
          "S(6)  1.216e-03",
          "PURE 1",
          "END")


## Biotite or "black mica": it's the Mg-Endmember of Mica, called
## Phlogopite. It is included as "Phlogopite" in the llnl.dat
## database. For the phreeqc.dat you can use this expression:
## --------------- cut here --------------------------
## Phlogopite
##     KMg3(AlSi3)O10(OH)2 + 10 H+ = 1 Al+3 + 1 K+ + 3 Mg+2 + 3 H4SiO4
##     log_k        41.082
##     -delta_H        -360.123    kJ/mol    # References :92cir/nav
##     -analytic    -1.7201279e+03    -2.6579576e-01    1.0718208e+05    6.1999929e+02    -4.7275095e+06
## #   Temperature range (K) :  273 - 1000
## #   References DGf or LogK n/a ; DHf 92cir/nav ; S 95rob/hem ; Cp 95rob/hem
## --------------- cut here --------------------------
## This definition is taken from the "thermoddem" database
## (https://thermoddem.brgm.fr/), primary literature is cited there.
## Muscovite or "white mica": it's the Al-endmember of mica, called
## "Muscovite" in llnl.dat and "K-mica" in phreeqc.dat, so you can
## just use these two minerals for each of the two databases.


## This is the standard "phreeqc.dat" database without unused stuff at
## the end and with added Phlogopite
phreeqc::phrLoadDatabase("phreeqc_invPra.dat")

## I sprlitted the groundwater models into 4 single ones: we read them
## scripts, we run them and import the results
fil <- list.files(pattern="Sol.*pqi")
inps <- lapply(fil, RPhreeFile)

phreeqc::phrSetOutputStringsOn(TRUE)
outl <- lapply(inps, function(x) {
    phreeqc::phrRunString(x)
    phreeqc::phrGetOutputStrings()
})

res <- lapply(outl, function(x) ReadOut(x)[[1]])

## Our target (N region) is the second element of the res list
target <- res[[2]]$tot
cvec <- target$molal
names(cvec) <- rownames(target)
## We also need pH
cvec <- c(cvec, pH=res[[2]]$desc["pH",])
## Remove the parentheses from names for consistency
names(cvec) <- sub("\\(.\\)","",names(cvec))


res[[1]]

## Vector of evaporated rainwater concentrations
init <- c(Al=NA, Ba=NA, C=1.091e-04, Ca=1.430e-04, Cl=1.321e-03,
          K=1.559e-04, Mg=1.081e-04, Na=1.692e-04, S=1.216e-03, Si=NA,
          pH=4.68 )

primary <- c("Albite", "Chlorite(14A)", "K-mica", "Phlogopite", "Calcite", "Anorthite", "K-feldspar",
             "Barite", "Pyrite", "Gibbsite")
secondary <- c("Kaolinite","Gypsum")


tot <- c(DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=3, procs=4),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=4, procs=5),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=5, procs=7),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=6, procs=7),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=7, procs=7),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=8, procs=7))

filtered <- Filter(tot, delta=0.5)




## Now we extract the same vector from all computed results and create
## a data.frame
dat <- sapply(names(cvec), function(x) RPinfo(filtered, ifelse(x=="pH", "desc", "tot"), x))



## Compute RMSLE between simulations and evaporated rainwater
res_rmse <- apply(dat, 1, function(x) rrmse(cvec,x))
res_mape <- apply(dat, 1, function(x) rmape(cvec,x))


## Count occurrences of minerals within the top 50

## Find out the best 50 simulations given a certain metric, in this
## case rrmse
inds_best <- order(res_rmse)
inds_best50 <- which(inds_best < 51)
## Extract all phase names from the best 50
AllPhases <- unlist(lapply(filtered[inds_best50], function(x) rownames(x$pphases)))
## Count
table(AllPhases)


## Best match: depends on the criterium!
which.min(res_rmse)
which.min(res_mape)

viz <- rbind(dat[which.min(res_rmse),], cvec)

out <- barplot(viz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10), col=c("orange", "light blue", "light green", "grey"))
legend("topleft", c("initial","rmsle/rmse", "mape", "target"), fill=c("orange", "light blue", "light green", "grey"), bty="n")



######### New visualization, specific components
comps <- c("Fe", "K")

## Look at argument "na.rm": if it is FALSE, and one of F and K is NA,
## then NA is returned
b1 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=comps, na.rm=FALSE)

## This is the previous behavior (na.rm=TRUE): just compute the RRMSE
## using the only components which are not NA!
b2 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=comps, na.rm=TRUE)

viz1 <- dat[order(b1),][1:5,]
viz2 <- dat[order(b2),][1:5,]

PlotComb <- function(res, samples, comp, ...) {
    if (missing(comp)) {
        comp <- intersect(colnames(res), colnames(samples))
    }
    tmp <- subset(res, select=comp)
    sam <- subset(samples, select=comp)

    mins <- apply(sam, 2, min, na.rm=TRUE)
    maxs <- apply(sam, 2, max, na.rm=TRUE)
    meds <- apply(sam, 2, median, na.rm=TRUE)
    meas <- apply(sam, 2, mean, na.rm=TRUE)

    colors <- heat.colors(nrow(tmp))
    
    out <- barplot(tmp, beside=TRUE, ylab="mol/kgw", log="y",
                   col=colors, las=1, ...)
    for (i in seq_along(comp)) {
        rect(out[1,i]-0.6, mins[i],out[nrow(out),i]+0.6, maxs[i],col=rgb(0,0,1.0,alpha=0.5))
        segments(out[1,i]-0.6, meds[i], out[nrow(out),i]+0.6, meds[i],col="red", lwd=2, lty="dashed")
        segments(out[1,i]-0.6, meas[i], out[nrow(out),i]+0.6, meas[i],col="grey", lwd=2, lty="dotted")
    }
}




library(readxl)
samples <- read_excel("Pra_data_M.xlsx", sheet=1)
samples[samples==0] <- NA

par(mfrow=c(2,1))
PlotComb(res=viz1, samples=samples)
PlotComb(res=viz2, samples=samples)

