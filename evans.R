# Time-stamp: "Last modified 2021-12-13 12:28:06 delucia"
library(RedModRphree)
#source(rfun)
options(width = 116)
library(splitstackshape)
library("ggplot2")
library(readxl)

## This script contains the evaporated rainwater as initial solution and the Northern zone as the target solution

######################################################################################
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
        cat(":: Registered default doParallel cluster with ", procs, "nodes")
        parallel::clusterCall(cl=ThisRunCluster, phreeqc::phrLoadDatabase, db)
        msg(":: Database loaded on each worker")
        
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


#####################################################################################################################
##Code to implement a function that shows only selected components to be included in the final output
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
## root mean square logarithmic error
##rmsle <- function(y_true, y_pred)
##    sqrt(mean((log1p(y_true) - log1p(y_pred))^2, na.rm = TRUE))

## relative root mean square error
##rrmse <- function(y_true, y_pred)
##    sqrt(mean(((y_true - y_pred)/y_true)^2, na.rm = TRUE))

## relative mean absolute percent error
##rmape <- function(y_true, y_pred)
##    mean(abs((y_true- y_pred)/y_true), na.rm = TRUE) * 100

##Additional errors relative Mean absolute error (mae)
##rmae <- function(y_true, y_pred)
##    mean(abs((y_true- y_pred)/y_true), na.rm = TRUE)

rrmse <- function(y_true, y_pred, na.rm=TRUE)
    sqrt(mean(((y_true - y_pred)/y_true)^2, na.rm = na.rm))

## mean absolute percent error
rmape <- function(y_true, y_pred, na.rm=TRUE)
    mean(abs((y_true- y_pred)/y_true), na.rm = na.rm) * 100

##Additional errors relative Mean absolute error (mae)
rmae <- function(y_true, y_pred, na.rm=TRUE)
    mean(abs((y_true- y_pred)/y_true), na.rm = na.rm)



## from "Evaporated_Rainwater_Conc.txt which was achieved by 47.2% evaporation"
ewat <- c("SOLUTION 1",
          "units mol/kgw",
          "temp 25",
          "pressure 1",
          "pH 3.749",
          "C(4)  1.147e-03",
          "Ca  1.504e-03",
          "Cl 1.321e-03",
          "K  1.559e-04",
          "Mg 1.081e-04",
          "Na 1.692e-04",
          "S(6)  1.216e-03",
          "PURE 1",
          "END")
## Northern zone concentration
##ewat <- c("SOLUTION 1",
##         "units mol/kgw",
##         "temp 25",
##         "pressure 1",
##          "pH 3.749",
##          "C(4)  1.778e-03",
##          "Ca 4.248e-05",
##          "Cl 2.370e-04",
##          "K  1.791e-05",
##          "Mg 2.633e-04",
##          "Na 5.177e-04",
##          "S(6)  1.177e-05",
##          "Al 1.638e-07",
##          "Ba  1.457e-07",
##          "Si  3.978e-04",
##          "PURE 1",
##          "END")


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
## the end and with added Phlogopite and Plagioclase
phreeqc::phrLoadDatabase("phreeqc_invPra.dat")

## I splitted the groundwater models into 4 single ones: we read them
## scripts, we run them and import the results
fil <- list.files(pattern="Sol.*pqi")
inps <- lapply(fil, RPhreeFile)

phreeqc::phrSetOutputStringsOn(TRUE)
outl <- lapply(inps, function(x) {
    phreeqc::phrRunString(x)
    phreeqc::phrGetOutputStrings()
})

res <- lapply(outl, function(x) ReadOut(x)[[1]]) ## what is the meaning of 1 here

## Our target (N region) is the second elemental concentration of the res list
target <- res[[2]]$tot
cvec <- target$molal
names(cvec) <- rownames(target)
## We also need pH
cvec <- c(cvec, pH=res[[2]]$desc["pH",])
##cvec <- c(cvec, pH=res[[3]]$desc["pH",])
## Remove the parentheses from names for consistency
names(cvec) <- sub("\\(.\\)","",names(cvec))

## Vector of evaporated rainwater concentrations
##init <- c(Al=NA, C=2.065e-04, Ca=2.707e-04, Cl=2.378e-04,
##          K=2.807e-05, Mg=1.946e-05, Na=3.045e-05, S=2.188e-04, Si=NA, Fe=NA,
##          pH=4.431)
## Newly computed concentration from evaporated rainwater to match chloride concentration of groundwater

init <- c(Al=NA, C=1.147e-03, Ca=1.504e-03, Cl=1.321e-03, Fe=NA,
          K=1.559e-04, Mg=1.081e-04, Na=1.692e-04, S=1.216e-03, Si=NA,
          pH=3.749 )
## Read the entire chemical data representing the Northern zone
#PraData <- read_excel("C:/Users/Asus/Desktop/Marco_combinv-main/CombInv/Pra_data_M.xlsx", sheet=1)

## Vector of Northern zone concentrations
##init <- c(Al=1.638e-07, Ba=1.457e-07, C=1.778e-03, Ca=4.248e-05, Cl=2.370e-04,
##          K=1.791e-05, Mg=2.633e-04, Na=5.177e-04, S=1.177e-05, Si=3.978e-04,
##          pH=6.33 )


##primary <- c("Albite", "Chlorite(14A)", "K-mica", "Phlogopite", "Calcite", "K-feldspar")
## Original phases
##primary <- c("Albite", "Chlorite(14A)", "Phlogopite", "Calcite", "K-feldspar",
##            "Barite", "Pyrite", "Quartz", "K-mica")

## modified phases
primary <- c("Albite", "Phlogopite", "K-mica", "K-feldspar", "Quartz", 
             "Chalcedony", "Pyrite", "Anorthite", "Plagioclase")
secondary <- c("Kaolinite", "Gibbsite", "Barite", "Illite", 
               "Ca-Montmorillonite", "Goethite", "Hematite", "Chlorite(14A)", "Calcite")
## working phases
##secondary <- c("Kaolinite", "Gibbsite", "Barite", "Illite", 
##              "Ca-Montmorillonite", "Goethite", "Hematite", "Chlorite(14A)", "Gypsum", "Calcite", "Dolomite")

## A full list of possible secondary mineral phases
##secondary <- c("Kaolinite", "Gibbsite", "Barite", "Illite", "Hydroxyapatite", 
##"Ca-Montmorillonite", "Goethite", "Hematite", "Chlorite(14A)", "Gypsum", "Calcite", "Dolomite")

##DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, 
##                       "Ca-Montmorillonite" len=3, procs=4)

tot <- c(DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=3),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=4),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=5),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=7),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=8))

##tot <- c(DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=5),
##         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=6),
##         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=7),
##         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=8))

######################  Results
filtered <- Filter(tot, delta=0.5)



## Now we extract the same vector from all computed results and create
## a data.frame
dat <- sapply(names(cvec), function(x) RPinfo(filtered, ifelse(x=="pH", "desc", "tot"), x))

################################################################################################################
## If we want to save the outputs then we can activate the commented codes below
## saveRDS(file="20220111_MyFile.rds", object=list(dat=dat, tot=tot, primary=primary, secondary=secondary))
# Read back the results
##tmp <- readRDS("20220111_MyFile.rds")
##dat <- tmp$dat
##tot <- tmp$tot
##primary <- tmp$primary
##secondary <- tmp$secondary
################################################################################################################


## Compute RMSLE between simulations and evaporated rainwater
res_rmse <- apply(dat, 1, function(x) rrmse(cvec,x))
##res_rmsle <- apply(dat, 1, function(x) rmsle(cvec,x))
res_mape <- apply(dat, 1, function(x) rmape(cvec,x))
## Additional mean absolute error
res_mae <- apply(dat, 1, function(x) rmae(cvec,x))

## Count occurrences of minerals within the top 50

## Find out the best 50 simulations given a certain metric, in this
## case rrmse
##inds_best <- order(res_rmse)
##inds_best50 <- which(inds_best < 51)
##Extract all phase names from the best 50
##AllPhases <- unlist(lapply(filtered[inds_best50], function(x) rownames(x$pphases)))
## Count
##table(AllPhases)

## Best match: depends on the criterion!
##which.min(res_rmsle)
which.min(res_rmse)
which.min(res_mape)
## Additional
which.min(res_mae)


### Test
## res_all <- ComputeMetric(dat, target = cvec, FUN = "rrmse")
## all.equal(res_rmse, res_all)
colnames(dat)
restot <- as.data.frame(lapply(colnames(dat), function(x) ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=x)),USE.NAMES=FALSE)
colnames(restot)<- colnames(dat)
##comps<- c("Ca")
##comps<- c("Mg")
##comps<- c("Na")
comps<- c("K")
##comps<- c("Al")
##comps<- c("Si")
##comps<- c("S")
##comps<- c("C")
##comps<- c("Ba")
##comps<- c("Fe", "K")
##comps<- c("Fe")
##comps<- c("pH")
##comps <- c("Ca", "Mg", "Na")
##comps <- c("Na", "K")
##comps <- c("Na", "Ca", "Al")
##comps <- c("Ca", "Mg", "Na")

res_CaMg <- ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=comps)
##res_Al <- ComputeMetric(dat, target = cvec, FUN = "rrmse", comp="Na")
##res_1 <- ComputeMetric(dat, target = cvec, FUN = "rrmse")

## Count occurrences of minerals within the top 50
## Find out the best 50 simulations given a certain metric, in this
## case rrmse
inds_best <- order(res_CaMg)
inds_best50 <- which(inds_best < 51)
## Extract all phase names from the best 50
AllPhases <- unlist(lapply(filtered[inds_best50], function(x) rownames(x$pphases)))
## Count
table(AllPhases)


viz <- rbind(init, dat[which.min(res_CaMg),], cvec)
##viz <- rbind(init, dat[which.min(res_1),], cvec)

##viz <- rbind(init, dat[which.min(res_Al),], cvec)
par(mfrow=c(1,1))
out <- barplot(viz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
               col=c("orange", "light green", "grey"))
## how to plot error bars
##standard_error<- 6.08
##out <- ggplot(viz, aes(restot, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
##              ) + geom_bar(stat = "identity") + 
##               geom_errorbar(aes(ymin=value-standard_error, ymax=value+standard_error), width=.2)
############################################################################################################
legend("topleft", c("initial","rmse", "target"),
       fill=c("orange", "light green", "grey"), bty="n")

## making a combinatorial plots for the best matched elemental solutions

##dev.new(width=10, height=6)
textCa <- filtered[[which.min(res_CaMg)]]$pphases
barplot(textCa$delta, names.arg=row.names(textCa), las=1, ylab="mol/kgw")
legend("bottomleft", c("Ca"))

textCa <- filtered[[which.min(res_1)]]$pphases
barplot(textCa$delta, names.arg=row.names(textCa), las=1, ylab="mol/kgw")
##legend("bottomleft", c("All"))

textMg <- filtered[[which.min(restot$Mg)]]$pphases
barplot(textMg$delta, names.arg=row.names(textMg), las=2, ylab="mol/kgw")
legend("topleft", c("Mg"))

textNa <- filtered[[which.min(restot$Na)]]$pphases
barplot(textNa$delta, names.arg=row.names(textNa), las=2, ylab="mol/kgw")
legend("topleft", c("Na"))

textK <- filtered[[which.min(restot$K)]]$pphases
barplot(textK$delta, names.arg=row.names(textK), las=2, ylab="mol/kgw")
legend("topleft", c("K"))

textAl <- filtered[[which.min(restot$Al)]]$pphases
barplot(textAl$delta, names.arg=row.names(textAl), las=2, ylab="mol/kgw")
legend("topleft", c("Al"))

textSi <- filtered[[which.min(restot$Si)]]$pphases
barplot(textSi$delta, names.arg=row.names(textSi), las=2, ylab="mol/kgw")
legend("topleft", c("Si"))
###############################################################################################################
##par(mfrow=c(3,2))
##textC <- tot[[5888]]$pphases
##barplot(textC$delta, names.arg=row.names(textC), ylab="mol/kgw")
##legend("topleft", c("C"))
##textBa <- tot[[12495]]$pphases
##barplot(textBa$delta, names.arg=row.names(textBa), ylab="mol/kgw")
##legend("topleft", c("Ba"))
##textS <- tot[[13617]]$pphases
##barplot(textS$delta, names.arg=row.names(textS), ylab="mol/kgw")
##legend("topleft", c("S"))
##textCl <- tot[[34811]]$pphases
##barplot(textCl$delta, names.arg=row.names(textCl), ylab="mol/kgw")
##legend("topleft", c("Cl"))
##textpH <- tot[[87]]$pphases
##barplot(textpH$delta, names.arg=row.names(textpH), ylab="mol/kgw")
##legend("topleft", c("pH"))


##viz <- rbind(init, dat[which.min(res_rmse), ], cvec) ##Original code
##viz <- rbind(init, dat[which.min(res_rmsle)[21]], dat[order(res_mape)[21]], cvec)
##viz <- rbind(init, dat[order(res_rmsle)[4],], cvec)
##viz <- rbind(init, dat[order(res_mape)[242],], cvec)
##viz <- cSplit(init, dat[which.min(res_rmse),], cvec, "|")

## Original codes
##out <- barplot(viz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10), col=c("orange", "light blue", "grey"))
##legend("topleft", c("initial","rmse/mape/mae", "target"), fill=c("orange", "light blue", "grey"), bty="n")

## mOdified code
##out <- barplot(viz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10), col=c("orange", "light blue", "grey"))
##legend("topleft", c("initial", "mape", "target"), fill=c("orange", "light blue", "grey"), bty="n")


#############################################################
## some codes to remember
## which.min(abs(dat[,"Cl"]-target["Cl",1])) to define the best matching results between observed and predicted
## order(res_rmsle) and order(res_mape) to list the best matching results in order of least to highest
## length(tot) to define the number of simulations
## which.min(abs(dat[,"Mg"]-target["Mg",1]))

## Codes to comit and push to the server
## git commit -a -m "Added the additional code to evans.R"
## git push ---helps to update the serve with changes made in the local directory
## git pull ---helps to retrieve updated files from the server
## cp CombInv.R evans.R ---to copy and paste a file into another file
## git add evans.R ---to add the newly created file to the git server
## git remote -v ---to check if the directory is a git repository
## textMg <- tot[[55784]]$pphases
##plot(c(1:8),textMg$delta)
## barplot(c(1:8),textMg$delta)
## barplot(textMg$delta, names.arg=row.names(textMg))
## ?barplot
## order(res_CaMg) ## This orders the best matched simulations in order of best to worst
## filtered[[which.min(res_CaMg)]]$tot to show the simulated conc for the overall best solution
## filtered[[which.min(res_CaMg)]]$tot["Ca",1] show the best simulated conc for Ca

###################################################################################################################P
PraData[PraData==0] <- NA
apply(PraData[,-c(1:3)],2,median,na.rm=TRUE)

####################################################################################################################
######### New visualization, specific components
##comps <- c("Fe", "K")
comps <- c("Fe")

######### New visualization modified

## Look at argument "na.rm": if it is FALSE, and one of F and K is NA,
## then NA is returned
b1 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", na.rm=FALSE) ##comp=comps,

## This is the previous behavior (na.rm=TRUE): just compute the RRMSE
## using the only components which are not NA!
b2 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", na.rm=TRUE) ##comp=comps,

viz1 <- dat[order(b1),][1:5,]
viz2 <- dat[order(b2),][1:5,]
##
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

par(mfrow=c(2,2))
PlotComb(res=viz1, samples=samples)
##PlotComb(res=viz2, samples=samples)

#########################################################################################################################################








######### New visualization

bm <- ComputeMetric(dat, target = cvec, FUN = "rrmse")

viz <- dat[order(bm),][1:5,]

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

## Plotting
samples <- read_excel("Pra_data_M.xlsx", sheet=3)
samples[samples==0] <- NA
PlotComb(res=viz, samples=samples)
##library(readxl)

apply(samples[,-c(1:3)],2,median,na.rm=TRUE)


## Simple regression  model with r code

#library(ehaGoF)

#input <- 0:4

# Target vector, observed values, dependent variable
#target <- c(1.9, 4.1, 5.89, 7.9, 10.01)

# Simple linear regression, target across input like: target = a * input + b,
# where a and b are coefficients.
#model <- lm(target~input)

# Information about the model
#summary(model)

# Values predicted by the model
#predicted <- predict(model)

# using library ehaGoF for goodness of fit
#library(ehaGoF)

# Goodness of fit - relative root mean square error (RRMSE)
#gofRRMSE(target, predicted)
