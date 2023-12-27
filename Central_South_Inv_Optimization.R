# Time-stamp: "Last modified 2022-03-02 18:00:03 delucia"
# This is a combinatorial inverse modeling using modeled concentration as the initial solution and Northern groundwater as the target solution
library(RedModRphree)
library(magrittr)
source("rfun.R")
options(width = 116)
library(splitstackshape)
library("ggplot2")
library(readxl)

## This helper function allows to load other files from the same directory.

## source_here <- function(x, ...) {
##     dir <- "."
##     if(sys.nframe()>0) {
##         frame <- sys.frame(1)
##         if (!is.null(frame$ofile)) {
##             dir <- dirname(frame$ofile)
##         }
##     }
##     source(file.path(dir, x), ...)
## }
## source_here("rfun.R")


## from median Central waters
ewat <- c("SOLUTION 1",
          "units mol/kgw",
          "temp 25",
          "pressure 1",
          "pH 6.05",
          "Al 5.162e-08",
          "C(4)  6.802e-04",
          "Ca  1.798e-04",
          "Cl 1.139e-03",
          "Fe  2.314e-07",
          "K  9.677e-05",
          "Mg 1.010e-04",
          "Na 7.45e-04",
          "S(6)  1.638e-05",
          "Si 3.893e-04",
          "PURE 1",
          "END")


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
temp <- InputFromList(res[[13]])
## temp contains the input solution in phreeqc script for the starting solution

## Our target (N region) is the second elemental concentration of the res list
target <- res[[4]]$tot
cvec <- target$molal
names(cvec) <- rownames(target)
## We also need pH
cvec <- c(cvec, pH=res[[4]]$desc["pH",])
##cvec <- c(cvec, pH=res[[3]]$desc["pH",])
## Remove the parentheses from names for consistency
names(cvec) <- sub("\\(.\\)","",names(cvec))



## Newly computed concentration from evaporated rainwater to match chloride concentration of groundwater
## median solution concentrations
init <- c(Al=5.162e-08, C=6.802e-04, Ca=1.798e-04, Cl=1.134e-03, Fe=2.314e-07,
          K=9.677e-05, Mg=1.010e-04, Na=7.45e-04, S=1.638e-05, Si=3.893e-04,
          pH=6.05)

## Read the entire chemical data representing the Northern zone
##PraData <- read_excel("C:/Users/Asus/Desktop/Marco_combinv-main/CombInv/Pra_data_M.xlsx", sheet=1)
PraData <- read_excel("Pra_data_M.xlsx", sheet=1)


## modified phases
#primary <- c("Albite", "Phlogopite", "K-mica", "K-feldspar", "Quartz", 
#             "Chalcedony", "Pyrite")
#secondary <- c("Kaolinite", "Illite", 
#               "Fe(OH)3(a)", "Chlorite(14A)", "Calcite")

##Original Phases
#primary <- c("Albite", "Phlogopite", "K-mica", "K-feldspar", "Quartz", 
#             "Chalcedony", "Pyrite", "Anorthite", "Plagioclase")
#secondary <- c("Kaolinite", "Gibbsite", "Barite", "Illite", 
#              "Ca-Montmorillonite", "Goethite", "Hematite", "Chlorite(14A)", "Calcite")
#new phases used in the rainwater / southern zone inverse
primary <- c("Albite", "K-mica", "Phlogopite", "K-feldspar", "Quartz",
             "Chalcedony", "Fe(OH)3(a)")
secondary <- c("Kaolinite", "Illite", "Calcite", 
               "Goethite", "Chlorite(14A)")

##Tobias
#primary <- c("Albite", "Phlogopite", "K-mica", "K-feldspar", 
#            "Pyrite")
#secondary <- c("Kaolinite", "Chlorite(14A)", "Jarosite-K", "Fe(OH)3(a)",
#              "Calcite")

## we try to call all the 8 individual initial solutions in a phreeqc format and automate for use
#allsamp <- list(solA=c(".."),
#                 solB=c("..."),
#                solC=c("..."))
#totA <- vector(mode="list", length=lenght(allsamp))
#filteredA <- vector(mode="list", length=lenght(allsamp))
#datA<- vector(mode="list",length=lenght(allsamp))

#for (i in seq_along(allsamp)) {
#    totA[[i]] <- c(DoCombSim(allsamp[[i]], db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=3, procs=6),
#                   DoCombSim(initsol = allsamp[[i]], db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=4, procs=6),
#                   DoCombSim(initsol = allsamp[[i]], db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=5, procs=6),
#                   DoCombSim(initsol = allsamp[[i]], db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=6, procs=6),
#                   DoCombSim(initsol = allsamp[[i]], db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=7, procs=6),
#                   DoCombSim(initsol = allsamp[[i]], db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=8, procs=6))
#    filteredA[[i]] <- Filter(totA[[i]],delta=0.5)

#    datA[[i]] <- sapply(names(cvec), function(x) RPinfo(filteredA[[i]], ifelse(x=="pH", "desc", "tot"), x))
#}                             
################################################################################################################
## If we want to save the outputs then we can activate the commented codes below
#saveRDS(file="20220301_Northern_evolution_optimization_median.rds", object=list(dat=dat, tot=tot, primary=primary, secondary=secondary, filtered=filtered))
#saveRDS(file=paste0("20220301_Northern_evolution_optimization", i , ".rds", list(dat=datA[[i]], tot=totA[[i]], 
#                                                                                 primary=primary, secondary=secondary, filtered=filteredA[[i]])))


tot <- c(DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=3, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=4, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=5, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=6, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=7, procs=6),
         DoCombSim(initsol = ewat, db="phreeqc_invPra.dat", primary = primary, secondary = secondary, len=8, procs=6))

######################  Results
filtered <- Filter(tot, delta=0.5)

dat <- sapply(names(cvec), function(x) RPinfo(filtered, ifelse(x=="pH", "desc", "tot"), x))

################################################################################################################
## If we want to save the outputs then we can activate the commented codes below
saveRDS(file="20220301_Central_optimization_Final_Mean.rds", object=list(dat=dat, tot=tot, primary=primary, secondary=secondary, filtered=filtered))


## Read back the results
#tmp <- readRDS("20220301_Northern_evolution_optimization_SolA.rds")
#dat <- tmp$dat
#tot <- tmp$tot
#primary <- tmp$primary
#secondary <- tmp$secondary
#filtered <- tmp$filtered

################################################################################################################

## Find out the best 50 simulations given a certain metric, in this
## ## case rrmse
## ##inds_best <- order(res_rmse)
## ##inds_best50 <- which(inds_best < 51)
## ##Extract all phase names from the best 50
## ##AllPhases <- unlist(lapply(filtered[inds_best50], function(x) rownames(x$pphases)))
## ## Count
## ##table(AllPhases)


### Test

### Compute metric for all concentrations separately
## restot <- as.data.frame(lapply(colnames(dat), function(x) ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=x)),USE.NAMES=FALSE)
## colnames(restot)<- colnames(dat)

ind_all <- ComputeMetric(dat, target = cvec, FUN = "rrmse")

##Final best simulated concentration automated
#SimCop <- filtered[[which.min(ind_all)]]$tot$molal
#SimCop <- c(SimCop, pH=filtered[[which.min(ind_all)]]$desc["pH",], pe=filtered[[which.min(ind_all)]]$desc["pe",])

## mineral combinations of the best matched simulation
#MinCom <- filtered[[which.min(ind_all)]]$pphases$delta




########################################
##   2022 03 02: Marco's alternative  ##
########################################
ind <- which.min(ind_all) # ind produces the phreeqc calculated list of all parameters

## Inspect the data structure: "filtered" is a list of list 
filtered[[ind]]      ## this just prints it
str(filtered[[ind]]) ## this brings out the structure of it

## Use the function InputFromList() from the RedModRphree package: you
## have documentation!
?InputFromList


## filtered[[ind]] is the best matching combination we found until
## now. From these computed results, form a new PHREEQC input script:
NewInput <- InputFromList(filtered[[ind]])


## workhorse function to form a vector of all total elements + pH
ExtractComponents <- function(lin) {
  concs <- as.numeric(lin$tot[,1])
  names(concs) <- rownames(lin$tot)
  pH <- as.numeric(lin$desc["pH",1])
  final <- c(concs, pH=pH)
  return(final)
}

## workhorse function to run simulations. NOTE: I hard coded the
## database into it!
.runPQC <- function(input) {
  phreeqc::phrLoadDatabase("phreeqc_invPra.dat")
  phreeqc::phrSetOutputStringsOn(TRUE)
  phreeqc::phrRunString(input)
  tmpout <- phreeqc::phrGetOutputStrings()
  res <- RedModRphree::ReadOut(tmpout)[[1]]
  return(res)
}

## For calibration we need to define a function which takes as first
## argument the numerical value of the variable we want to calibrate
## for and the input script/template, and returns a metric (error)
## which we minimize.
#ToOptimizePE <- function(pe, input, target) {
#    tmp <- Distribute(input, prop="pe", values=pe)
## Now we use the initial evaporated solution as the input
ToOptimizePE <- function(pe, input, target) {
  tmp <- Distribute(input, prop="pe", values=pe)
  tmp <- Distribute(tmp, prop="pe_Fix", values=paste(-pe, " O2(g)  0.5")) 
  
  ## call PHREEQC and put the results into a list
  res <- .runPQC(tmp)
  
  ## we extract the total concentrations and pH and put them in a
  ## vector called "final"
  final <- ExtractComponents(res)
  
  ## make sure the "target" vector does not contain ()
  names(target) <- gsub("\\(.*$","", names(target))
  
  ## This function is written so that we can exclude some components
  ## from "target" by just removing them. However here we need to
  ## make sure that we compare the same components - compute the
  ## "intersection" of the vector names!
  components <- intersect(names(target), names(final))
  
  ## We compute the single numeric metric value. NB: I hard coded
  ## "rrmse" here but it can be changed
  ans <- rrmse(target[components], final[components], na.rm=FALSE)
  return(ans)
}

## We now use the R standard function "optimize()". Allowed range for
## pe is set to c(-8, 10)
#calib <- optimise(f=ToOptimizePE, input = NewInput, target = cvec,
#                  interval=c(-8,10), maximum = FALSE)## are the numbers subjective?
##Here we use the initial solution as the input
calib <- optimise(f=ToOptimizePE, input=inps[[13]], target = cvec,
                  interval=c(-7,7), maximum = FALSE)

## Look at the results
calib
str(calib)

temp <- Distribute(inps[[13]], prop="pe", values=calib$minimum) %>%
  Distribute(prop="pe_Fix", values=paste(-calib$minimum," O2(g)  0.5"))


## our final solution:
## require(magrittr) ## to use the %>% operator
## Calibrated <- Distribute(NewInput, prop="pe", values=calib$minimum) %>% .runPQC #this is originally commented out
#Calibrated <- .runPQC(Distribute(NewInput, prop="pe", values=calib$minimum))
Calibrated <- .runPQC(temp)
#InputFromList(Calibrated) this collects the parameter concentrations into a phreeqc script form

## comparison of the script before calibration and after

cviz <- rbind(ExtractComponents(Calibrated), ExtractComponents(filtered[[ind]]), target=cvec)
out <- barplot(cviz, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
               col=c("yellow", "light green", "grey"), las=1)
legend("topleft", c("Calibrated", "rrmse", "target"),
       fill=c("yellow", "light green", "grey"), bty="n")



## Now a full working Pyrite optimization achieved

#########################################################################
## another exemplary implementation varying Pyrite SI


## For calibration we need to define a function which takes as
## argument the numerical value of the variable we want to calibrate
## for and the input script/template, and returns a metric (error)
## which we minimize. 
ToOptimizePyriteSI <- function(sipy, input, target) {
  tmp <- Distribute(input, prop="Pyrite", values=sipy)
  tmp <- Distribute(tmp, prop="Pyrite", values=paste(-sipy, "2"))
  
  ## call PHREEQC and put the results into a list
  res <- .runPQC(tmp)
  
  ## we extract the total concentrations and pH and put them in a
  ## vector called "final"
  final <- ExtractComponents(res)
  
  ## make sure the "target" vector does not contain ()
  names(target) <- gsub("\\(.*$","", names(target))
  
  ## This function is written so that we can exclude some components
  ## from "target" by just removing them. However here we need to
  ## make sure that we compare the same components - compute the
  ## "intersection" of the vector names!
  components <- intersect(names(target), names(final))
  
  ## We compute the single numeric metric value. NB: I hardcoded
  ## "rrmse" here but it can be changed
  ans <- rrmse(target[components], final[components], na.rm=FALSE)
  return(ans)
}


## We now use the R standard function "optimise()". Allowed range for
## pe is set to c(-8, 10)
calibpyr <- optimise(f=ToOptimizePyriteSI, input = inps[[13]], target = cvec,
                     interval=c(-2,2), maximum = FALSE)

## Look at the results
calibpyr
str(calibpyr)

temp <- Distribute(inps[[13]], prop="Pyrite", values=calib$minimum) %>%
  Distribute(prop="Pyrite", values=paste(-calib$minimum,"2"))

## our final solution:
## require(magrittr) ## to use the %>% operator
## Calibrated <- Distribute(NewInput, prop="pe", values=calib$minimum) %>% .runPQC
#CalibPyr <- .runPQC(Distribute(inps[[13]], prop="Pyrite", values=paste(calibpyr$minimum, "2")))
CalibPyr <- .runPQC(temp)

cviz2 <- rbind(ExtractComponents(CalibPyr), ExtractComponents(filtered[[ind]]), cvec)
out <- barplot(cviz2, beside=TRUE, ylab="mol/kgw", log="y", ylim = c(1E-9,10),
               col=c("yellow", "light green", "grey"), las=1)


temp
Distribute(inps[[13]], prop="Pyrite", values=paste(calibpyr$minimum, "2"))




################## old file resumes here





## Count occurrences of minerals within the top 50
## Find out the best 50 simulations given a certain metric, in this
## case rrmse
##inds_best <- order(res_CaMg)
inds_best <- order(ind_all)[1:10]
## Extract all phase names from the best 50
AllPhases <- unlist(lapply(filtered[inds_best], function(x) rownames(x$pphases)))
#write.csv(AllPhases,"C:\\Users\\Asus\\Desktop\\Marco_combinv-main\\CombInv\\viz\\Frequency.csv", row.names = FALSE)
## Count
table<-table(AllPhases)
out <- barplot(table, ylab="Frequency", las=2, col=c("cyan"))


##viz <- rbind(init, dat[which.min(res_CaMg),], cvec)
viz <- rbind(init, dat[which.min(ind_all),], cvec)
## Plottingg the calibrated concentrations with initial and target solutions
#viz <- rbind(init, ExtractComponents(Calibrated), cvec)
## this writes the output of the viz table into csv and save it into the folder called viz
write.csv(viz, "/home/emanu/Documents/combinv/SolH.csv", row.names = FALSE)

##viz <- rbind(init, dat[which.min(res_Al),], cvec)
par(mfrow=c(2,2))
out <- barplot(viz, beside=TRUE, log="y", ylim = c(1E-9,10),
               col=c("orange", "light green", "grey"), cex.axis=1.2, cex=1)
mtext(side=2, line=3, "mol/kgw", font=1)

############################################################################################################
legend("topleft", c("initial","rrmse", "target"),
       fill=c("orange", "light green", "grey"), bty="n", cex=0.8)

## making a combinatorial plots for the best matched elemental solutions

##dev.new(width=10, height=6)
##textCa <- filtered[[which.min(res_CaMg)]]$pphases
##barplot(textCa$delta, names.arg=row.names(textCa), las=1, ylab="mol/kgw")
##legend("bottomleft", c("Ca"))

textCa <- filtered[[which.min(ind_all)]]$pphases
barplot(textCa$delta, names.arg=row.names(textCa), las=2, cex.axis=1.2, cex=1.3)
## this helps to offset the labels on the axis
mtext(side=2, line=4, "mol/kgw", font=1)
#mtext(side=2, line=4, "mol/kgw", col="black", font=2)
##legend("bottomleft", c("All"))

data_Plot <- filtered[[which.min(ind_all)]]$pphases$delta
write.csv(data_Plot, "/home/emanu/Documents/combinv/data_plot.csv", row.names = FALSE)

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

library(readxl)
samples <- read_excel("Pra_data_M.xlsx", sheet=1)
samples[samples==0] <- NA

par(mfrow=c(2,2))
PlotComb(res=viz1, samples=samples, cex.axis=1, cex=1)
#mtext(side=2, line=4, font=1)
#mtext(side=2, line=3, "mol/kgw", font=1)
##PlotComb(res=viz2, samples=samples)

