# Time-stamp: "Last modified 2022-03-02 18:00:03 delucia"
# This is a combinatorial inverse modeling using modeled concentration as the initial solution and Northern groundwater as the target solution
library(RedModRphree)
library(magrittr)
source("rfun.R")
options(width = 116)
library(splitstackshape)
library("ggplot2")
library(readxl)


## from median concentration of all samples in the northern zone
ewat <- c("SOLUTION 1",
          "units mol/kgw",
          "temp 25",
          "pressure 1",
          "pH 6.04",
          "Al     2.076e-07",
          "C(4)  1.082e-03",
          "Ca  3.197e-04",
          "Cl 2.384e-04",
          "Fe   1.683e-07",
          "K  2.123e-05",
          "Mg 2.143e-04",
          "Na 5.059e-04",
          "S(6)  2.270e-05",
          "Si   2.750e-04",
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
temp <- InputFromList(res[[1]])
## temp contains the input solution in phreeqc script for the starting solution

## Our target (N region) is the second elemental concentration of the res list
target <- res[[6]]$tot
cvec <- target$molal
names(cvec) <- rownames(target)
## We also need pH
cvec <- c(cvec, pH=res[[6]]$desc["pH",])
##cvec <- c(cvec, pH=res[[3]]$desc["pH",])
## Remove the parentheses from names for consistency
names(cvec) <- sub("\\(.\\)","",names(cvec))


## Precipitation solution evaporated match chloride concentration of groundwater
init <- c(C=1.082e-03, Ca=3.197e-04, Cl=2.384e-04, Fe=1.683e-07,
          K=2.123e-05, Mg=2.143e-04, Na=5.059e-04, S=2.270e-05, Si=2.750e-04,
          pH=6.040)

## Read the entire chemical data representing the Northern zone
##PraData <- read_excel("C:/Users/Asus/Desktop/Marco_combinv-main/CombInv/Pra_data_M.xlsx", sheet=1)
PraData <- read_excel("Pra_data_M.xlsx", sheet=2)

primary <- c("Albite", "Plagioclase", "Anorthite", "Phlogopite", "K-mica", "K-feldspar", 
             "Fe(OH)3(a)", "Calcite")
secondary <- c("Kaolinite", "Ca-Montmorillonite", "Chlorite(14A)", "Quartz", "Chalcedony", "Pyrite")

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
#saveRDS(file="20220301_Northern_evolution_optimization_SolH.rds", object=list(dat=dat, tot=tot, primary=primary, secondary=secondary, filtered=filtered))


################################################################################################################

ind_all <- ComputeMetric(dat, target = cvec, FUN = "rrmse")

colnames(dat)
restot <- as.data.frame(lapply(colnames(dat), function(x) ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=x)),USE.NAMES=FALSE)
colnames(restot)<- colnames(dat)


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

## Count occurrences of minerals within the top 50
## Find out the best 50 simulations given a certain metric, in this
## case rrmse
inds_best <- order(ind_all)[1:50]
# The corresponding rrmse values are:
best50_rrmse <- ind_all[order(ind_all)[1:50]]
## Extract all phase names from the best 50
AllPhases <- unlist(lapply(filtered[inds_best], function(x) rownames(x$pphases)))
#write.csv(AllPhases,"C:\\Users\\Asus\\Desktop\\Marco_combinv-main\\CombInv\\viz\\Frequency.csv", row.names = FALSE)
## Count
table<-table(AllPhases)
out <- barplot(table, ylab="Frequency", las=2, col=c("cyan"))

table<-table(AllPhases)
par(mar=c(12.4,5,1,1))
out <- barplot(table, ylab="Frequency",cex.axis=1.6, cex=1.6, cex.lab=1.6, las=2, col=c("cyan"))


##viz <- rbind(init, dat[which.min(res_CaMg),], cvec)
viz <- rbind(init, dat[which.min(ind_all),], cvec)

par(mfrow=c(1,1))
out <- barplot(viz, beside=TRUE, log="y", ylim = c(1E-9,10),
               col=c("orange", "light green", "grey"), cex.axis=1.2, cex=1)
mtext(side=2, line=3, "mol/kgw", font=1)

legend("topleft", c("initial","rrmse", "target"),
       fill=c("orange", "light green", "grey"), bty="n", cex=0.8)


textCa <- filtered[[which.min(ind_all)]]$pphases
barplot(textCa$delta, names.arg=row.names(textCa), las=2, cex.axis=1.2, cex=1.3)
## this helps to offset the labels on the axis
mtext(side=2, line=4, "mol/kgw", font=1)

######### New visualization modified

## Look at argument "na.rm": if it is FALSE, and one of F and K is NA,
## then NA is returned
b1 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", na.rm=FALSE) ##comp=comps,
#b1 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", comp=comps, na.rm=FALSE) ##comp=comps,

## This is the previous behavior (na.rm=TRUE): just compute the RRMSE
## using the only components which are not NA!
b2 <- ComputeMetric(dat, target = cvec, FUN = "rrmse", na.rm=TRUE) ##comp=comps,

viz1 <- dat[order(b1),][1:50,]
viz2 <- dat[order(b2),][1:50,]

library(readxl)
samples <- read_excel("Pra_data_M_Good.xlsx", sheet=2)
samples[samples==0] <- NA

par(mar=c(2.5,6,1,1))
PlotComb(res=viz1, samples=samples, cex.axis=1.7, cex=1.7, cex.lab=1.4)
mtext(side=2, line=4.5, "mol/kgw", cex=2)


### MDL 2023 02 03 Filter

## Extract the ranges from the "samples" tibble excluding the first 2
## columns
ranges <- sapply(samples[, 3:12], range)
## dat is a matrix holding the simulation results

IsInRange <- function(matsim, ranges) {
  ## We check concentration-wise (column-wise) if we are within the
  ## range. This return a matrix of the same dimension as "matsim"
  ## filled with TRUE, FALSE or NA
  sapply(colnames(ranges), function(conc) ifelse(ranges[1, conc]< matsim[,conc] & ranges[2,conc] > matsim[,conc], TRUE, FALSE))
}

inrange <- IsInRange(matsim=dat, ranges=ranges)

## How to use: TRUE equals 1 and FALSE 0, so we want all the rows
## whose sum is equal the number of columns in "ranges":
which(rowSums(inrange)==ncol(ranges)) ## 0: NO SIMULATION FALLS WITHIN THE RANGE

## Repeat excluding Fe and K (resp. column 4 and 5 in "ranges")
inrange_nokfe <- IsInRange(matsim=dat, ranges=ranges[, -c(4,5)])

## How to use: TRUE equals 1 and FALSE 0:
which(rowSums(inrange_nokfe)==ncol(ranges[, -c(4,5)])) ## 229

