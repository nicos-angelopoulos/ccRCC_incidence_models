library(ggplot2)
library(reshape2)
library(MASS)
library(lme4)
library(pander)
library(survival)
library(MCMCpack)
library(extraDistr)
library(abind)

# Set-up some basic script parameters
if (!exists("rcc.iters")) { rcc.iters <- 50000 }        # number of iterations for the Monte Carlo Markov chain (MCMC)
if (!exists("rcc.burn.in")) { rcc.burn.in <- 20000 }      # number of burn ins for the MCMC
if (!exists("chrono.seed"))  { chrono.seed <- 31 }

set.seed(chrono.seed)

print(rcc.iters)
print(rcc.burn.in)
print(chrono.seed)

# Read in patient data
pt.df <- read.table("data/summary/patient_summary.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
pt.df <- pt.df[pt.df$include == "yes",]

# Read in phylogenetic trees
branch.df <- data.frame(Sample=character(0), Name=character(0), Num.clonal=integer(0), Num.mutations=integer(0))

branch.files <- system("ls params/Tree_structures/*optima*", intern = TRUE)
branch.files <- branch.files[substring(branch.files, first = 24, last = 30) %in% pt.df$sanger.ID]
branch.files <- branch.files[!(substring(branch.files, first = 24, last = 30) %in% c("PD24238", "PD24242"))] 
# This previous line removes a patient for which only one sample gave good sequence; and a tumour with two completely independent clones (ie two separate cancers)

for (i in branch.files)
{
  temp <- suppressWarnings(read.table(i, header=TRUE, sep="\t", stringsAsFactors = FALSE))
  sample <- substring(i, 24, 30)
  lri.sample <- pt.df[pt.df$"sanger.ID"==sample,2]
  
  temp.term <- data.frame(Sample = rep(sample, sum(temp$Terminal == "Yes")), 
                          LRISample = rep(lri.sample, sum(temp$Terminal == "Yes")), 
                          Name = temp$Name[temp$Terminal == "Yes"],
                          Num.clonal = rep(temp$no.of.mutations.assigned[temp$Name=="T"], sum(temp$Terminal=="Yes")),
                          Num.mutations = rep(0, sum(temp$Terminal == "Yes")),
                          Age = rep(pt.df$Age[pt.df$sanger.ID == sample], sum(temp$Terminal == "Yes")))
  
  # This next loop calculates the number of mutations from root of phylogenetic tree to tip of each terminal branch
  for (j in temp.term$Name)
  {
    constituents <- unlist(strsplit(j, ">"))
    for (k in 2:length(constituents)) {constituents[k] <- paste(constituents[k-1], constituents[k], sep=">")}
    temp.term$Num.mutations[temp.term$Name == j] <- sum(temp$no.of.mutations.assigned[temp$Name %in% constituents])
  }
  branch.df <- rbind(branch.df, temp.term)  
}

# Now plot the relationship between number of mutations on each branch and age
ggplot(data = branch.df, aes(x=Age, y=Num.mutations, color=LRISample)) + 
  geom_point(shape=16) +
  xlim(c(0,80)) + ylim(c(0,10000)) +
  xlab("Age (years)") + ylab("Number of mutations per subclone")

# Test fits with and without intercepts

muts.per.year.lmer.with.pt.intercepts <- lmer(Num.mutations ~ Age + (Age | Sample), data=branch.df, REML=FALSE)
muts.per.year.lmer.with.pop.intercept <- lmer(Num.mutations ~ Age -1 + (Age | Sample), data=branch.df, REML=FALSE)
muts.per.year.lmer <- lmer(Num.mutations ~ Age - 1 + (Age - 1 | Sample), data=branch.df, REML=FALSE)

kable(anova(muts.per.year.lmer, muts.per.year.lmer.with.pop.intercept, muts.per.year.lmer.with.pt.intercepts))

# Print summary of model without intercepts
print(summary(muts.per.year.lmer))

 # Generate the bootstrap dataset for this LME
 mySumm.muts.per.yr <- function(.) {
   unlist(fixef(.) + ranef(.)$Sample[,1]) 
 }
 
 boot2 <- bootMer(muts.per.year.lmer, mySumm.muts.per.yr, nsim=1000, use.u=TRUE, type="parametric")
 colnames(boot2$t) <- row.names(ranef(muts.per.year.lmer)$Sample)
 
branch.df$Scaled.num.muts <- scale(branch.df$Num.mutations, center=FALSE)
scale.factor <- attr(branch.df$Scaled.num.muts, which = "scaled:scale")

age.lmer <- lmer(Age ~ Scaled.num.muts - 1 + (Scaled.num.muts - 1 | Sample), data=branch.df)
branch.df$Pt.fitted.num.muts <- scale.factor / (ranef(age.lmer)$Sample[branch.df$Sample,] + fixef(age.lmer)["Scaled.num.muts"]) * branch.df$Age

# Visual check of model's fit
ggplot(data = branch.df, aes(x=Age, y=Num.mutations, color=LRISample)) + 
  geom_segment(aes(x=0, y = 0, xend=Age, yend=Pt.fitted.num.muts, colour=LRISample), data=branch.df) +
  geom_abline(intercept=0,slope=scale.factor / fixef(age.lmer), size=3) +
  geom_point(shape=16) +
  xlim(c(0,80)) + ylim(c(0,10000))

# Generate point estimates for the timing of the MRCA emergence
new.dat <- data.frame(Sample=branch.df$Sample, LRISample=branch.df$LRISample, Scaled.num.muts=branch.df$Num.clonal / scale.factor, Age=branch.df$Age)
new.dat <- new.dat[!duplicated(new.dat$Sample),]
new.dat$MRCA.pred.age <- predict(age.lmer, newdata = new.dat)

# Now generate 95% CIs on the predictions for MRCA timing
# Note that CIs for lme models are challenging - bootstrapping appears to be most robust

mySumm <- function(.) {
  predict(., newdata=new.dat, re.form=NULL)
}

####Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

boot1 <- bootMer(age.lmer, mySumm, nsim=1000, use.u=FALSE, type="parametric")

PI.boot1 <- sumBoot(boot1)
PI.boot1$Age <- new.dat$Age
PI.boot1$pred.fit <- new.dat$MRCA.pred.age
PI.boot1$Time.lag <- PI.boot1$Age - PI.boot1$fit
PI.boot1$Time.lag[PI.boot1$Time.lag < 0] <- PI.boot1$Age[PI.boot1$Time.lag < 0] - PI.boot1$pred.fit[PI.boot1$Time.lag < 0]
PI.boot1$lag.lwr <- PI.boot1$Age - PI.boot1$upr
PI.boot1$lag.lwr[PI.boot1$lag.lwr < 0] <- 0
PI.boot1$lag.upr <- PI.boot1$Age - PI.boot1$lwr
PI.boot1$ID <- new.dat$Sample
PI.boot1$LRIID <- new.dat$LRISample
PI.boot1 <- PI.boot1[order(PI.boot1$Time.lag),]

plot(PI.boot1$Time.lag, 1:nrow(PI.boot1), pch=20, xlim=c(-10,max(PI.boot1$lag.upr)), axes=FALSE, xlab="Estimated time between MRCA and diagnosis (years)", ylab="")
segments(x0 = PI.boot1$lag.lwr, x1 = PI.boot1$lag.upr, y0 = 1:nrow(PI.boot1))
axis(side=1, at=(0:5)*10, labels=(0:5)*10)
text(x = PI.boot1$lag.lwr-1, y = 1:nrow(PI.boot1), labels = PI.boot1$LRIID, adj = 1, cex=0.75)

# First we define a function that reads in the relevant mutations and decides whether a variant is clonal or subclonal; and if clonal whether it occurs pre- or post-duplication

sa.df <- read.table("data/summary/sample_summary.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
sa.df <- sa.df[sa.df$include == "yes",]

timing <- function(pt = "PD24236", sample.suffixes = c("a", "c", "d", "e"), normal.contam, chr, lower, upper) 
{
  cn <- list()

  # Read in the mutation data for that patient
  file.name <- system(paste("ls data/patient_data/merged/", pt, "/1152_*_consolidated_snp.tsv", sep=""), intern=TRUE)
  header.line.num <- system(paste("grep -n VariantID data/patient_data/merged/", pt, "/1152_*_consolidated_snp.tsv", sep=""), intern=TRUE)
  header.line.num <- as.double(strsplit(header.line.num[2], split = ":")[[1]][1])
  muts <- read.table(file.name, header=TRUE, blank.lines.skip = TRUE, stringsAsFactors = FALSE, sep="\t", comment.char="#", row.names=NULL, skip=header.line.num-1)
  muts <- muts[muts$Chrom == chr & muts$Pos >= lower & muts$Pos <= upper,]
  
  # For each sample in that patient, calculate VAF and get CN at each mutated base-pair
  for (j in 1:length(sample.suffixes)) 
  {
    k <- sample.suffixes[j]
    muts[,paste(k,"VAF",sep="_")] <- muts[,paste(pt, k, "_MTR", sep="")] / muts[,paste(pt, k, "_DEP", sep="")]

    # Get CN file for that patient
    cn.file <- system(paste("ls data/sample_data/BB_output/", pt, k, "*subclones.txt", sep=""), intern=TRUE)
    cn <- read.table(cn.file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    
    # Get total and minor copy number estimates for each base in the mutations file where the region shows no or minimal (<20%) subclonal CN variation
    muts[, paste(k,"T.Minor.CN",sep="_")]   <- sapply(1:nrow(muts), function(i,x,y) {
            temp.CN <- y$nMin1_A[y$chr == x$Chrom[i] & y$startpos <= x$Pos[i] & y$endpos >= x$Pos[i] & y$frac1_A > 0.8];
            if (length(temp.CN)>0) {return(temp.CN[1])} else {return(NA)}}, 
            x=muts, y=cn)
    muts <- muts[!is.na(muts[,paste(k,"T.Minor.CN",sep="_")]),]
    muts[, paste(k,"T.Total.CN",sep="_")] <- muts[,paste(k,"T.Minor.CN",sep="_")] + sapply(1:nrow(muts), function(i,x,y) {
            temp.CN <- y$nMaj1_A[y$chr == x$Chrom[i] & y$startpos <= x$Pos[i] & y$endpos >= x$Pos[i] & y$frac1_A > 0.8];
            if (length(temp.CN)>0) {return(temp.CN[1])} else {return(NA)}}, 
            x=muts, y=cn)
    
    # Calculate the CCF (fraction of cancer cells carrying mutation) from VAF and sample-specific normal cell contamination fraction     
    muts[,paste(k,"CCF",sep="_")] <- muts[,paste(k,"VAF",sep="_")] * ((1-normal.contam[j]) * muts[,paste(k,"T.Total.CN",sep="_")] + normal.contam[j] * 2) / (1-normal.contam[j])

    # From CCF, infer whether mutation occurred pre- or post-duplication in that sample
    muts[,paste(k,"Status",sep="_")] <- rep("Uninformative", nrow(muts))
    muts[muts[,paste(k,"T.Total.CN",sep="_")] >= 2 & muts[,paste(k,"T.Minor.CN",sep="_")] == 0 & muts[,paste(k,"CCF",sep="_")] > 1.5, paste(k,"Status",sep="_")] <- "Preduplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] >= 2 & muts[,paste(k,"T.Minor.CN",sep="_")] == 0 & muts[,paste(k,"CCF",sep="_")] <= 1.5, paste(k,"Status",sep="_")] <- "Post-duplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] == 3 & muts[,paste(k,"T.Minor.CN",sep="_")] == 1 & muts[,paste(k,"CCF",sep="_")] > 1.5, paste(k,"Status",sep="_")] <- "Preduplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] == 3 & muts[,paste(k,"T.Minor.CN",sep="_")] == 1 & muts[,paste(k,"CCF",sep="_")] <= 1.5, paste(k,"Status",sep="_")] <- "Post-duplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] == 6 & muts[,paste(k,"T.Minor.CN",sep="_")] == 2 & muts[,paste(k,"CCF",sep="_")] > 2.5, paste(k,"Status",sep="_")] <- "Preduplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] == 6 & muts[,paste(k,"T.Minor.CN",sep="_")] == 2 & muts[,paste(k,"CCF",sep="_")] <= 2.5, paste(k,"Status",sep="_")] <- "Post-duplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] == 5 & muts[,paste(k,"T.Minor.CN",sep="_")] == 2 & muts[,paste(k,"CCF",sep="_")] > 2.5, paste(k,"Status",sep="_")] <- "Preduplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] == 5 & muts[,paste(k,"T.Minor.CN",sep="_")] == 2 & muts[,paste(k,"CCF",sep="_")] <= 2.5, paste(k,"Status",sep="_")] <- "Post-duplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] >= 4 & muts[,paste(k,"T.Minor.CN",sep="_")] == 1 & muts[,paste(k,"CCF",sep="_")] > 1.5, paste(k,"Status",sep="_")] <- "Preduplication"
    muts[muts[,paste(k,"T.Total.CN",sep="_")] >= 4 & muts[,paste(k,"T.Minor.CN",sep="_")] == 1 & muts[,paste(k,"CCF",sep="_")] <= 1.5, paste(k,"Status",sep="_")] <- "Post-duplication"
    
  }

  # Now, from all statuses across different samples, get a consensus status
  # First, read in the branch assignments of all mutations
  tree.file <- system(paste("ls params/Tree_structures/", pt, "*optima*", sep=""), intern = TRUE)
  tree.struct <- read.table(file=tree.file, header=TRUE, stringsAsFactors = FALSE, sep="\t")
  branch.assign <- read.table(paste("data/patient_data/clusters/", pt, "/", pt, "_assigned_mutations.txt", sep=""), header=TRUE, sep="\t", stringsAsFactors = FALSE)
  branch.assign <- branch.assign[branch.assign$chr == chr & branch.assign$pos >= lower & branch.assign$pos <= upper,]
  
  muts <- muts[muts$Pos %in% branch.assign$pos,]
  muts$Branch <- sapply(muts$Pos, function(i,x) {x$cluster[x$pos == i]}, x=branch.assign)

  # Identify clonal mutations (assigned to trunk of phylogenetic tree & >1 read reporting variant in all samples)
  muts$Timing <- rep("Subclonal", nrow(muts))
  muts$Timing[muts$Branch %in% (1:nrow(tree.struct))[tree.struct$Name == "T"] &
              sapply(1:nrow(muts), function(i,x,pt.name,samps) {all(x[i,paste(pt.name, samps, "_MTR", sep="")] > 1)}, x=muts, pt.name=pt, samps=sample.suffixes)] <- "Uncertain"
  muts$Timing[sapply(1:nrow(muts), function(i,x) {sum(x[i,grepl("_Status", names(x))] == "Preduplication") > length(sample.suffixes)/2 & x$Timing[i] == "Uncertain"}, x=muts)] <- "Clonal Preduplication"
  muts$Timing[sapply(1:nrow(muts), function(i,x) {sum(x[i,grepl("_Status", names(x))] == "Post-duplication") > length(sample.suffixes)/2 & x$Timing[i] == "Uncertain"}, x=muts)] <- "Clonal Post-duplication"
  
  # Generate plot of CCF per sample
  for (i in sample.suffixes)
  {
    point.col <- as.character(factor(muts$Timing, levels=c("Clonal Post-duplication", "Clonal Preduplication", "Subclonal", "Uncertain"), labels=c("steelblue2", "springgreen3", "plum3", "black")))
    point.col[muts[,paste(pt, i, "_MTR", sep="")] > 1 & muts$Timing == "Subclonal"] <- "orange"
    plot(muts$Pos, muts[,paste(i, "CCF", sep="_")], pch=20, col=point.col, xlab=paste("Chr", chr, "position"), ylab="Number of mutation copies per cancer cell", las=1, ylim=c(0,max(muts[,paste(i, "CCF", sep="_")])+1))
    legend("topright", pch=20, col=c("steelblue2", "springgreen3", "orange", "plum3", "black"), legend=c("Clonal, post-duplication", "Clonal, pre-duplication", "Subclonal, present in this sample", "Subclonal, absent in this sample", "Uncertain"))
    sanger.sid <- paste(pt, i, sep="")
    lri.sid <- sa.df[sa.df$"sanger.ID"==sanger.sid,2]
    title(main = lri.sid)
  }
  
  return(muts)
}

  # The first method of calculating timing uses MRCA estimates and fraction of mutations at ploidy 1 and ploidy 2 to estimate length 
 triploid.2.1.with.mrca.est <- function(ploidy.1, ploidy.2, bs, ID, age) 
 {
   # bs is bootstrap output from LME
   point.est <- ploidy.2 / (ploidy.2 + (ploidy.1 - ploidy.2)/3)
   iter <- bs$R
   total.muts <- ploidy.1 + ploidy.2
   prob.2 <- ploidy.2 / total.muts
     
   resamps <- rbinom(n=iter, size=total.muts, p=prob.2)
   ests.bs <- resamps / (resamps + (total.muts - 2*resamps)/3)
   ests.bs[ests.bs > 1] <- 1
   
   age.mrca.bs <- bs$t[,ID]
   age.mrca.bs[age.mrca.bs > age] <- age
   
   ests.bs <- ests.bs * age.mrca.bs
   
   return(quantile(ests.bs, c(0.5,0.025,0.975)))
 }
 
 
 # The second method estimates the timing directly from the estimated mutation rate and number of mutations acquired before duplication
  triploid.direct.est <- function(ploidy.2, bs, ID, dup.region.size.in.Gb, callable.genome.size=5.32) 
 {
   # bs is bootstrap output from above
   iter <- bs$R
   mut.rate.bs <- bs$t[,ID] / callable.genome.size
   resamps <- rpois(n=iter, lambda = ploidy.2)
   ests.bs <- resamps / dup.region.size.in.Gb / mut.rate.bs 

   return(quantile(ests.bs, c(0.5,0.025,0.975)))
 }

 
 # Set up the data-frame for the results
 tloc.ids = c("PD21872", "PD21874", "PD21877", "PD23583", "PD23585", "PD24231", "PD24232", "PD24234", "PD24235", "PD24236", "PD24237", "PD24239")
 tloc.ids.len <- length(tloc.ids)

 lri.ids <- rep("", tloc.ids.len)
 for (x in 1:tloc.ids.len) {
  lri.ids[x] <- pt.df[pt.df$"sanger.ID"==tloc.ids[x],2]
 }

 t3_5.abs.time <- data.frame(ID = tloc.ids, triploid.lwr = rep(0, tloc.ids.len), triploid.est = rep(0, tloc.ids.len), triploid.upr = rep(0, tloc.ids.len), stringsAsFactors = FALSE)
 row.names(PI.boot1) <- PI.boot1$ID
 t3_5.abs.time$Age <- PI.boot1[t3_5.abs.time$ID, "Age"]
 t3_5.abs.time$MRCA.est <- PI.boot1[t3_5.abs.time$ID, "Age"] - PI.boot1[t3_5.abs.time$ID, "Time.lag"]
 colnames(boot1$t) <- new.dat$Sample
 t3_5.abs.time$Est.mut.rate.per.Gb.per.year <- (fixef(muts.per.year.lmer) + ranef(muts.per.year.lmer)$Sample[t3_5.abs.time$ID,"Age"]) / 6
 t3_5.abs.time$alt.trip.upr <- t3_5.abs.time$alt.trip.est <- t3_5.abs.time$alt.trip.lwr <- rep(0,tloc.ids.len)
 t3_5.abs.time$LRIID <- lri.ids
 row.names(t3_5.abs.time) <- t3_5.abs.time$ID
 
 # Run the patients through the functions in turn ([params/timing_inputs.csv](file:///params/timing_inputs.csv));
time.inputs <- read.table( "params/timing_inputs.csv", header=F, sep=",", stringsAsFactors=F, strip.white=T, as.is=T )

for (i in 1:nrow(time.inputs)) {
   branch.names <- unlist(strsplit(time.inputs[i,2],";"))
   branch.parms <- as.numeric(unlist(strsplit(time.inputs[i,3],";")))
   branches.pid <- time.inputs[i,1]     # patient ID
   branches.age <- time.inputs[i,8]     # patient age
   timing(branches.pid, branch.names, branch.parms, time.inputs[i,4], time.inputs[i,5], time.inputs[i,6]) -> a
   t3_5.abs.time[t3_5.abs.time$ID == branches.pid, c("triploid.est", "triploid.lwr", "triploid.upr")] <- 
   triploid.2.1.with.mrca.est(ploidy.1 = table(a$Timing)["Clonal Post-duplication"], 
                              ploidy.2 = table(a$Timing)["Clonal Preduplication"], 
                              bs=boot1, ID=branches.pid, age=branches.age)
   t3_5.abs.time[branches.pid, c("alt.trip.est", "alt.trip.lwr", "alt.trip.upr")] <- triploid.direct.est(table(a$Timing)["Clonal Preduplication"], boot2, branches.pid, time.inputs[i,7]/1000)
}

# Report the results
 # kable(t3_5.abs.time,row.names=FALSE)
 t3_5.abs.time.rep <- t3_5.abs.time
 row.names(t3_5.abs.time.rep) <- t3_5.abs.time$LRIID
 kable(t3_5.abs.time.rep)

 # Generate plot of estimated age of occurrence of t(3;5) - method 1
 t3_5.abs.time <- t3_5.abs.time[order(t3_5.abs.time$triploid.est),]
 plot(t3_5.abs.time$triploid.est, 1:nrow(t3_5.abs.time), pch=20, xlim=c(0,max(t3_5.abs.time$triploid.upr)), axes=FALSE, xlab="Estimated age at t(3;5) (years)", ylab="")
 segments(x0 = t3_5.abs.time$triploid.lwr, x1 = t3_5.abs.time$triploid.upr, y0 = 1:nrow(t3_5.abs.time))
 axis(side=1, at=(0:7)*5, labels=(0:7)*5)
 text.x.poss <- t3_5.abs.time$triploid.lwr-0.5
 text.x.poss[1] <- text.x.poss[1] + 0.5
 text(x = text.x.poss, y = 1:nrow(t3_5.abs.time), labels = t3_5.abs.time$LRIID, adj = 1, cex=0.75)
 title("Estimated age of t(3;5) - method 1")
 
 
  # Generate plot of estimated age of occurrence of t(3;5) - method 2
 plot(t3_5.abs.time$alt.trip.est, 1:nrow(t3_5.abs.time), pch=20, xlim=c(0,max(t3_5.abs.time$alt.trip.upr)), axes=FALSE, xlab="Estimated age at t(3;5) (years)", ylab="")
 segments(x0 = t3_5.abs.time$alt.trip.lwr, x1 = t3_5.abs.time$alt.trip.upr, y0 = 1:nrow(t3_5.abs.time))
 axis(side=1, at=(0:7)*5, labels=(0:7)*5)
 text.x.poss <- t3_5.abs.time$alt.trip.lwr-0.5
 text.x.poss[1] <- text.x.poss[1] + 0.5
 text(x = text.x.poss, y = 1:nrow(t3_5.abs.time), labels = t3_5.abs.time$LRIID, adj = 1, cex=0.75)
 title("Estimated age of t(3;5) - method 2")

rev.comp <- function(x) {
  x <- toupper(x)
  x <- rev(unlist(strsplit(x,"")))
  x <- paste(sapply(x, switch,  "A"="T", "T"="A","G"="C","C"="G"),collapse="")
  return(x)
}

# Define the per year per genome substitution rate
rate.lmer <- lmer(Num.mutations ~ Age - 1 + (Age - 1 | Sample), data=branch.df)
SUB.RATE <- fixef(rate.lmer)["Age"]
CALLABLE.GENOME.SIZE <- 5.32E9

# Read in the data for the gene and mutation signatures
# gene.cds has to include 1 non-coding base before CDS and one after for mutation signature analysis
gene.cds <- substring(read.table(file="../../data/sampler_data/VHL_transcript.txt", header=FALSE, stringsAsFactors = FALSE, fill=TRUE)$V1[1], 213, 856)
cosmic.muts <- read.csv("../../data/sampler_data/COSMIC_VHL_systematic_somatic.csv", quote="", stringsAsFactors = FALSE)

# Read in and manipulate the file for random positions in the genome
rand.pos <- read.table("../../data/sampler_data/Random_build37_positions_with_context.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rand.pos <- rand.pos[!grepl("N", rand.pos$CONTEXT, fixed=TRUE),]
rand.pos$Pyrimidine.WT <- rand.pos$WT
rand.pos$Pyrimidine.WT[rand.pos$WT == "A"] <- "T"
rand.pos$Pyrimidine.WT[rand.pos$WT == "G"] <- "C"
rand.pos$Pyrimidine.context <- substring(rand.pos$CONTEXT, 10, 12)
rand.pos$Pyrimidine.context[rand.pos$WT %in% c("A", "G")] <- sapply(rand.pos$Pyrimidine.context[rand.pos$WT %in% c("A", "G")], rev.comp)

# Now, calculate the per trinucleotide per mutation type mutation rate
observed.mut.spectrum <- read.table("../../data/sampler_data/t3_5_pts_chr5_muts.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
observed.mut.spectrum$Pyrimidine.WT <- observed.mut.spectrum$WT
observed.mut.spectrum$Pyrimidine.WT[observed.mut.spectrum$WT == "A"] <- "T"
observed.mut.spectrum$Pyrimidine.WT[observed.mut.spectrum$WT == "G"] <- "C"
observed.mut.spectrum$Pyrimidine.Alt <- observed.mut.spectrum$Alt
observed.mut.spectrum$Pyrimidine.Alt[observed.mut.spectrum$WT %in% c("A", "G")] <- sapply(observed.mut.spectrum$Pyrimidine.Alt[observed.mut.spectrum$WT %in% c("A", "G")], rev.comp)
observed.mut.spectrum$Pyrimidine.context <- substring(observed.mut.spectrum$CONTEXT, 10, 12)
observed.mut.spectrum$Pyrimidine.context[observed.mut.spectrum$WT %in% c("A", "G")] <- sapply(observed.mut.spectrum$Pyrimidine.context[observed.mut.spectrum$WT %in% c("A", "G")], rev.comp)
observed.mut.spectrum$Mut.type <- paste(observed.mut.spectrum$Pyrimidine.WT, ">", observed.mut.spectrum$Pyrimidine.Alt, sep="")
obs.mut.table <- table(observed.mut.spectrum$Mut.type, observed.mut.spectrum$Pyrimidine.context) / nrow(observed.mut.spectrum)

trinucs.in.genome <- table(rand.pos$Pyrimidine.context) / nrow(rand.pos)
mut.rate.per.type.per.trinuc <- data.frame(Mut.type = rep(sort(unique(observed.mut.spectrum$Mut.type)), each=16), 
                                           Trinuc = c(rep(paste(rep(c("A", "C", "G", "T"), each=4), "C", rep(c("A", "C", "G", "T"), 4), sep=""), 3), rep(paste(rep(c("A", "C", "G", "T"), each=4), "T", rep(c("A", "C", "G", "T"), 4), sep=""), 3)),
                                           stringsAsFactors = FALSE)
Prob.per.mut.type.per.trinuc <- sapply(1:nrow(mut.rate.per.type.per.trinuc), function(i,x,y,z) {y[x$Mut.type[i], x$Trinuc[i]] / (z[x$Trinuc[i]] * CALLABLE.GENOME.SIZE * 2)}, x=mut.rate.per.type.per.trinuc, y=obs.mut.table, z=trinucs.in.genome) # Need to double callable genome size because diploid and looking for VHL mutation on specific parental haplotype
names(Prob.per.mut.type.per.trinuc) <- paste(mut.rate.per.type.per.trinuc$Mut.type, mut.rate.per.type.per.trinuc$Trinuc, sep="::")

# Read in the genetic code
genetic.code.df <- read.table("../../data/sampler_data/Genetic_code.txt", sep=" ", stringsAsFactors = FALSE, col.names=c("codon", "aa", "aa.code"))
genetic.code <- genetic.code.df$aa.code
names(genetic.code) <- genetic.code.df$codon
genetic.code <- substring(genetic.code, 2, 2)

# Generate all possible mutations
poss.drivers <- data.frame(Ref = rep(unlist(strsplit(substring(gene.cds,2,nchar(gene.cds)-1), split = "")), each=4), 
                           Alt = rep(c("A", "C", "G", "T"), nchar(gene.cds)-2), 
                           cDNA.pos = rep(1:(nchar(gene.cds)-2), each=4),
                           aa.pos = rep(1:((nchar(gene.cds)-2)/3), each=12),
                           Pos.in.codon = rep(rep(1:3, each=4), (nchar(gene.cds)-2)/3),
                           Base.5prime = rep(unlist(strsplit(substring(gene.cds,1,nchar(gene.cds)-2), split = "")), each=4),
                           Base.3prime = rep(unlist(strsplit(substring(gene.cds,3,nchar(gene.cds)), split = "")), each=4),
                           Codon = rep(sapply(seq(2,nchar(gene.cds)-3,3), function(i,x) substring(x, i, i+2), x=gene.cds), each=12), stringsAsFactors = FALSE)

poss.drivers$Alt.codon <- sapply(1:nrow(poss.drivers), function(i,x) {cod <- unlist(strsplit(x$Codon[i],"")); cod[x$Pos.in.codon[i]] <- x$Alt[i]; return(paste(cod, collapse=""))}, x=poss.drivers)
poss.drivers <- poss.drivers[genetic.code[poss.drivers$Codon] != genetic.code[poss.drivers$Alt.codon],]
poss.drivers$Mut.string <- paste("p.", genetic.code[poss.drivers$Codon], poss.drivers$aa.pos, genetic.code[poss.drivers$Alt.codon], sep="")

poss.drivers$Pyrimidine.WT <- poss.drivers$Ref
poss.drivers$Pyrimidine.WT[poss.drivers$Ref == "A"] <- "T"
poss.drivers$Pyrimidine.WT[poss.drivers$Ref == "G"] <- "C"
poss.drivers$Pyrimidine.Alt <- poss.drivers$Alt
poss.drivers$Pyrimidine.Alt[poss.drivers$Ref %in% c("A", "G")] <- sapply(poss.drivers$Pyrimidine.Alt[poss.drivers$Ref %in% c("A", "G")], rev.comp)
poss.drivers$Pyrimidine.context <- paste(poss.drivers$Base.5prime, poss.drivers$Ref, poss.drivers$Base.3prime, sep="")
poss.drivers$Pyrimidine.context[poss.drivers$Ref %in% c("A", "G")] <- sapply(poss.drivers$Pyrimidine.context[poss.drivers$Ref %in% c("A", "G")], rev.comp)
poss.drivers$Mut.type <- paste(poss.drivers$Pyrimidine.WT, ">", poss.drivers$Pyrimidine.Alt, sep="")

poss.drivers <- poss.drivers[poss.drivers$aa.pos == 1 |
  (genetic.code[poss.drivers$Alt.codon] == "X" & genetic.code[poss.drivers$Codon] != "X") |
  (genetic.code[poss.drivers$Alt.codon] != "X" & genetic.code[poss.drivers$Codon] == "X") |
  poss.drivers$Mut.string %in% cosmic.muts$AA.Mutation,]

# Calculate estimated VHL driver substitution rate per year
vhl.driver.sub.rate <- sum(Prob.per.mut.type.per.trinuc[paste(poss.drivers$Mut.type, poss.drivers$Pyrimidine.context, sep="::")]) * SUB.RATE

print(paste("Average driver substitution rate in VHL:", signif(vhl.driver.sub.rate, digits = 2), "/cell/year"))

# Now we need to consider indel drivers
# The indel rate is ~5% the SUB.RATE, and we can consider all indels occurring within VHL CDS as drivers
vhl.driver.indel.rate <- (1/8.5) * SUB.RATE / CALLABLE.GENOME.SIZE * 644
print(paste("Average driver indel rate in VHL:", signif(vhl.driver.indel.rate, digits = 2), "/cell/year"))

# Therefore
vhl.driver.mutation.rate <- vhl.driver.indel.rate + vhl.driver.sub.rate
print(paste("Average driver point mutation rate in VHL:", signif(vhl.driver.mutation.rate, digits = 2), "/cell/year"))

NUM.VHL.CASES <- 100
NUM.SPORADIC <- 1000

# Read in age-incidence curves
vhl.penetrance <- read.table("../../data/sampler_data/Penetrance of RCC in vHL disease.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
cruk.sporadic.incidence <- read.table("../../data/sampler_data/Age-incidence of sporadic RCC.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

# Generate random sample of NUM.VHL.CASES age-incidences and censoring values for inherited RCCs
# Note that the vhl.penetrance data are cumulative Kaplan-Meier incidences and occur in nearly yearly increments
vhl.penetrance <- rbind(vhl.penetrance, c(63,1))
vhl.penetrance$lower <- c(0, vhl.penetrance$Penetrance[1:(nrow(vhl.penetrance)-1)])

vhl.incidence <- data.frame(Type = rep("vHL", NUM.VHL.CASES), stringsAsFactors = FALSE)
draws <- runif(n = NUM.VHL.CASES, min=0, max=1)
vhl.incidence$Age <- sapply(draws, function(i,x) x$Age[x$lower <= i & x$Penetrance > i], x=vhl.penetrance) 
vhl.incidence$Censored <- as.double(vhl.penetrance$lower[nrow(vhl.penetrance)] <= draws)

# Generate random sample of NUM.CASES age-incidences and censoring values for sporadic RCCs
# Note that the cruk.sporadic.incidence data are age-incidence rates per year per 100,000 population in 5 year blocks
cruk.sporadic.incidence$Rate.per.year <- (cruk.sporadic.incidence$Male.Rates + cruk.sporadic.incidence$Female.Rates)/2/100000
cruk.by.year <- data.frame(Age = 1:95, 
                           Incidence = rep(cruk.sporadic.incidence$Rate.per.year, each=5))
cruk.by.year$Penetrance <- 1 - cumprod(1-cruk.by.year$Incidence)
cruk.by.year <- rbind(cruk.by.year, c(95, 1, 1))
cruk.by.year$lower <- c(0, cruk.by.year$Penetrance[1:(nrow(cruk.by.year)-1)])


sporadic.incidence <- data.frame(Type = rep("Sporadic", NUM.SPORADIC), stringsAsFactors = FALSE)
draws <- runif(n = NUM.SPORADIC, min=0, max=1)
sporadic.incidence$Age <- sapply(draws, function(i,x) x$Age[x$lower <= i & x$Penetrance > i], x=cruk.by.year) 
sporadic.incidence$Censored <- as.double(cruk.by.year$lower[nrow(cruk.by.year)] <= draws)

# Generate MLE estimators for gamma distribution of age of 3p loss
ages.3p.loss <- t3_5.abs.time$triploid.est
s <- log(mean(ages.3p.loss)) - mean(log(ages.3p.loss))
alpha.1.est <- (3 - s + sqrt((s-3)^2 + 24*s)) / (12*s)
beta.1.est <- alpha.1.est / mean(ages.3p.loss)

alpha.1.est
beta.1.est

# Set-up the Gibbs sampler
rcc.age.gs <- function(spor.ages, vhl.ages, chr3p.ages, d.sf, g.sf, iterations, alpha.1.init = alpha.1.est, beta.1.init = beta.1.est, lambda.init=1e-3) {
  # Gibbs sampler for the three waiting times model of RCC development
  # spor.ages is a data.frame of ages of incidence / censoring of sporadic RCC
  # vhl.ages is a data.frame of ages of incidence / censoring of RCC in vHL disease
  # chr3p.ages is a vector of ages of incidence of chr3p loss inferred from t(3;-) cases
  # d.sf is the scaling factor for the Dirichlet proposal distribution - larger values mean smaller jumps
  # g.sf is the scaling factor for the Gamma proposal distribution - larger values mean smaller jumps
  # iterations is the number of iterations of the Gibbs sampler

  # Split age-incidences into censored and non-censored
  spor.ages.cens.T <- spor.ages[spor.ages$Censored == 1,]
  spor.ages.cens.F <- spor.ages[spor.ages$Censored == 0,]
  vhl.ages.cens.T <- vhl.ages[vhl.ages$Censored == 1,]
  vhl.ages.cens.F <- vhl.ages[vhl.ages$Censored == 0,]
  
  age.of.censoring.vhl <- unique(vhl.ages.cens.T$Age)
  age.of.censoring.sporadic <- unique(spor.ages.cens.T$Age)
  if (length(c(age.of.censoring.vhl,age.of.censoring.sporadic)) != 2) {stop("Age at censoring should be unique")} 
  
  # Set up initial values and hyperparameters
  alpha.1 <- beta.1 <- alpha.3 <- beta.3 <- accepted.1 <- accepted.3 <- lambda <- rep(NA, iterations)
  alpha.1[1] <- alpha.1.init
  beta.1[1] <- beta.1.init
  
  alpha.3[1] <- 20
  beta.3[1] <- 0.5
  
  accepted.1[1] <- accepted.3[1] <- TRUE # Acceptance rates for the M-H algorithm for the Z1 & Z3 gamma parameters
  
  lambda[1] <- lambda.init # The VHL driver mutation rate per year
  
  # Values for hyperparameters on gamma distribution for Z1
  log.hyperparameter.p.1 <- sum(log(chr3p.ages))
  hyperparameter.q.1 <- sum(chr3p.ages)
  hyperparameter.r.1 <- hyperparameter.s.1 <- length(chr3p.ages)
  
  # Set up matrices of waiting times
  vhl.wait.times.cens.F <- array(data = NA, dim=c(nrow(vhl.ages.cens.F), 2, iterations))
  vhl.wait.times.cens.F[,1,1] <- rep(15, nrow(vhl.ages.cens.F))
  vhl.wait.times.cens.F[,2,1] <- vhl.ages.cens.F$Age - 15

  vhl.wait.times.cens.T <- array(data = NA, dim=c(nrow(vhl.ages.cens.T), 2, iterations))
  vhl.wait.times.cens.T[,1,1] <- rep(15, nrow(vhl.ages.cens.T))
  vhl.wait.times.cens.T[,2,1] <- vhl.ages.cens.T$Age - 5

  spor.wait.times.cens.F <- array(data = NA, dim=c(nrow(spor.ages.cens.F), 3, iterations))
  spor.wait.times.cens.F[,1,1] <- rep(15, nrow(spor.ages.cens.F))
  spor.wait.times.cens.F[,2,1] <- (spor.ages.cens.F$Age - 15) / 2
  spor.wait.times.cens.F[,3,1] <- (spor.ages.cens.F$Age - 15) / 2
  
  spor.wait.times.cens.T <- array(data = NA, dim=c(nrow(spor.ages.cens.T), 3, iterations))
  spor.wait.times.cens.T[,1,1] <- rep(15, nrow(spor.ages.cens.T))
  spor.wait.times.cens.T[,2,1] <- (spor.ages.cens.T$Age - 5) / 2
  spor.wait.times.cens.T[,3,1] <- (spor.ages.cens.T$Age - 5) / 2
  
  sampled.wait.times.less0 <- sampled.wait.times.less0.25 <- sampled.wait.times.less0.5 <- sampled.wait.times.less0.75 <- matrix(NA, nrow=nrow(spor.ages), ncol=iterations)
  
  # Run the Gibbs sampler
  for (iter in 2:iterations) {

    ########################################################## 
    # Update of Z1 and Z3 for VHL cases that were not censored
    # Done by Metropolis-Hastings with proposal distribution as Dirichlet centred on curr Z1 & Z3 and scaled
    # Choice of Dirichlet enforces constraint that waiting times must sum to the age of incidence (with scaling)
    vhl.wait.times.cens.F[,,iter] <- censored.wait.times.GS(vhl.wait.times.cens.F[,,iter-1], d.sf, c(alpha.1[iter-1], alpha.3[iter-1]), c(beta.1[iter-1], beta.3[iter-1]))
    
    ######################################################
    # Update of Z1 and Z3 for vHL cases that were censored
    # Done with rejection sampling
    counter <- 1
    while(counter <= nrow(vhl.ages.cens.T)) {
      z1.proposed <- rgamma(1, alpha.1[iter-1], beta.1[iter-1])
      z3.proposed <- rgamma(1, alpha.3[iter-1], beta.3[iter-1])
      age.proposed <- z1.proposed + z3.proposed
      
      if (age.proposed > age.of.censoring.vhl) {
        vhl.wait.times.cens.T[counter,1,iter] <- z1.proposed
        vhl.wait.times.cens.T[counter,2,iter] <- z3.proposed
        counter <- counter+1
      } 
    }
    
    ###################################################################
    # Update of Z1, Z2 and Z3 for sporadic cases that were not censored
    # Done by Metropolis-Hastings with proposal distribution as Dirichlet centred on curr Z1, Z2 & Z3 and scaled
    # Choice of Dirichlet enforces constraint that waiting times must sum to the age of incidence (with scaling)
    spor.wait.times.cens.F[,,iter] <- censored.wait.times.GS(spor.wait.times.cens.F[,,iter-1], d.sf, c(alpha.1[iter-1], 1, alpha.3[iter-1]), c(beta.1[iter-1], lambda[iter-1], beta.3[iter-1]))
    
    ##########################################################
    # Update of Z1, Z2 and Z3 for vHL cases that were censored
    # Done with rejection sampling
    counter <- 1
    while(counter <= nrow(spor.ages.cens.T)) {
      z1.proposed <- rgamma(1, alpha.1[iter-1], beta.1[iter-1])
      z2.proposed <- rexp(1, lambda[iter-1])
      z3.proposed <- rgamma(1, alpha.3[iter-1], beta.3[iter-1])
      age.proposed <- z1.proposed + z2.proposed + z3.proposed
      
      if (age.proposed > age.of.censoring.sporadic) {
        spor.wait.times.cens.T[counter,1,iter] <- z1.proposed
        spor.wait.times.cens.T[counter,2,iter] <- z2.proposed
        spor.wait.times.cens.T[counter,3,iter] <- z3.proposed
        counter <- counter+1
      } 
    }
       
    ##############################################
    # Update age-incidence distribution parameters
    update.1 <- gamma.parameter.GS(alpha.1[iter-1], beta.1[iter-1], log.hyperparameter.p.1, hyperparameter.q.1, hyperparameter.r.1, hyperparameter.s.1, c(vhl.wait.times.cens.T[,1,iter], vhl.wait.times.cens.F[,1,iter]), g.sf)
    alpha.1[iter] <- update.1[1]
    beta.1[iter] <- update.1[2]
    accepted.1[iter] <- update.1[3]
    
    update.3 <- gamma.parameter.GS(alpha.3[iter-1], beta.3[iter-1], 0,1,1,1, c(vhl.wait.times.cens.T[,2,iter], vhl.wait.times.cens.F[,2,iter]), g.sf)
    alpha.3[iter] <- update.3[1]
    beta.3[iter] <- update.3[2]
    accepted.3[iter] <- update.3[3]

    # Update rate of VHL driver mutations contributing to Z2
    # Taken from conjugate prior to the exponential (ie gamma) with uninformative priors
    lambda[iter] <- rgamma(1, shape = 0.01 + length(c(spor.wait.times.cens.F[,2,iter], spor.wait.times.cens.T[,2,iter])), rate = 0.01 + sum(c(spor.wait.times.cens.F[,2,iter], spor.wait.times.cens.T[,2,iter])))
    
    # Sample alternate waiting times for clones of different sizes
    sampled.wait.times.less0[,iter] <- sample.alternate.wait.times(spor.wait.times.cens.F[,,iter], spor.wait.times.cens.T[,,iter], lambda[iter], 1)
    sampled.wait.times.less0.25[,iter] <- sample.alternate.wait.times(spor.wait.times.cens.F[,,iter], spor.wait.times.cens.T[,,iter], lambda[iter], 0.75)
    sampled.wait.times.less0.5[,iter] <- sample.alternate.wait.times(spor.wait.times.cens.F[,,iter], spor.wait.times.cens.T[,,iter], lambda[iter], 0.5)
    sampled.wait.times.less0.75[,iter] <- sample.alternate.wait.times(spor.wait.times.cens.F[,,iter], spor.wait.times.cens.T[,,iter], lambda[iter], 0.25)
    
    # print(lambda[iter])
  }
  
  return(list(alpha.1=alpha.1, beta.1=beta.1, alpha.3=alpha.3, beta.3=beta.3, accepted.1=accepted.1, accepted.3=accepted.3, vhl.wait.times.cens.F=vhl.wait.times.cens.F, vhl.wait.times.cens.T=vhl.wait.times.cens.T, spor.wait.times.cens.T=spor.wait.times.cens.T, spor.wait.times.cens.F=spor.wait.times.cens.F, lambda=lambda, sampled.wait.times.less0=sampled.wait.times.less0, sampled.wait.times.less0.25=sampled.wait.times.less0.25, sampled.wait.times.less0.5=sampled.wait.times.less0.5, sampled.wait.times.less0.75=sampled.wait.times.less0.75))
}

######################################################################
censored.wait.times.GS <- function(curr.WT, d.sf, alphas, betas) {
  # Function to provide updates of estimated waiting times for censored cases (ie observed ages of incidence)
  # Done by Metropolis-Hastings with proposal distribution as Dirichlet centred on curr Z1, Z2 & Z3 and scaled
  # Choice of Dirichlet enforces constraint that waiting times must sum to the age of incidence (with scaling)
  # curr.WT is matrix of the current estimates of waiting times
  # d.sf is the scaling factor for the Dirichlet proposal distribution
  # alphas is a vector of the current alpha values for the gamma distributions of waiting times
  # betas is a vector of the current beta values for the gamma distributions of waiting times
  
  # Generate proposed ages
  proposal <- t(sapply(1:nrow(curr.WT), 
                            function(i,x) rdirichlet(1, d.sf * x[i,] / sum(x[i,])),
                            x=curr.WT))
  obs.ages <- rowSums(curr.WT)
  proposal.scaled <- proposal * obs.ages
  
  # Calculate importance ratios
  log.proposal.dist.ratio <- ddirichlet(curr.WT / obs.ages, d.sf * proposal, log=TRUE) - 
      ddirichlet(proposal, d.sf * curr.WT / obs.ages, log=TRUE)
  
  log.density.ratio <- rowSums( t((alphas-1) * t(log(proposal.scaled) - log(curr.WT)) - betas * t(proposal.scaled - curr.WT)) )
  
  importance.ratio <- exp(log.proposal.dist.ratio + log.density.ratio)
  
  if (!all(!is.na(importance.ratio))) {
    importance.ratio[is.na(importance.ratio)] <- 0
  }
  
  draws <- runif(n=nrow(curr.WT))
  curr.WT[draws < importance.ratio,] <- proposal.scaled[draws < importance.ratio,] 
  return(curr.WT)
}

######################################################################
gamma.parameter.GS <- function(curr.alpha, curr.beta, log.hyper.p, hyper.q, hyper.r, hyper.s, x, gamma.sf) {
  # Function to return updates on gamma parameters from its conjugate prior
  # curr.alpha and curr.beta are the current values of the alpha and beta parameters
  # log.hyper.p, hyper.q, hyper.r and hyper.s are hyperparameters for the conjugate prior for curr.alpha & curr.beta
  # x is a vector of the data values
  # gamma.sf is the scaling factor for the M-H algorithm

  # Up-date alpha & beta from the conjugate prior 
  # Due to the unnormalised nature of the conjugate, we sample using Metropolis-Hastings since the normalisation factor will cancel out
  # Use independent gamma proposal distributions for alpha and beta
  alpha.proposed <- rgamma(n=1, curr.alpha * gamma.sf, gamma.sf)
  beta.proposed <- rgamma(n=1, curr.beta * gamma.sf, gamma.sf)
  
  proposal.dist.ratio <- dgamma(x = curr.alpha, shape = alpha.proposed * gamma.sf, rate = gamma.sf) * 
    dgamma(x = curr.beta, shape = beta.proposed * gamma.sf, rate = gamma.sf) /
    dgamma(x = alpha.proposed, shape = curr.alpha * gamma.sf, rate = gamma.sf) /
    dgamma(x = beta.proposed, shape = curr.beta * gamma.sf, rate = gamma.sf)
  
  log.post.p <- log.hyper.p + sum(log(x))
  post.q <- hyper.q + sum(x)
  post.r <- hyper.r + length(x)
  post.s <- hyper.s + length(x)
   
  log.density.ratio <- (alpha.proposed - curr.alpha) * log.post.p - 
    post.q * (beta.proposed - curr.beta) - 
    post.r * (lgamma(alpha.proposed) - lgamma(curr.alpha)) +
    post.s * (alpha.proposed * log(beta.proposed) - curr.alpha * log(curr.beta))
  
  importance.ratio <- proposal.dist.ratio * exp(log.density.ratio)
  if (runif(n=1) < importance.ratio) {
    new.alpha <- alpha.proposed
    new.beta <- beta.proposed
    accepted <- TRUE
  } else {
    new.alpha <- curr.alpha
    new.beta <- curr.beta
    accepted <- FALSE
  }
  
  return(c(new.alpha, new.beta, accepted))
}

######################################################################
sample.alternate.wait.times <- function(WT.cens.F, WT.cens.T, lambda, lambda.fraction) {
  # Function to sample alternate wait times as if number of cells were lambda.fraction x current estimate
  WT.cens.F[,2] <- rexp(nrow(WT.cens.F), rate = lambda.fraction * lambda)
  WT.cens.T[,2] <- rexp(nrow(WT.cens.T), rate = lambda.fraction * lambda)
  return(rowSums(rbind(WT.cens.F, WT.cens.T)))
}

######################################################################
# Run the Gibbs sampler

age.inc.gs <- rcc.age.gs(spor.ages = sporadic.incidence, vhl.ages = vhl.incidence, chr3p.ages = t3_5.abs.time$triploid.est, d.sf=100, iterations=rcc.iters, g.sf = 10)

save(age.inc.gs, file="age_inc_gs.RData")
aig.lambda <- age.inc.gs$lambda
save(aig.lambda, file="age_inc_gs-lambda.RData")
rm(age.inc.gs)

load("age_inc_gs.RData")

par(mfrow=c(2,2))
plot(1:rcc.iters, age.inc.gs[[1]], type="l", xlab="Iteration", ylab="Value drawn", main = "alpha_1")
plot(1:rcc.iters, age.inc.gs[[2]], type="l", xlab="Iteration", ylab="Value drawn", main = "beta_1")
plot(1:rcc.iters, age.inc.gs[[3]], type="l", xlab="Iteration", ylab="Value drawn", main = "alpha_3")
plot(1:rcc.iters, age.inc.gs[[4]], type="l", xlab="Iteration", ylab="Value drawn", main = "beta_3")

par(mfrow=c(1,1))
plot(1:rcc.iters, age.inc.gs[[11]] / vhl.driver.mutation.rate, type="l", xlab="Iteration", ylab="Value drawn", main = "Convergence plot of number of cells with 3p loss")

hist(age.inc.gs[[11]][rcc.burn.in:rcc.iters] / vhl.driver.mutation.rate, breaks=20, probability = TRUE, xlab = "Number of cells with 3p loss before VHL mutation", ylab="Density", main = "Posterior distribution of number of cells with 3p loss")

age.in.dens <- density(age.inc.gs[[11]][rcc.burn.in:rcc.iters] / vhl.driver.mutation.rate)
plot(age.in.dens, xlab = "Number of cells with 3p loss before VHL mutation", ylab="Density", main = "Posterior distribution of number of cells with 3p loss")
# polygon(age.in.dens, col="lightblue", border="black")
polygon(age.in.dens, col=adjustcolor("lightblue",alpha.f=0.2), border=adjustcolor("lightblue",alpha.f=0.8))

post.less0 <- apply(age.inc.gs[[12]][,(rcc.burn.in+1):rcc.iters], MARGIN = 2, FUN=sort)
post.less0 <- t(apply(post.less0, MARGIN = 1, FUN=quantile, probs=c(0.025,0.5,0.975)))
baseline.survfit <- survfit(Surv(post.less0[,2], rep(1, length(post.less0[,2]))) ~ 1)

post.less25 <- apply(age.inc.gs[[13]][,(rcc.burn.in+1):rcc.iters], MARGIN = 2, FUN=sort)
post.less25 <- t(apply(post.less25, MARGIN = 1, FUN=quantile, probs=c(0.025,0.5,0.975)))
less25.survfit.0.025 <- survfit(Surv(post.less25[,1], rep(1, length(post.less25[,1]))) ~ 1)
less25.survfit.0.5 <- survfit(Surv(post.less25[,2], rep(1, length(post.less25[,1]))) ~ 1)
less25.survfit.0.975 <- survfit(Surv(post.less25[,3], rep(1, length(post.less25[,1]))) ~ 1)
plot(NULL, xlim=c(0,90), ylim=c(0.95,1), xlab="Age (years)", ylab="Fraction free of kidney cancer", las=1, main="Decrease in clone size of 25%")
polygon(x=c(0,less25.survfit.0.025$time, rev(less25.survfit.0.975$time),0), y=c(1,less25.survfit.0.025$surv, rev(less25.survfit.0.975$surv),1), col="lightblue", border=NA)
lines(c(0,baseline.survfit$time), c(1,baseline.survfit$surv), lwd=3)
lines(c(0,less25.survfit.0.5$time), c(1,less25.survfit.0.5$surv), lwd=3, col="darkblue")
legend(10,0.97, lwd=c(3,3,0), col=c("black", "darkblue", NULL), fill=c(NA, NA, "lightblue"), border = "white", legend=c("Baseline incidence", "Incidence if clone size decreased by 25%", "95% posterior interval"), bty="n", merge=TRUE)


post.less50 <- apply(age.inc.gs[[14]][,(rcc.burn.in+1):rcc.iters], MARGIN = 2, FUN=sort)
post.less50 <- t(apply(post.less50, MARGIN = 1, FUN=quantile, probs=c(0.025,0.5,0.975)))
less50.survfit.0.025 <- survfit(Surv(post.less50[,1], rep(1, length(post.less50[,1]))) ~ 1)
less50.survfit.0.5 <- survfit(Surv(post.less50[,2], rep(1, length(post.less50[,1]))) ~ 1)
less50.survfit.0.975 <- survfit(Surv(post.less50[,3], rep(1, length(post.less50[,1]))) ~ 1)
plot(NULL, xlim=c(0,90), ylim=c(0.95,1), xlab="Age (years)", ylab="Fraction free of kidney cancer", las=1, main="Decrease in clone size of 50%")
polygon(x=c(0,less50.survfit.0.025$time, rev(less50.survfit.0.975$time),0), y=c(1,less50.survfit.0.025$surv, rev(less50.survfit.0.975$surv),1), col="lightblue", border=NA)
lines(c(0,baseline.survfit$time), c(1,baseline.survfit$surv), lwd=3)
lines(c(0,less50.survfit.0.5$time), c(1,less50.survfit.0.5$surv), lwd=3, col="darkblue")
legend(10,0.97, lwd=c(3,3,0), col=c("black", "darkblue", NULL), fill=c(NA, NA, "lightblue"), border = "white", legend=c("Baseline incidence", "Incidence if clone size decreased by 50%", "95% posterior interval"), bty="n", merge=TRUE)

post.less75 <- apply(age.inc.gs[[15]][,(rcc.burn.in+1):rcc.iters], MARGIN = 2, FUN=sort)
post.less75 <- t(apply(post.less75, MARGIN = 1, FUN=quantile, probs=c(0.025,0.5,0.975)))
less75.survfit.0.025 <- survfit(Surv(post.less75[,1], rep(1, length(post.less75[,1]))) ~ 1)
less75.survfit.0.5 <- survfit(Surv(post.less75[,2], rep(1, length(post.less75[,1]))) ~ 1)
less75.survfit.0.975 <- survfit(Surv(post.less75[,3], rep(1, length(post.less75[,1]))) ~ 1)
plot(NULL, xlim=c(0,90), ylim=c(0.95,1), xlab="Age (years)", ylab="Fraction free of kidney cancer", las=1, main="Decrease in clone size of 75%")
polygon(x=c(0,less75.survfit.0.025$time, rev(less75.survfit.0.975$time),0), y=c(1,less75.survfit.0.025$surv, rev(less75.survfit.0.975$surv),1), col="lightblue", border=NA)
lines(c(0,baseline.survfit$time), c(1,baseline.survfit$surv), lwd=3)
lines(c(0,less75.survfit.0.5$time), c(1,less75.survfit.0.5$surv), lwd=3, col="darkblue")
legend(10,0.97, lwd=c(3,3,0), col=c("black", "darkblue", NULL), fill=c(NA, NA, "lightblue"), border = "white", legend=c("Baseline incidence", "Incidence if clone size decreased by 75%", "95% posterior interval"), bty="n", merge=TRUE)


Sys.time()
sessionInfo()
