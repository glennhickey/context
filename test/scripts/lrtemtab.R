# compare probabilites obtained using context elements and not
# call like :  R --vanilla --slave --args ../results/em/something/emtab 10 30 < lrtemtab.R > lrtemtab.Rout
argv <- commandArgs()
offset <- which(argv == "--args")
baseName <- argv[offset + 1]
count <- as.numeric(argv[offset + 2])

prvec <- vector(mode="numeric", length=count)
prvecNC <- vector(mode="numeric", length=count)
pvals <- vector(mode="numeric", length=count)
bic <- vector(mode="numeric", length=count)
bicNC <- vector(mode="numeric", length=count)
lenvec <- vector(mode="numeric", length=count)

for (i in 1:count) {
  ithTable <- read.table(paste(baseName, i, ".txt", sep=""), header=FALSE)
  ithTableNC <- read.table(paste(baseName, i, "_nc.txt", sep=""), header=FALSE)
  lenvec[i] <- ithTable$V2[nrow(ithTable)]
  prvec[i] <- ithTable$V4[nrow(ithTable)]
  prvecNC[i] <- ithTableNC$V4[nrow(ithTableNC)]
# pval is computed for symmetric model! change 3 to 5 for full model
  pvals[i] <- 1-pchisq(-2 * (prvecNC[i] - prvec[i]), df=3)
# bic is computed for symmetric model!  change 7/4 to 12/7 for full mode!
  bic[i] <- -2 * prvec[i] + 7 * log(lenvec[i])
  bicNC[i] <- -2 * prvecNC[i] + 4 * log(lenvec[i])
}

results=data.frame(trial=1:count, prBasic=prvecNC, prContext=prvec, bicBasic=bicNC, bicContext=bic, lrtPVal=pvals)

write.table(results, paste(baseName, "Stats.txt", sep=""), row.names=FALSE, quote=FALSE)
