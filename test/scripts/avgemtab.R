# compute average of emtabs for different trials and save in single  table
# call like :  R --vanilla --slave --args ../results/em/something/emtab 10 < avgemtab.R > avgemtab.Rout
argv <- commandArgs()
offset <- which(argv == "--args")
baseName <- argv[offset + 1]
count <- as.numeric(argv[offset + 2])

ithTable <- read.table(paste(baseName, "1", ".txt", sep=""), header=FALSE)
maxSize <- nrow(ithTable)
for (i in 2:count) {
  ithTable <- read.table(paste(baseName, i, ".txt", sep=""), header=FALSE)
  maxSize <- max(maxSize, nrow(ithTable))
}

ithTable <- read.table(paste(baseName, "1", ".txt", sep=""), header=FALSE)
ithSize <- nrow(ithTable)
if (ithSize < maxSize) {
  for (j in (ithSize+1):maxSize) {
    ithTable <- rbind(ithTable, ithTable[nrow(ithTable),])
    ithTable$V1[j] <- j-1
    rownames(ithTable)[j] <- j
  }
}
total <- ithTable

for (i in 2:count) {
  print(i)
  ithTable <- read.table(paste(baseName, i, ".txt", sep=""), header=FALSE)
  ithSize <- nrow(ithTable)
  if (ithSize < maxSize) {
    for (j in (ithSize+1):maxSize) {
      ithTable <- rbind(ithTable, ithTable[nrow(ithTable),])
      ithTable$V1[j] <- j-1
      rownames(ithTable)[j] <- j
    }
  }
  total <- total + ithTable
}

total <- total / count
write.table(total, paste(baseName, "Avg.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE)
