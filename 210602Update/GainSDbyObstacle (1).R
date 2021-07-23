library(tidyverse)
gainSD <- read_csv("dataf_FinalGain_Variance.csv")
colnames(gainSD)[1] <- "Name"
gainSD <- gainSD[order(gainSD$varE, gainSD$Ne),]
plot(gainSD$Mean, gainSD$Sd, pch=16)
pdf("gainSDbyInterventionLines.pdf", height=8, width=8)
op <- par(mfrow=c(2,2))
for (f in 4:7){
  sb <- (4:7)[-(f-3)]
  print(c(sb, f))
  gainSD <- gainSD[order(gainSD[,8], gainSD[,9], gainSD[,sb[1]], gainSD[,sb[2]], gainSD[,sb[3]], gainSD[,f]),]
  print(gainSD)
  curNxt <- c(2, 1, 2, 1)[f-3]
  mainTitle <- c("Selection on SP", "Number of Plots Eval.", "Cycle Time in Years", "nGP per Parental SP")[f-3]
  plot(gainSD$Mean, gainSD$Sd, pch=16, col=rep(curNxt:(3-curNxt), 32), xlab="Genetic Gain", ylab="Final Genetic StdDev", main=mainTitle)
  for (i in 0:31){
    lines(gainSD$Mean[i*2 + 1:2], gainSD$Sd[i*2 + 1:2], lwd=0.5, col="gray")
  }
  text(6.6, 0.93, labels=c("A", "B", "C", "D")[f-3], cex=1.5)
}
dev.off()
