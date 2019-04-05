install.packages("palettetown")
library(palettetown)

setwd("~/Downloads/")
# tab delemited text file
# make sure no decimals for average number of alleles
dat <- read.table("test2.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)
#dat <- dat[1,]


OUT <- NULL
for(dl in 1:length(dat[, 1])){
  print(dl)
  SSIZE <- dat[dl, 2]
  NALLELES <- dat[dl, 3]
  
  allele.set <- (NALLELES-10):(NALLELES+35)   
  allele.set  <- allele.set[allele.set >= 2]
  
  for(al in allele.set){
    NALLELES2 <- al
    
    #create allele freqs============================================================================#
    NLOCI <- 1
    Nloci=as.numeric(NLOCI)
    wantednumberalleles=as.numeric(NALLELES2)
    a=c(1:wantednumberalleles)
    b=rev(a)
    c=b/a
    d=sum(c)
    freqs=c/d
    lowestallele=140
    alleles2=cbind(seq(lowestallele,lowestallele+wantednumberalleles-1,1),freqs)
    afreqs = rep(alleles2[,2],Nloci)
    afreqs=cbind(sort(rep(1:Nloci,wantednumberalleles)),afreqs)
    afreqs=cbind(afreqs[,1],afreqs)
    
    #create genotypes
    afreqs <- afreqs
    #OUT=NULL
    #sims=function(sims) {
    z = 1 
    alleles2=afreqs[which(afreqs[,1]==z),]
    alleles3=cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   #table allele frequencies
    alleles3 <- alleles3[, -2]
    homos=(alleles3[,2])^2                                                          #create homozygote allele frequencies
    homos2=cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
    hets=t(combn(alleles3[,2],2))                                                   #create heterozygote allele frequencies
    hetfreq=2*(hets[,1]*hets[,2])
    hetvals=t(combn(as.character(alleles3[,1]),2))                                  #create heterozygote allele names
    hets2=cbind(hetvals,hetfreq)
    gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
    n=1000000                                                                       #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
    gfreqs1=rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            #create genotypes(by coloumn, 1 for each allele)
    gfreqs2=rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
    gtypes=cbind(gfreqs1,gfreqs2)
    
    #begin maximum likelihood
    OUT2 <- NULL
    for(i in 1:100){
      nsamp <- gtypes[sample(1:length(gtypes[, 1]), SSIZE), ]
      nsamp2 <- c(nsamp[, 1], nsamp[, 2])
      nalleles <- length(unique(nsamp2))
      out <- cbind(i, NALLELES, NALLELES2, nalleles)
      OUT2 <- rbind(OUT2, out)
    }
    mean.allele <- mean(OUT2[, 4])
    mean.diff <- mean(OUT2[, 4]- OUT2[, 2])
    
    out <- cbind(dl, NALLELES, NALLELES2, mean.allele, mean.diff)
    OUT <- rbind(OUT, out)
    write.table(out,file="Richness_output.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
  }
}

colnames(OUT) <- c("Study.No", "Reported_NA", "Simulated_NA", "Mean_sampled_NA", "Mean_difference_NA")
# doing this by brute force; need to check that no final are at max or min of alleleset

#head(OUT)

#OUT3 <- NULL
#for(r in 1:length(unique(OUT[, 1]))){
# dats <- OUT[OUT[, 1] == r, ]
# out <- dats[which.min(abs(dats[, 5])), ] 
# out <- c(out, which.min(abs(dats[, 5]))/length(dats[, 5]))
# OUT3 <- rbind(OUT3, out)
#}
  
#colnames(OUT3)[6] <- "allele.likelihood.position"
#OUT3 #anything >90% and less than 10% should be rerun
#write.table(OUT3,file="Richness_final_output.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
par(mfrow = c(1,2))
plot.cols <- ichooseyou(pokemon = "carvanha", spread = 3)

#for(i in 1:10){
  i=9
  plot.cols <- ichooseyou(pokemon = i, spread = 3)
  o1 <- subset(OUT, OUT[,1] == 1)
  o2 <- subset(OUT, OUT[,1] == 2)
  
  plot(o1[, 3], abs(o1[, 5]), cex = 2, pch = 21, bg = plot.cols[1], xlab = "Tested number of alleles", ylab="Absolute value of mean difference", main="Chaetodon capistratus")
  abline(v = dat[1,3], col = plot.cols[2], lwd = 2.5)
  row.num1 <- which(abs(o1[,5]) == min(abs(o1[,5]))) #the row in OUT corresponding to min mean difference
  rarefied.all1 <- o1[row.num1,3] #number of rarefied alleles corresponding to min mean diff
  abline(v = rarefied.all1, col = plot.cols[3], lwd = 2.5)
  
  plot(o2[, 3], abs(o2[, 5]), cex = 2, pch = 21, bg = plot.cols[1], xlab = "Tested number of alleles", ylab="Absolute value of mean difference", main="Chaetodon ornatissimus")
  abline(v = dat[2,3], col = plot.cols[2], lwd = 2.5)
  row.num2 <- which(abs(o2[,5]) == min(abs(o2[,5]))) #the row in OUT corresponding to min mean difference
  rarefied.all2 <- o2[row.num2,3] #number of rarefied alleles corresponding to min mean diff
  abline(v = rarefied.all2, col = plot.cols[3], lwd = 2.5)
  
#}


o1 <- subset(OUT, OUT[,1] == 1)
o2 <- subset(OUT, OUT[,1] == 2)

plot(o1[, 3], abs(o1[, 5]), cex = 2, pch = 21, bg = plot.cols[1], xlab = "Tested number of alleles", ylab="Absolute value of mean difference", main="")
abline(v = dat[1,3], col = plot.cols[2])
row.num1 <- which(abs(o1[,5]) == min(abs(o1[,5]))) #the row in OUT corresponding to min mean difference
rarefied.all1 <- o1[row.num1,3] #number of rarefied alleles corresponding to min mean diff
abline(v = rarefied.all1, col = plot.cols[3])

plot(o2[, 3], abs(o2[, 5]), cex = 2, pch = 21, bg = plot.cols[1], xlab = "Tested number of alleles", ylab="Absolute value of mean difference", main="")
abline(v = dat[2,3], col = plot.cols[2])
row.num2 <- which(abs(o2[,5]) == min(abs(o2[,5]))) #the row in OUT corresponding to min mean difference
rarefied.all2 <- o2[row.num2,3] #number of rarefied alleles corresponding to min mean diff
abline(v = rarefied.all2, col = plot.cols[3])

points(OUT[, 3], abs(OUT[, 5]), cex = 2, pch = 21, bg = "blue")
