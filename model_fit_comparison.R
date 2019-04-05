#____________________________________________________________________________________________________#
#script created on 5/21/18 by ASM
#script plots heterozygosity vs both mean num alleles and mean rarefied num alleles, fits models to both 
#____________________________________________________________________________________________________#

library(palettetown)

setwd("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/bootstrapping/bootstrap_simulated_NA/Model_wo_cons_status/Model_wo_cons_status_minus_high_ss/")

#read in data; modify as necessary
newdata = read.table("../../../full_dataset_simulated_nalleles_5_16_18.csv", header=TRUE, sep=",")
alldata = newdata 
alldata = subset(newdata, newdata$numIndv < 1000) #commented out after calculating rarefied number of alleles to account for differences in sample size
alldata$G_sIUCN <- factor(alldata$G_sIUCN)

attach(alldata)

#create plots of heteroz ~ allelesM; ~ simulated_nAlleles
par(mfrow=c(2,1))
with(alldata, plot(hetM, allelesM, ylim = c(0,100), main = "Heterozygosity vs. mean number of alleles", pch = 21, bg="blue" ))
with(alldata, plot(hetM, simulated_nAlleles, main = "Heterozygosity vs. mean rarefied number of alleles", pch = 21, bg="blue"))

#Fit a quasipoisson family regression; plot trend-line & Rsquared
m_all <- glm(allelesM ~ hetM, family = quasipoisson(link = "log"))
m_rare <- glm(simulated_nAlleles ~ hetM, quasipoisson(link = "log"))

summ_all <- summary(m_all)
sum_rare <- summary(m_rare)

#create plots of heteroz ~ allelesM; ~ simulated_nAlleles
par(mfrow=c(2,1))
pokecols <- ichooseyou(pokemon = 9, spread = 1)

#plot quasipoisson trendline for mean num alleles data
with(alldata, plot(hetM, allelesM, ylim = c(0,110), main = "Heterozygosity vs. mean number of alleles", ylab = "Mean num. alleles",pch = 21, bg=pokecols ))
hetaxis <- seq(0,1,0.01)
y <- predict(m_all, list(hetM = hetaxis))
lines(hetaxis, exp(y), lwd = 2, col = "red")

#plot quasipoisson trendline for rarefied mean num alleles data
with(alldata, plot(hetM, simulated_nAlleles, ylim = c(0,110), main = "Heterozygosity vs. mean rarefied number of alleles",ylab = "Rarefied mean num. alleles", pch = 21, bg=pokecols))
y_rare <- predict(m_rare, list(hetM = hetaxis))
lines(hetaxis, exp(y_rare), lwd = 2, col = "red")

#_________________________________Repeat above; change predictor variable to sample size____________________________________________#

#plot sample size vs. num alleles and rarefied num alleles
with(alldata, plot(numIndv, allelesM, ylim = c(0,100), main = "Sample Size vs. mean number of alleles", pch = 21, bg=pokecols, xlab = "Number of individuals sampled", ylab = "Mean number of alleles" ))
with(alldata, plot(numIndv, simulated_nAlleles, main = "Sample size vs. mean rarefied number of alleles", pch = 21, bg=pokecols, xlab = "Number of individuals sampled", ylab = "Rarefied mean number of alleles"))

#Fit a quasipoisson family regression; plot trend-line & Rsquared
ss_all <- glm(allelesM ~ numIndv, family = quasipoisson(link = "log"))
ss_rare <- glm(simulated_nAlleles ~ numIndv, quasipoisson(link = "log"))

summ_ss_all <- summary(ss_all)
summ_ss_rare <- summary(ss_rare)

#plot quasipoisson trendline for mean num alleles data
with(alldata, plot(numIndv, allelesM, ylim = c(0,110), main = "Sample size vs. mean number of alleles", ylab = "Mean num. alleles", xlab = "Number of individuals sampled",pch = 21, bg=pokecols ))
hetaxis <- seq(0,1000,1)
y <- predict(ss_all, list(numIndv = hetaxis))
lines(hetaxis, exp(y), lwd = 2, col = "red")

#plot quasipoisson trendline for rarefied mean num alleles data
with(alldata, plot(numIndv, simulated_nAlleles, ylim = c(0,110), main = "Sample size vs. mean rarefied number of alleles",ylab = "Rarefied mean num. alleles", xlab = "Number of individuals sampled", pch = 21, bg=pokecols))
y_rare <- predict(ss_rare, list(numIndv = hetaxis))
lines(hetaxis, exp(y_rare), lwd = 2, col = "red")

#_________________________________Repeat abvoe; change GLM to model Gaussian family rather than quasi poisson____________________________________________#

#subset data based on sample sizes
d50 <- subset(newdata, newdata$numIndv < 51)
d100 <- subset(newdata, newdata$numIndv < 101)
d1000 <- subset(newdata, newdata$numIndv < 1001)
dsets = list(d50,d100,d1000)
ssizes = c("50","100","1000")
ss_num = c(50,100,1000)

for(i in 1:length(dsets)){
  df = dsets[[i]]
  attach(df)  
  
  #plot sample size vs. num alleles and rarefied num alleles
  with(df, plot(numIndv, allelesM, ylim = c(0,110), main = paste("Mean NA: SS =",ssizes[i], sep=" "), pch = 21, bg=pokecols, xlab = "Number of individuals sampled", ylab = "Mean number of alleles" ))
  with(df, plot(numIndv, simulated_nAlleles, ylim = c(0,110), main = paste("Rarefied mean NA: SS =",ssizes[i], sep=" "), pch = 21, bg=pokecols, xlab = "Number of individuals sampled", ylab = "Rarefied mean number of alleles"))
  
  #Fit a gaussian family regression; plot trend-line & Rsquared
  ss_all <- with(df, glm(allelesM ~ numIndv, family = "gaussian"))
  ss_rare <- with(df, glm(simulated_nAlleles ~ numIndv, family = "gaussian"))
  
  summ_ss_all <- summary(ss_all)
  summ_ss_rare <- summary(ss_rare)
  print(summ_ss_all)
  print(summ_ss_rare)
  
  #plot gaussian trendline for mean num alleles data
  pdf(file = paste("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/allele_rarefaction/alleles_ss",ssizes[i],".pdf",sep=""))
  with(df, plot(numIndv, allelesM, ylim = c(5,110), main = paste("Sample size vs. mean NA:", ssizes[i], sep = " " ), ylab = "Mean num. alleles", xlab = "Number of individuals sampled",pch = 21, bg=pokecols, type = "n"))
  hetaxis <- seq(0,ss_num[i],ss_num[i]/100)
  y <- predict(ss_all, list(numIndv = hetaxis))
  lines(hetaxis, y, lwd = 3, col = "red")
  with(df, points(numIndv,allelesM, pch = 21, bg = pokecols))
  #lines(hetaxis, exp(y), lwd = 2, col = "red") # use this is family = quasipoisson
  dev.off()
  
  #plot gaussian trendline for rarefied mean num alleles data
  pdf(file = paste("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/allele_rarefaction/raref_alleles_ss",ssizes[i],".pdf",sep=""))
  plot(numIndv, simulated_nAlleles, ylim = c(0,110), xlim = c(8,ss_num[i]), main = paste("Sample size vs. rarefied mean NA:", ssizes[i], sep = " " ),ylab = "Rarefied mean num. alleles", xlab = "Number of individuals sampled", pch = 21, bg=pokecols, type = "n")
  y_rare <- predict(ss_rare, list(numIndv = hetaxis))
  lines(hetaxis, y_rare, lwd = 3, col = "red")
  points(numIndv, simulated_nAlleles, pch = 21, bg = pokecols)
  #lines(hetaxis, exp(y), lwd = 2, col = "red") # use this is family = quasipoisson
  dev.off()
  
  detach(df)
}




# par(mfrow = c(1,2))
# plot(m_all, 2)
# plot(m_rare, 2)
# 
# bc1 <- boxcox(m_all)
# r1 <- which(bc1$y == max(bc1$y))
# exp1 <- bc1$x[r1]
# 
# bc2 <- boxcox(m_rare)
# r2 <- which(bc2$y == max(bc2$y))
# exp2 <- bc2$x[r2]
# 
# #transform data
# mn_all <- alldata$allelesM^(exp1)
# rare_all <- alldata$simulated_nAlleles^(exp2)
# 
# 
# shapiro.test(mn_all)
# shapiro.test(rare_all)