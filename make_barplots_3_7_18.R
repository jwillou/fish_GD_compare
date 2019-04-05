setwd("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/bootstrapping/bootstrap_simulated_NA/Model_wo_cons_status/Model_wo_cons_status_minus_high_ss/")

class   = ("ACTINOPTERYGII")     #choose: "ACTINOPTERYGII", "CHONDRICHTHYES"
filename = c("lifehistory")  #choose: lifehistory  habitat  constatus  lifehistorynolog  ##see below for column/data choice implications##

combined  = FALSE

modeldata = read.table(paste(class, "_", filename, "_coefficients.txt", sep=""), header=TRUE, sep="\t", row.names = NULL)
colnames(modeldata) = c("estlab", "estimate", "mean", "l95CI", "u95CI", "sd", "yvalue", "xvalue")

intercept    = c("hab_shortf", "hab_shortm")
slopes       = c("pred:hab_shortf", "pred:hab_shortm")

#plots
if(combined==TRUE){
  pdf(paste(class, filename, "combinedUPDATED.pdf", sep=""), width=10, height=10, useDingbats=FALSE, onefile=TRUE)
  par(mfrow=c(2,2))
}else{
  pdf(paste(class, filename, "UPDATED.pdf", sep=""), width=3, height=5, useDingbats=FALSE, onefile=TRUE)
  par(mfrow=c(1,1))
}

#datasets = unique(modeldata$dataset)
responsecols = unique(modeldata$yvalue)
predictcols  = unique(modeldata$xvalue)

#plotting parameters
cols = c("goldenrod1","dodgerblue1")
xl = c(0.25)
xu = c(1.75)
yl = c(-0.04,-2.25,-0.04,-2.25,-0.05,-0.35,-0.005,-0.35)
yu = c( 0.02, 0.25, 0.02,0.25, 0.05, 0.35, 0.005, 0.05)
ss = c(0.01, 0.25, 0.01,0.25, 0.0025, 0.05, 0.0025, 0.05)

#variable to index y limits
limit = 0

for(p in 1:length(predictcols)){
  for(r in 1:length(responsecols)){
    #track colors
    c = 1
    limit = limit + 1
    d=1
    #for(d in 1:length(datasets)){
      #isoalte needed data
      toplot = modeldata[modeldata$yvalue==responsecols[r] & modeldata$xvalue==predictcols[p],,drop=FALSE ]
      if(nrow(toplot)<1){
        c = 3
        next
      }
      toplot = toplot[3:4,]
      
      #start plot
      space = 0.05
      
      #set up plot
      plot(-100, -100, xlim=c(xl[1], xu[1]), ylim=c(yl[limit], yu[limit]), xlab=predictcols[p], ylab=paste("coefficient estimate: ", responsecols[r], sep=""), axes=FALSE)
      segments(x0=0.25, x1=2, y0=0, y1=0)
      axis(side=2, at=seq(yl[limit], yu[limit], ss[limit]), labels=seq(yl[limit], yu[limit], ss[limit]), tick=TRUE, pos=0.25)
    
      #generate plot
      increments = c(0.5, 1)
      width = 0.2
      
      for(a in 1:length(increments)){
        polygon(x=c(increments[a]-width, increments[a]+width, increments[a]+width, increments[a]-width), y=c(0,0,toplot$estimate[a],toplot$estimate[a]), col=cols[c])
        segments(x0=increments[a], y0=toplot$l95CI[a], x1=increments[a], y1=toplot$u95CI[a], lwd=1.5)
        segments(x0=increments[a]-0.1, y0=toplot$l95CI[a], x1=increments[a]+0.1, y1=toplot$l95CI[a], lwd=1.5)
        segments(x0=increments[a]-0.1, y0=toplot$u95CI[a], x1=increments[a]+0.1, y1=toplot$u95CI[a], lwd=1.5)
        
        c = c + 1
      #}#d
    }#r
  }#p
}
dev.off()
