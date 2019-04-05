setwd("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/bootstrapping/bootstrap_simulated_NA/Model_wo_cons_status/Model_wo_cons_status_minus_high_ss/")

classes   = c("ACTINOPTERYGII", "CHONDRICHTHYES")     #choose: "ACTINOPTERYGII", "CHONDRICHTHYES"
filenames = c( "habitat", "constatus")  #choose: lifehistory  habitat  constatus  lifehistorynolog  ##see below for column/data choice implications##

combined  = FALSE

for(cc in 1:length(classes)){ #iterate over classes
  class = classes[cc]
  for(f in 1:length(filenames)){ #iterate over models/predictors
    filename  = filenames[f]
    
    #read in original data
    data = read.table("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/bootstrapping/full_dataset_simulated_nalleles_2_21_18.csv", header=TRUE, sep=",")
    data = data[data$class %in% class,,drop=FALSE]
    
    #set up t/nt variable
    data$Istatus = rep(NA, nrow(data))
    data$Istatus[data$status=="LC" | data$status=="NT"] = "not_threatened"
    data$Istatus[data$status=="EN" | data$status=="CR" | data$status=="VU" | data$status=="EX"] = "threatened"
    
    if(filename=="lifehistory"){
      modeldata = read.table(paste(class, "_", filename, "_coefficients.txt", sep=""), header=TRUE, sep="\t", row.names = NULL)
      colnames(modeldata) = c("estlab", "estimate", "mean", "l95CI", "u95CI", "sd", "dataset", "yvalue", "xvalue")
      
      data$Max_Fec = log(data$Max_Fec)
      data$AAAdultWeight_g = log(data$AAAdultWeight_g)
      datasets     = list(freshwater = data[data$habitat=="freshwater",], marine = data[data$habitat=="marine",])
      responsecols = c(9,37)               #hetM, simulated_nAlleles
      predictcols  = c(20,22,30,31)        #Max_Fec, MatAge_Min, AAAdultWeight_g , AAMaxLongevity
      intercept    = c("Istatusnot_threatened", "Istatusthreatened")
      slopes       = c("pred:Istatusnot_threatened", "pred:Istatusthreatened")
      
      #plotting parameters
      cols = c("goldenrod1", "darkorange", "dodgerblue1", "navyblue")
      xl = c(0, 0, 0, 0)
      xu = c(20, 20, 20,125)
      yl = c(0,0)
      yu = c(1,20)
      
      #plots
      if(combined==TRUE){
        pdf(paste(class, filename, "combined.pdf", sep=""), width=10, height=10, useDingbats=FALSE, onefile=TRUE)
        par(mfrow=c(2,2))
      }else{
        pdf(paste(class, filename, ".pdf", sep=""), width=5, height=5, useDingbats=FALSE, onefile=TRUE)
        par(mfrow=c(1,1))
      }
      for(p in 1:length(predictcols)){
        for(r in 1:length(responsecols)){
          #track colors
          c = 1
          for(d in 1:length(datasets)){
            #start plot
            
            plot(-100, -100, xlim=c(xl[p], xu[p]), ylim=c(yl[r], yu[r]), xlab=colnames(data)[predictcols[p]], ylab=colnames(data)[responsecols[r]])
            
            #isoalte needed data
            temp = datasets[[d]]
            nt   = temp[temp$Istatus=="not_threatened",,drop=FALSE]
            t    = temp[temp$Istatus=="threatened",,drop=FALSE]
            linedata = modeldata[modeldata$dataset==names(datasets)[[d]] & modeldata$yvalue==colnames(data)[responsecols[r]] & modeldata$xvalue==colnames(data)[predictcols[p]],,drop=FALSE ]
            
            #nt
            if(nrow(linedata)>0){
              if((linedata$l95CI[linedata$estlab==slopes[1]] > 0 & linedata$u95CI[linedata$estlab==slopes[1]] > 0) |
                 (linedata$l95CI[linedata$estlab==slopes[1]] < 0 & linedata$u95CI[linedata$estlab==slopes[1]] < 0) ){
                setlwd = 4
              }else{
                setlwd = 2
              }
            }else{
              setlwd = 1
            }
            points(x=nt[,predictcols[p]], y=nt[,responsecols[r]], pch=21, bg=cols[c], col="grey30", cex=1.25)
            if(nrow(linedata)>0){
              abline(a=linedata$estimate[linedata$estlab==intercept[1]], b=linedata$estimate[linedata$estlab==slopes[1]], col=cols[c], lty=1, lwd=setlwd)
            }
            c = c + 1
            
            #t
            if(nrow(linedata)>0){
              if((linedata$l95CI[linedata$estlab==slopes[2]] > 0 & linedata$u95CI[linedata$estlab==slopes[2]] > 0) |
                 (linedata$l95CI[linedata$estlab==slopes[2]] < 0 & linedata$u95CI[linedata$estlab==slopes[2]] < 0) ){
                setlwd = 4
              }else{
                setlwd = 2
              }
            }else{
              setlwd = 1
            }
            points(x=t[,predictcols[p]], y=t[,responsecols[r]], pch=21, bg=cols[c], col="grey30", cex=1.25)
            if(nrow(linedata)>0){
              abline(a=linedata$estimate[linedata$estlab==intercept[2]], b=linedata$estimate[linedata$estlab==slopes[2]], col=cols[c], lty=1, lwd=setlwd)
            }
            c = c + 1
          }#d
        }#r
      }#p
      dev.off()
    }
    
    if(filename=="habitat"){
      if(nrow(data[data$habitat=="freshwater",,drop=FALSE]) > 5 & nrow(data[data$habitat=="marine",,drop=FALSE]) > 5 ){
        datasets     = list(data=data)
        responsecols = c(9,37)         #hetM, simulated_NAlleles
        predictcols  = c(26)           #habitat
        estimates    = c("predfreshwater", "predmix", "predmarine")
      }else{
        next
      }
      modeldata = read.table(paste(class, "_", filename, "_coefficients.txt", sep=""), header=TRUE, sep="\t", row.names = NULL)
      colnames(modeldata) = c("estlab", "estimate", "mean", "l95CI", "u95CI", "sd", "dataset", "yvalue", "xvalue")
      
      #plotting parameters
      cols = c("goldenrod1", "chartreuse3", "dodgerblue1")
      xl = c(0.25)
      xu = c(1.75)
      yl = c(0,0)
      yu = c(1,30)
      ss = c(0.2, 5)
      
      #plots
      if(combined==TRUE){
        pdf(paste(class, filename, "combined.pdf", sep=""), width=6, height=5, useDingbats=FALSE, onefile=TRUE)
        par(mfrow=c(1,2))
      }else{
        pdf(paste(class, filename, ".pdf", sep=""), width=3, height=5, useDingbats=FALSE, onefile=TRUE)
        par(mfrow=c(1,1))
      }
    
    for(p in 1:length(predictcols)){  
      for(r in 1:length(responsecols)){
          #track colors
          c = 1
          
          #start plot
          plot(-100, -100, xlim=c(xl[1], xu[1]), ylim=c(yl[r], yu[r]), xlab=colnames(data)[predictcols[p]], ylab=paste("coefficient estimate: ", colnames(data)[responsecols[r]], sep=""), axes=FALSE)
          axis(side=2, at=seq(yl[r], yu[r], ss[r]), labels=seq(yl[r], yu[r], ss[r]), tick=TRUE, pos=0.25)
 
          #isoalte needed data
          modeloutput = modeldata[modeldata$yvalue==colnames(data)[responsecols[r]],,drop=FALSE ]
          
          #reorder to plot correctly
          modeloutput = as.data.frame(rbind(modeloutput[1,], modeloutput[3,], modeloutput[2,]))
          
          #generate plot
          increments = c(0.5, 1, 1.5)
          width = 0.2
          for(a in 1:length(increments)){
            polygon(x=c(increments[a]-width, increments[a]+width, increments[a]+width, increments[a]-width), y=c(0,0,modeloutput$estimate[a],modeloutput$estimate[a]), col=cols[a])
            segments(x0=increments[a], y0=modeloutput$l95CI[a], x1=increments[a], y1=modeloutput$u95CI[a], lwd=1.5)
            segments(x0=increments[a]-0.1, y0=modeloutput$l95CI[a], x1=increments[a]+0.1, y1=modeloutput$l95CI[a], lwd=1.5)
            segments(x0=increments[a]-0.1, y0=modeloutput$u95CI[a], x1=increments[a]+0.1, y1=modeloutput$u95CI[a], lwd=1.5)
          }
       }#r
    }#p
      dev.off()
    }
    
    if(filename=="constatus"){
      if(nrow(data[data$habitat=="freshwater",,drop=FALSE]) >5 & nrow(data[data$habitat=="marine",,drop=FALSE]) > 5 ){
        datasets     = list(freshwater = data[data$habitat=="freshwater",], marine = data[data$habitat=="marine",])
      }else if(nrow(data[data$habitat=="freshwater",,drop=FALSE]) > 5 & nrow(data[data$habitat=="marine",,drop=FALSE]) < 5 ){
        datasets     = list(freshwater = data[data$habitat=="freshwater",])
      }else if(nrow(data[data$habitat=="freshwater",,drop=FALSE]) < 5 & nrow(data[data$habitat=="marine",,drop=FALSE]) > 5 ){
        datasets     = list(marine = data[data$habitat=="marine",])
      }else{
        next
      }
      modeldata = read.table(paste(class, "_", filename, "_coefficients.txt", sep=""), header=TRUE, sep="\t", row.names = NULL)
      colnames(modeldata) = c("estlab", "estimate", "mean", "l95CI", "u95CI", "sd", "dataset", "yvalue", "xvalue")
      responsecols = c(9,37)         #hetM, simulated_nAlleles
      predictcols  = c(39)   #Istatus
      estimates    = c("prednot_threatened", "predthreatened")
      
      #plotting parameters
      cols = c("goldenrod1", "darkorange", "dodgerblue1", "navyblue")
      xl = c(0)
      xu = c(1)
      yl = c(0,0)
      yu = c(1,30)
      ss = c(0.2, 5)
      
      if(length(names(datasets))==1){
        cols = cols[3:4]
      }
      
      if(combined==TRUE){
        pdf(paste(class, filename, "combined.pdf", sep=""), width=6, height=10, useDingbats=FALSE, onefile=TRUE)
        par(mfrow=c(2,2))
      }else{
        pdf(paste(class, filename, ".pdf", sep=""), width=3, height=5, useDingbats=FALSE, onefile=TRUE)
        par(mfrow=c(1,1))
      }
      #plots
      for(r in 1:length(responsecols)){
        #isoalte needed data
        modeloutput = modeldata[modeldata$yvalue==colnames(data)[responsecols[r]],,drop=FALSE ]
        
        c = 1
        for(d in 1:length(names(datasets))){
          #start plot
          plot(-100, -100, xlim=c(xl[1], xu[1]), ylim=c(yl[r], yu[r]), xlab=colnames(data)[predictcols[1]], ylab=paste("coefficient estimate: ", colnames(data)[responsecols[r]], sep=""), axes=FALSE)
          axis(side=2, at=seq(yl[r], yu[r], ss[r]), labels=seq(yl[r], yu[r], ss[r]), tick=TRUE, pos=0)
          
          #isolate data
          toplot = modeloutput[modeloutput$dataset==names(datasets)[[d]],,drop=FALSE]
          
          #generate plot
          increments = c(0.25, 0.75)
          width = 0.2
          
          for(a in 1:length(increments)){
            polygon(x=c(increments[a]-width, increments[a]+width, increments[a]+width, increments[a]-width), y=c(0,0,toplot$estimate[a],toplot$estimate[a]), col=cols[c])
            segments(x0=increments[a], y0=toplot$l95CI[a], x1=increments[a], y1=toplot$u95CI[a], lwd=1.5)
            segments(x0=increments[a]-0.1, y0=toplot$l95CI[a], x1=increments[a]+0.1, y1=toplot$l95CI[a], lwd=1.5)
            segments(x0=increments[a]-0.1, y0=toplot$u95CI[a], x1=increments[a]+0.1, y1=toplot$u95CI[a], lwd=1.5)
            c = c + 1
          }
        }
      }#r
      dev.off()
    }
  }
}
traceback()
