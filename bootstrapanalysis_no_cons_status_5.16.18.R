library(nlme)
library(ape)
library(phytools)
library(caper)
library(geiger)
setwd("/Users/alexandermartinez/Google Drive/Fish_GeneticDiversity_Proj/bootstrapping/bootstrap_simulated_NA/Model_wo_cons_status/Model_wo_cons_status_minus_high_ss/")

#set number of bootstrapps and model choice
boots     = 10000
classes   = c("OSTEICHTHYES", "CHONDRICHTHYES")     #choose: "OSTEICHTHYES", "CHONDRICHTHYES"
filenames = c("lifehistory", "habitat", "constatus")  #choose: lifehistory  lmratio  habitat  constatus  lifehistorynolog  ##see below for column/data choice implications##

#read in data
newdata = read.table("../../../full_dataset_simulated_nalleles_5_16_18.csv", header=TRUE, sep=",")
alldata = newdata 
alldata = subset(newdata, newdata$numIndv < 1000) #commented out after calculating rarefied number of alleles to account for differences in sample size
alldata$G_sIUCN <- factor(alldata$G_sIUCN)
rownames(alldata) = as.character(alldata$G_sIUCN)

#set up t/nt variable
alldata$Istatus = rep(NA, nrow(alldata))
alldata$Istatus[alldata$status=="LC" | alldata$status=="NT"] = "not_threatened"
alldata$Istatus[alldata$status=="EN" | alldata$status=="CR" | alldata$status=="VU" | alldata$status=="EX"] = "threatened"

#set up lifespan/age at maturity variable
alldata$lifeage = alldata$AAMaxLongevity/alldata$MatAge_Min

for(c in 1:length(classes)){ #iterate over classes
  class = classes[c]
  data  = alldata[alldata$class %in% class,,drop=FALSE]
  for(f in 1:length(filenames)){ #iterate over models/predictors
    filename = filenames[f]
    
    #set up model and divide data appropriately
    if(filename=="lifehistory"){
      #log variables
      data$Max_Fec = log(data$Max_Fec)
      data$AAAdultWeight_g = log(data$AAAdultWeight_g)
      
      temp_fresh = data[data$habitat=="freshwater",]
        temp_fresh$hab_short = rep("f", length(temp_fresh$G_sIUCN))
        
      temp_mar = data[data$habitat=="marine",]
        temp_mar$hab_short = rep("m", length(temp_mar$G_sIUCN))
        
      data = rbind(temp_fresh,temp_mar)
      
      datasets     = list(data)
      tomodel      = resp ~ pred*hab_short-1-pred
      responsecols = c(9,37)               #hetM, simulated_nAlleles
      predictcols  = c(20,22,30,31)        #Max_Fec, MatAge_Min, AAAdultWeight_g, AAMaxLongevity
    }
    if(filename=="lmratio"){
      temp_fresh = data[data$habitat=="freshwater",]
        temp_fresh$hab_short = rep("f", length(temp_fresh$G_sIUCN))
      
      temp_mar = data[data$habitat=="marine",]
        temp_mar$hab_short = rep("m", length(temp_mar$G_sIUCN))
      
      data = rbind(temp_fresh,temp_mar)
      
      datasets     = data
      tomodel      = resp ~ pred*hab_short-1-pred
      responsecols = c(9,37)               #hetM, simulated_nAlleles
      predictcols  = c(20,22,30,31)        #Max_Fec, MatAge_Min, AAAdultWeight_g, AAMaxLongevity
    }
    if(filename=="lifehistorynolog"){
      temp_fresh = data[data$habitat=="freshwater",]
      temp_fresh$hab_short = rep("f", length(temp_fresh$G_sIUCN))
      
      temp_mar = data[data$habitat=="marine",]
      temp_mar$hab_short = rep("m", length(temp_mar$G_sIUCN))
      
      data = rbind(temp_fresh,temp_mar)
      
      datasets     = data
      tomodel      = resp ~ pred*hab_short-1-pred
      responsecols = c(9,37)               #hetM, simulated_nAlleles
      predictcols  = c(40)                 #lifespan:maturity
    }
    if(filename=="habitat"){
      if(nrow(data[data$habitat=="freshwater",,drop=FALSE]) > 5 & nrow(data[data$habitat=="marine",,drop=FALSE]) > 5 ){
        datasets     = list(data=data)
        tomodel      = resp ~ pred-1
        responsecols = c(9,37)         #hetM, simulated_nAlleles
        predictcols  = c(26)           #family_habitat
      }else{
        next
      }
     
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
      
      tomodel      = resp ~ pred-1
      responsecols = c(9,37)         #hetM, simulated_nAlleles
      predictcols  = c(39)   #Istatus
    }
    
    #set up output file
    write.table(date(), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(paste("model formula:", as.character(tomodel[3]), "\n", sep=" "), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE,  col.names=FALSE, row.names=FALSE, quote=FALSE)
    
    #set up output object
    allout = NULL
    
    #set up phylo tree from taxonomy, set up comp data object for analysis
    as.phylo.formula2 = function (x, data = parent.frame(), ...){
      err <- "Formula must be of the kind \"~A1/A2/.../An\"."
      if (length(x) != 2) 
        stop(err)
      if (x[[1]] != "~") 
        stop(err)
      f <- x[[2]]
      taxo <- list()
      while (length(f) == 3) {
        if (f[[1]] != "/") 
          stop(err)
        if (!is.factor(data[[deparse(f[[3]])]])) 
          stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
        taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
        if (length(f) > 1) 
          f <- f[[2]]
      }
      if (!is.factor(data[[deparse(f)]])) 
        stop(paste("Variable", deparse(f), "must be a factor."))
      taxo[[deparse(f)]] <- data[[deparse(f)]]
      taxo.data <- as.data.frame(taxo)
      leaves.names <- as.character(taxo.data[, 1])
      taxo.data[, 1] <- 1:nrow(taxo.data)
      f.rec <- function(subtaxo) {
        u <- ncol(subtaxo)
        levels <- unique(subtaxo[, u])
        if (u == 1) {
          if (length(levels) != nrow(subtaxo)) 
            warning("Error, leaves names are not unique.")
          return(as.character(subtaxo[, 1]))
        }
        t <- character(length(levels))
        for (l in 1:length(levels)) {
          x <- f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)])
          t[l] <- paste("(", paste(x, collapse = ","), ")", sep = "")
        }
        return(t)
      }
      string <- paste("(", paste(f.rec(taxo.data), collapse = ","),");", sep = "")
      phy <- read.newick(text = string) ## so that singles will be read without error
      phy$edge.length <- rep(1,nrow(phy$edge))
      phy <- collapse.singles(phy)
      phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
      return(phy)
    }
    
    #Make phylogenetic tree in Newick format for all data
    atree  = as.phylo.formula2(~class/order/family/genus/G_sIUCN, data=data) 
    
    for(d in 1:length(datasets)){               #data devisions defined above, datasets
      #record dataset division
      write.table(paste("data:", names(datasets)[d], "\n", sep=" "), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE,  col.names=FALSE, row.names=FALSE, quote=FALSE)
      for(y in 1:length(responsecols)){         #y values, responsecols
        for(x in 1:length(predictcols)){        #x values, predictcols
          #record variables
          write.table(paste("variables: y=",  colnames(datasets[[d]])[responsecols[y]], " x=", colnames(datasets[[d]])[predictcols[x]], sep=""), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE,  col.names=FALSE, row.names=FALSE, quote=FALSE)
          
          #remove species with missing GD data
          tdata = datasets[[d]]
          temp  = tdata[!is.na(tdata[,responsecols[y]]), ,drop=FALSE]
          
          #remove species with missing x data
          temp   = temp[!is.na(temp[,predictcols[x]]), ,drop=FALSE]
          temp   = temp[!is.na(temp$hab_short), , drop=FALSE]
          tocomp = data.frame(G_sIUCN=temp$G_sIUCN, pred = temp[,predictcols[x]], resp=temp[,responsecols[y]], hab_short=temp$hab_short)
          rownames(tocomp) = tocomp$G_sIUCN
          
          #check for enough data
          if(nrow(tocomp)<1){
            next
            write.table(paste("not enough data", "\n", sep=""), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE,  col.names=FALSE, row.names=FALSE, quote=FALSE)
            if(table(tocomp$hab_short)[1]<5 | table(tocomp$hab_short)[2]<5){
              write.table(paste("not enough data", "\n", sep=""), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE,  col.names=FALSE, row.names=FALSE, quote=FALSE)
              next
            }
          }

          #remove species that are not going to be analysed from phylo tree
          toremove = as.character(data$G_sIUCN[!(data$G_sIUCN %in% tocomp$G_sIUCN)])
          
          #update tree
          temptree = drop.tip(atree, toremove)
          temptree = multi2di(temptree)
          
          #run pgls, save output
          outpgls  = gls(tomodel, correlation = corBrownian(phy = temptree), data = tocomp, method = "ML")
          summary(outpgls)
          adata = c(coef(summary(outpgls))[,1], coef(summary(outpgls))[,2], coef(summary(outpgls))[,3], coef(summary(outpgls))[,4])
          
          #bootstrap coefficient estimates with pgls
          boot_out = NULL
          missingboot = 0
          fresh_tocomp  = tocomp[tocomp$hab_short=="f",,drop=FALSE]
          mar_tocomp   = tocomp[tocomp$hab_short=="m",,drop=FALSE]
          
          for(b in 1:boots){
            #take sample from full dataset
            rowstotakeNT = round((nrow(fresh_tocomp) * 0.8), 0)
            rowstotakeT  = round((nrow( mar_tocomp) * 0.8), 0)
            tocompb      = rbind(fresh_tocomp[sample(seq(1, nrow(fresh_tocomp), 1), rowstotakeNT, replace=FALSE), ,drop=FALSE],
                                  mar_tocomp[sample(seq(1, nrow( mar_tocomp), 1), rowstotakeT,  replace=FALSE), ,drop=FALSE])
            
            #check for enough data
            if(table(tocompb$hab_short)[1]<5 | table(tocompb$hab_short)[2]<5){
              missingboot = missingboot + 1
              next
            }
            
            #prune tree
            toremoveb = as.character(tocomp$G_sIUCN[!(tocomp$G_sIUCN %in% tocompb$G_sIUCN)])
            temptreeb = drop.tip(temptree, toremoveb)
            temptreeb = multi2di(temptreeb)
            
            #run model and save output
            outpgls  = gls(tomodel, correlation = corBrownian(phy = temptreeb), data = tocompb, method = "ML")
            boot_out = rbind(boot_out, c(coef(summary(outpgls))[,1], coef(summary(outpgls))[,2], coef(summary(outpgls))[,3], coef(summary(outpgls))[,4]))
          }
          #set up output matrix
          p = matrix(ncol=5, nrow=length(coef(summary(outpgls))[,1]))
          
          #calc CI for each coefficient
          for(r in 1:nrow(p)){
            p[r,1] = adata[r]
            p[r,2] = mean(boot_out[,r])
            p[r,3] = quantile(boot_out[,r], probs=0.08)
            p[r,4] = quantile(boot_out[,r], probs=0.92)
            p[r,5] = sd(boot_out[,r])
          }
          rownames(p) = names(adata)[1:nrow(p)]
          colnames(p) = c("estimate", "mean", "lSE", "uSE", "sd")
          
          #write output
          write.table(p, paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE, col.names=TRUE, row.names=TRUE, quote=FALSE)
          
          #note singular bootstrap samples
          if(missingboot != 0){
            write.table(paste(missingboot, " bootstrap sample(s) unusable - singular covariance", sep=""), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE, col.names=TRUE, row.names=TRUE, quote=FALSE)
          }
          write.table(paste("\n", sep=""), paste(class, "_", filename,"_bootpgls.txt", sep=""), sep="\t", append=TRUE, col.names=TRUE, row.names=TRUE, quote=FALSE)
          
          p = cbind(p, rep(names(datasets[d]), nrow(p)), rep(colnames(datasets[[d]])[responsecols[y]], nrow(p)), rep(colnames(datasets[[d]])[predictcols[x]], nrow(p)))
          allout = rbind(allout, p)
        }
      }
    }
    colnames(allout) = c("estimate", "mean", "l95CI", "u95CI", "sd", "yvalue", "xvalue")
    write.table(allout, paste(class, "_", filename,"_coefficients.txt", sep=""), sep="\t", append=FALSE, col.names=TRUE, row.names=TRUE, quote=FALSE)
  }
}

