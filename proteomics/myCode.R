logPlot <- function()
{
  plot(relF1_data[,2]/relF1_data[,10],type="p",col="black",ylim=c(0.013,1.25),xaxt="n",xlab="",log="y",ylab="Relative abundace to WT_R")
  axis(1, at=1:27, labels = FALSE)
  text(1:27, par("usr")[1] + 0.049, srt = 90, adj = 1,labels = labels, xpd = TRUE)
  #text(1:27, par("usr")[1] - 0.04, srt = 90, adj = 1,labels = labels, xpd = TRUE)
  lines(relF1_data[,3]/relF1_data[,10],type="p",col="red")
  lines(relF1_data[,4]/relF1_data[,10],type="p",col="orange")
  lines(relF1_data[,5]/relF1_data[,10],type="p",col="green")
  lines(relF1_data[,6]/relF1_data[,10],type="p",col="lightblue")
  lines(relF1_data[,7]/relF1_data[,10],type="p",col="blue")
  lines(relF1_data[,8]/relF1_data[,10],type="p",col="brown")
  lines(relF1_data[,9]/relF1_data[,10],type="p",col="purple")
  colvector = c("black","red","pink","orange","green","lightblue","blue","brown","purple")
  #legend(0.5,1.27,c("SMARCA2","SMARCA4.4","SMARCA4.6","SMARCC1","ARID1A.3","ARID1A.10","ARID1B","ARID2"),col=colvector,text.col=colvector,pch=1)
  legend(0.5,0.045,c("SMARCA2","SMARCA4.4","SMARCA4.6","SMARCC1","ARID1A.3","ARID1A.10","ARID1B","ARID2"),col=colvector,text.col=colvector,pch=1,cex=.75)
}

linPlot <- function()
{
  plot(relF1_data[,2]/relF1_data[,10],type="p",col="black",ylim=c(0.013,1.25),xaxt="n",xlab="",log="",ylab="Relative abundace to WT_R")
  axis(1, at=1:27, labels = FALSE)
  #text(1:27, par("usr")[1] + 0.049, srt = 90, adj = 1,labels = labels, xpd = TRUE)
  text(1:27, par("usr")[1] - 0.04, srt = 90, adj = 1,labels = labels, xpd = TRUE)
  lines(relF1_data[,3]/relF1_data[,10],type="p",col="red")
  lines(relF1_data[,4]/relF1_data[,10],type="p",col="orange")
  lines(relF1_data[,5]/relF1_data[,10],type="p",col="green")
  lines(relF1_data[,6]/relF1_data[,10],type="p",col="lightblue")
  lines(relF1_data[,7]/relF1_data[,10],type="p",col="blue")
  lines(relF1_data[,8]/relF1_data[,10],type="p",col="brown")
  lines(relF1_data[,9]/relF1_data[,10],type="p",col="purple")
  colvector = c("black","red","pink","orange","green","lightblue","blue","brown","purple")
  legend(0.5,1.27,c("SMARCA2","SMARCA4.4","SMARCA4.6","SMARCC1","ARID1A.3","ARID1A.10","ARID1B","ARID2"),col=colvector,text.col=colvector,pch=1,cex=.75)
  #legend(0.5,0.045,c("SMARCA2","SMARCA4.4","SMARCA4.6","SMARCC1","ARID1A.3","ARID1A.10","ARID1B","ARID2"),col=colvector,text.col=colvector,pch=1,cex=.75)
}

doAxis <- function(N,ylimV,logV)
{
   plot(c(-.5,N),c(.5,.5),ylim=ylimV,type="l",xaxt="n",yaxt="n",bty="n",lty=2,col="grey",xlim=c(0,N),ylab="",xlab="",log=logV)
   lines(c(-.5,N),c(1,1),lty=2,col="grey")
   lines(c(-.5,N),c(1.5,1.5),lty=2,col="grey")
   axis(2, pos=-.5, at=c(0.0,0.5,1.0,1.5), labels = TRUE)
}

errorBarPlot <- function(N, dataMean,dataSd,colorV,addV,ylimV,logV,axesV,label="",main="")
{
	barplot(dataMean,space=0,axes=axesV, log = logV, ylim = ylimV, col=colorV, add=addV,main=main)
	arrows(1:N - .5, dataMean - dataSd/1, 1:N - .5, dataMean + dataSd/1, lwd = 1.5, angle = 90, code = 3, length = 0.05)
	text(N/2 + 1,1.5,label,col="black",cex=1.5)
}

barPlot2 <- function(legendPos=1.0,labels,ddd="BRG1"){
	data <- BRG1.ab[which(BRG1.ab[,8]=="BAF_complex"),38:41]
	if (ddd != "BRG1")
		data <- ARID1A.ab[which(ARID1A.ab[,12]=="BAF_complex"),75:78]
	N <- dim(data)[1]
	dataMean <- data[,1]
	dataSd <- data[,1]
	for(i in 1:N){
            dataMean[i] = tryCatch( mean(as.numeric(data[i,]),na.rm=T),error = function(e) {0})
            dataSd[i] =   tryCatch( sd(as.numeric(data[i,]),na.rm=T), error = function(e) {0})
        }
	#doAxis(N, c(1e4,1e13), "y");
	ord <- order(dataMean,decreasing=T)
	errorBarPlot(N, dataMean[ord],dataSd[ord],"darkgrey",F,c(1e4,1e13),"y",T,"")
	text(1:N - .5, legendPos, srt = 90, adj = 1,labels = labels[ord], xpd = TRUE)
}

scatterPlot <- function(){
        data2 <- BRG1.ab[which(BRG1.ab[,8]=="BAF_complex"),38:41]
        data1 <- ARID1A.ab[which(ARID1A.ab[,12]=="BAF_complex"),75:78]
        N <- dim(data1)[1]
        dataMean1 <- data1[,1]
        dataSd1 <- data1[,1]
        for(i in 1:N){
            dataMean1[i] = tryCatch( mean(as.numeric(data1[i,]),na.rm=T),error = function(e) {0})
            dataSd1[i] =   tryCatch( sd(as.numeric(data1[i,]),na.rm=T), error = function(e) {0})
        }

	dataMean2 <- dataMean1
	for(i in 1:N){
	   x <- which(labels_BRG1_full==labels_A1A_full[i])
	   #print(x)
	   if (length(x)==0) dataMean2[i]=0
	   else dataMean2[i] = tryCatch( mean(as.numeric(data2[x[1],]),na.rm=T),error = function(e) {0})
	}
        #doAxis(N, c(1e4,1e13), "y");
        #ord <- order(dataMean,decreasing=T)
        #errorBarPlot(N, dataMean[ord],dataSd[ord],"darkgrey",F,c(1e4,1e13),"y",T,"")
        #text(1:N - .5, legendPos, srt = 90, adj = 1,labels = labels[ord], xpd = TRUE)
	plot(dataMean1,dataMean2,log="xy",xlim=c(1e5,5e10),ylim=c(1e7,1e13))
	sset <- c(1:3,6,8,10:16,19,22:24,26:29)
	#sset <- 1:29
	text(dataMean1[sset],dataMean2[sset],labels_A1A_full[sset])
	dataMean1 <- dataMean1[sset]
	dataMean2 <- dataMean2[sset]
	a1 <- log(dataMean1[which(dataMean1>0 & dataMean2>0)],10)
	a2 <- log(dataMean2[which(dataMean1>0 & dataMean2>0)],10)
	lm(a2~a1)
}


barPlot4 <- function(){
        par(mfrow=c(1,1))
	common <- intersect(BRG1.ab[,4],ARID1A.ab[,6])
        data2 <- BRG1.ab[which(BRG1.ab[,4] %in% common),c(4,38:41)]
        data1 <- ARID1A.ab[which(ARID1A.ab[,6] %in% common),c(6,75:78)]
        #print(data1[1,])
        N <- dim(data1)[1]
        dataMean1 <- data1[,2]
        #print(data1[28,])
        dataSd1 <- data1[,2]
        for(i in 1:N){

            # rescale MS1 abundances by MS3-based ratios                
            for(j in 1:4){
               colsj <- grep(paste("Abundance.F",j,sep=""),colnames(ARID1A.rr),value=F,perl=T)
               row <- which(ARID1A.rr[,6]==data1[i,1])
               if(sum(is.na(ARID1A.rr[row,colsj]))!=10)
                       data1[i,j+1] <- data1[i,j+1] * (ARID1A.rr[row,colsj[10]]/sum(ARID1A.rr[row,colsj],na.rm=T))
               else
                       data1[i,j+1] <- data1[i,j+1]/10
            }

            dataMean1[i] = tryCatch( mean(as.numeric(data1[i,2:5]),na.rm=T),error = function(e) {0})
            dataSd1[i] =   tryCatch( sd(as.numeric(data1[i,2:5]),na.rm=T), error = function(e) {0})
        }
        #dataMean1[which(is.nan(dataMean1))]=NA
        #print(dataMean1)

        dataMean2 <- dataMean1
        dataSd2 <- dataSd1

	if (1==0){
	brg_accs <- data2[,1]
        for(i in 1:N){
           x <- which(brg_accs==data1[i,1])
           #print(x)
           if (length(x)==0){
                   dataMean2[i]=NA
                   dataSd2[i] = NA
           }
           else {
                   # rescale MS1 abundances by MS3-based ratios
                   for(j in 1:4){
                       WRind = NA
                       if (j==1 || j==2) WRind = 10
                       if (j==3 || j==4) WRind = 2
                       colsj <- grep(paste("Abundance.F",j,sep=""),colnames(BRG1.rr),value=F,perl=T)
                       row <- which(BRG1.rr[,4]==data2[x[1],1])
                       if(sum(is.na(BRG1.rr[row,colsj]))!=10)
                           data2[x[1],j+1] <- data2[x[1],j+1] * (BRG1.rr[row,colsj[WRind]]/sum(BRG1.rr[row,colsj],na.rm=T))
                       else
                           data2[x[1],j+1] <- data2[x[1],j+1]/10
                   }

                   dataMean2[i] = tryCatch( mean(as.numeric(data2[x[1],2:5]),na.rm=T),error = function(e) {0})
                   dataSd2[i] =   tryCatch( sd(as.numeric(data2[x[1],2:5]),na.rm=T), error = function(e) {0})
           }
        }
	}

	oo <- order(dataMean1,decreasing=T)[1:50]
        #barplot(matrix(c(22*dataMean1[oo],dataMean2[oo]),nrow=2,byrow=T),log="y",col=c("red","blue"),beside=T,ylim=c(1e5,1e12),
        #        legend.text=c("ARID1A pulldowns","BRG1 pulldowns"),args.legend = list(x = "topright"),ylab="MS1 abundance",main="Abundance of BAF complex members (in the WT.R replicate)")
        #arrows(3*(1:N)-2 + .5, 22*(dataMean1 - dataSd1)[oo], 3*(1:N)-2 + .5, 22*(dataMean1 + dataSd1)[oo], lwd = 1.5, angle = 90, code = 3, length = 0.02)
        #arrows(3*(1:N)-1 + .5, (dataMean2 - dataSd2)[oo], 3*(1:N)-1 + .5, (dataMean2 + dataSd2)[oo], lwd = 1.5, angle = 90, code = 3, length = 0.02)

	barplot(22*dataMean1[oo],log="y",space=0)	

	if (1==1){

	accessions <- data1[oo,1]
	positions <- (1:length(oo))-.5
	pos <- c(); labels <- c();

	for(i in 1: length(oo)){
		k <- which(BRG1.ab[,4]==accessions[i])
		#if(BRG1.ab[k,8]=="BAF_complex"){
			descr <- BRG1.ab[k,5]
			res <- regexpr("GN=[^ ]*",descr)
			labels <- c(labels,substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1))
			pos <- c(pos,positions[i])
		#}
	}
         text(pos,1.2e10,labels,srt=90,cex=.7,xpd=T)
	}

	accessions[1:30]
}


barPlot3 <- function(KO="WT\\.R"){
	par(mfrow=c(1,1))
        #data2 <- BRG1.ab[which(BRG1.ab[,"Marked.as"]=="BAF_complex"),c(4,38:41)]
        #data1 <- ARID1A.ab[which(ARID1A.ab[,"Marked.as"]=="BAF_complex"),c(6,75:78)]
	data2 <- BRG1.ab2[which(BRG1.ab2[,"Marked.as"]=="BAF_complex"),c(6,106:111)]
        data1 <- ARID1A.ab2[which(ARID1A.ab2[,"Marked.as"]=="BAF_complex"),c(4,49:54)]
	N <- dim(data1)[1]
	
	labels <- rep("",N)
	seqL <- rep(0,N)
        for (i in 1:N){
          descr <- ARID1A.ab2[which(ARID1A.ab2[,"Accession"]==data1[i,1]),"Description"]
          res <- regexpr("GN=[^ ]*",descr)
          labels[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
	  seqL[i] <- length(observable.peptides(seq=BRG1.ab2[which(BRG1.ab2[,"Accession"]==data1[i,1]),"Sequence"],nmc=1)[[1]])
	  print(labels[i])
	  print(seqL[i])
        }
	#print(labels)	

        dataMean1 <- data1[,2]
        dataSd1 <- data1[,2]

	for(i in 1:N){
	    row <- which(ARID1A.rr2[,"Accession"]==data1[i,1])
	    # rescale MS1 abundances by MS3-based ratios 		
	    for(j in 1:6){
	       colsj <- grep(paste("Abundance.F",j,sep=""),colnames(ARID1A.rr2),value=F,perl=T)
	       WRind <- grep(paste("Abundance.F",j,".*",KO,sep=""),colnames(ARID1A.rr2)[colsj],value=F,perl=T)
	       if (length(WRind)==0){
		  data1[i,j+1] <- NA
	       }
	       else{
	         if(sum(is.na(ARID1A.rr2[row,colsj]))!=10) 
		       data1[i,j+1] <- data1[i,j+1] * (ARID1A.rr2[row,colsj[WRind]]/sum(ARID1A.rr2[row,colsj],na.rm=T))
	         else
		       data1[i,j+1] <- data1[i,j+1]/10
	         }
	    }

            dataMean1[i] = tryCatch( mean(as.numeric(data1[i,2:7]),na.rm=T)/seqL[i],error = function(e) {0})
            dataSd1[i] =   tryCatch( sd(as.numeric(data1[i,2:7]),na.rm=T)/seqL[i], error = function(e) {0})/2
        }

	print(dataMean1[29])
	print(dataSd1[29])

        dataMean2 <- dataMean1
	dataSd2 <- dataSd1
        for(i in 1:N){
           x <- which(data2[,1]==data1[i,1])
           if (length(x)==0){ 
		   dataMean2[i]=NA
		   dataSd2[i] = NA
	   }
           else {
		   row <- which(BRG1.rr2[,"Accession"]==data2[x[1],1])
		   # rescale MS1 abundances by MS3-based ratios
		   # file ID's in PD study : F1,F2,F3,F4,F9,F10
		   fileID <- c(1,2,3,4,9,10)
            	   for(j in 1:6){
                       colsj <- grep(paste("Abundance.F",fileID[j],"\\.",sep=""),colnames(BRG1.rr2),value=F,perl=T)
		       WRind <- grep(paste("Abundance.F",fileID[j],"\\..*",KO,sep=""),colnames(BRG1.rr2)[colsj],value=F,perl=T)
                       if (length(WRind)==0){
                  		data2[x[1],j+1] <- NA
               		}else{
		       	   if(sum(is.na(BRG1.rr2[row,colsj]))!=10)
                               data2[x[1],j+1] <- data2[x[1],j+1] * (BRG1.rr2[row,colsj[WRind]]/sum(BRG1.rr2[row,colsj],na.rm=T))
                           else
                               data2[x[1],j+1] <- data2[x[1],j+1]/10
		       }
                   }

		   dataMean2[i] = tryCatch( mean(as.numeric(data2[x[1],2:7]),na.rm=T)/seqL[i],error = function(e) {0})
		   dataSd2[i] =   tryCatch( sd(as.numeric(data2[x[1],2:7]),na.rm=T)/seqL[i], error = function(e) {0})/2
	   }
        }
	
	oo <- order(dataMean1,decreasing=T)
	cc <- 1.0
        barplot(matrix(c(cc*dataMean1[oo],dataMean2[oo]),nrow=2,byrow=T),log="y",col=c("red","blue"),beside=T,ylim=c(1e2,1e9),
		legend.text=c("ARID1A pulldowns","BRG1 pulldowns"),args.legend = list(x = "topright"),ylab="MS1 abundance",main="Abundance of BAF complex members (in the WT.R replicate)")
	arrows(3*(1:N)-2 + .5, cc*(dataMean1 - dataSd1)[oo], 3*(1:N)-2 + .5, cc*(dataMean1 + dataSd1)[oo], lwd = 1.5, angle = 90, code = 3, length = 0.02)
	arrows(3*(1:N)-1 + .5, (dataMean2 - dataSd2)[oo], 3*(1:N)-1 + .5, (dataMean2 + dataSd2)[oo], lwd = 1.5, angle = 90, code = 3, length = 0.02)
	text(3*(1:length(oo))-1,3e1,labels[oo],srt=90,cex=.7,xpd=T)
        return(1)


	bb <- c(1:2,4,6:10,12,15,17:25)
	#print(labels[bb])
	ss = sum(is.na(dataMean2))
	if (ss<29){
	print(ss)
	print(lm(log(dataMean2[bb],10) ~ log(dataMean1[bb],10)))
	
	title=paste("Abundance of BAF complex members in the ",KO," KO",sep="")
	if (KO=="WT\\.R") title = "Abundance of BAF complex members in WT.R replicate"
	#labelPlot(dataMean1, dataMean1/dataMean2,labels,title=title,myMin=1e4) #+ scale_y_log10() + scale_x_log10()
	labelPlot(dataMean1, dataMean2,labels,title=title,myMin=1e4) #+ scale_y_log10() + scale_x_log10()
	}
	
}

arid1APlot <- function(data=ARID1A.rr, ID="O14497", labelprefix=27,legendPos=-8000,name="ARID1A"){
  acID <- which(colnames(data)=="Accession")
  a1a <- which(data[,acID]==ID)
  c1 <- grep("Abundance.F1\\.",colnames(data),value=F)
  c2 <- grep("Abundance.F2",colnames(data),value=F)
  c3 <- grep("Abundance.F3",colnames(data),value=F)
  c4 <- grep("Abundance.F4",colnames(data),value=F)
  c5 <- c(grep("Abundance.F5",colnames(data),value=F),grep("Abundance.F9",colnames(data),value=F))
  c6 <- c(grep("Abundance.F6",colnames(data),value=F),grep("Abundance.F10",colnames(data),value=F))
  cols <- c(c1,c2,c3,c4,c5,c6)
  #return(cols)
  #print(colnames(data)[cols])
 
  labels <- colnames(data)[cols]
  for(i in 1:60){
    labels[i]=substring(labels[i],labelprefix)
    if(substring(labels[i],1,1)==".")
            labels[i]=substring(labels[i],2)
    if(substring(labels[i],1,2)=="B.")
            labels[i]=substring(labels[i],3)
    labels[i]=substring(labels[i],1,nchar(labels[i])-2)
  }

  dd <- as.numeric(data[a1a,cols])

  par(mfrow=c(3,1))
  plot(dd,xaxt="n",xlab="",ylab=paste(name," abundance (raw)",sep=""),xlim=c(2.5,58.5),pch=20)
  axis(1, at=1:60, labels = F)
  text(1:60,legendPos,labels,srt=90,cex=.7,xpd=T)
  lines(c(10.5,10.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(20.5,20.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(30.5,30.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(40.5,40.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(50.5,50.5),c(-1000,1e6),lty=2,col="grey")

  mm <- max(dd[c(1,11,21,31,41,51)])
  dd[1:10] <-  dd[1:10]  * (mm/dd[1])
  dd[11:20] <- dd[11:20] * (mm/dd[11])
  dd[21:30] <- dd[21:30] * (mm/dd[21])
  dd[31:40] <- dd[31:40] * (mm/dd[31])
  dd[41:50] <- dd[41:50] * (mm/dd[41])
  dd[51:60] <- dd[51:60] * (mm/dd[51])

  #par(mar=c(6,4,1,1))
  plot(dd,xaxt="n",xlab="",ylab=paste(name," abundance",sep=""),xlim=c(2.5,58.5),pch=20)
  axis(1, at=1:60, labels = F)
  text(1:60,legendPos,labels,srt=90,cex=.7,xpd=T)
  lines(c(10.5,10.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(20.5,20.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(30.5,30.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(40.5,40.5),c(-1000,1e6),lty=2,col="grey")
  lines(c(50.5,50.5),c(-1000,1e6),lty=2,col="grey")

  dd2 <- matrix(c(dd[2:10],dd[22:29],dd[43:50],dd[12:20],dd[32:39],dd[53:60],dd[30],rep(NA,24),dd[40],rep(NA,24),dd[42],rep(NA,24),dd[52],rep(NA,24)),nrow=6,byrow=T)
  dataMean <- rep(0,25)
  dataSd <- rep(0,25)
  for(i in 1:25){
     dataMean[i] = tryCatch( mean(dd2[,i],na.rm=T),error = function(e) {0})
     dataSd[i] =   tryCatch( sd(dd2[,i],na.rm=T), error = function(e) {0})
  }
  #doAxis(17, NULL, ""); 
  plot(dataMean+dataSd,type="l",xaxt="n",yaxt="l",bty="n",lty=2,col="white",xlim=c(0,25),ylab="",xlab="",ylim=c(0,max(dataMean+dataSd)))
  errorBarPlot(25, dataMean, dataSd, "darkgrey", T, NULL, "", F, "",main=paste(name," abundance",sep=""))
  text(-.5 +1:25,legendPos,labels[c(2:10,22:29,43:50)],srt=90,cex=.7,xpd=T)
  #dd2
  
}


barPlot <- function(data1, data2 ,ymax=1.2,logV="", start = 2, legendPos=1.0, legends, labels = labels_BRG1)
{
	N <- length(labels)
	ord <- 1:N #order(labels)
	addV = T
	axesV = F
	ylimV = c(0.0,ymax)  #ylimV = c(0.01,1.2)   #ylimV = NULL
	par(mfrow=c(9,1),mar=c(0,4,1,0)+.1)
	dataSd <- data1
	dataMean <- data1
	for(i in 1:10)
	  for(j in 1:N){
            dataMean[j,i] = tryCatch( mean(c(data1[j,i],data2[j,i]),na.rm=T),error = function(e) {0}) 
	    dataSd[j,i] =   tryCatch( sd(c(data1[j,i],data2[j,i]),na.rm=T), error = function(e) {0})
	  }
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+0],dataSd[ord,start+0],"darkgrey",addV,ylimV,logV,axesV,legends[start+0])
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+1],dataSd[ord,start+1],"red",addV,ylimV,logV,axesV,legends[start+1])
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+2],dataSd[ord,start+2],"orange",addV,ylimV,logV,axesV,legends[start+2]) 
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+3],dataSd[ord,start+3],"green",addV,ylimV,logV,axesV,legends[start+3])
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+4],dataSd[ord,start+4],"lightblue",addV,ylimV,logV,axesV,legends[start+4])
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+5],dataSd[ord,start+5],"blue",addV,ylimV,logV,axesV,legends[start+5])
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+6],dataSd[ord,start+6],"brown",addV,ylimV,logV,axesV,legends[start+6])
	doAxis(N, ylimV, logV); errorBarPlot(N, dataMean[ord,start+7],dataSd[ord,start+7],"purple",addV,ylimV,logV,axesV,legends[start+7])
	barplot(data1[ord,start+7],space=0,axes=F,col="white",border=NA)
	text(1:N - .5, legendPos, srt = 90, adj = 1,labels = labels[ord], xpd = TRUE)
}


labels_BRG1_full <- c("ARID1A","SS18L1","ACTL6A","SMARCA2","SMARCA4","ACTB","SMARCB1","SS18","BCL7A","ARID2","SMARCD3","PBRM1","ARID1B","SMARCC2","PHF10","BCL7C","DPF1","DPF3","DPF2","SMARCC1","SMARCD2","SMARCE1","SMARCD1","BCL7B","BCL11A","BRD9","BRD7")

exclude_list_BRG1 <- c("P60709","P51531","Q9H165","Q9C0K0")#"P51532")

labels_BRG1 <- c("ARID1A","SS18L1","ACTL6A","SMARCA4","SMARCB1","SS18","BCL7A","ARID2","SMARCD3","PBRM1","ARID1B","SMARCC2","PHF10","BCL7C","DPF1","DPF3","DPF2","SMARCC1","SMARCD2","SMARCE1","SMARCD1","BCL7B","BRD9","BRD7")

labels_A1A_full <- c("SMARCD3","SMARCA4","SMARCC2","PHF10","BRD9","SMARCB1","ARID1B","SMARCD1","PBRM1","SMARCC1","SMARCE1","SS18","BCL11B","SMARCD2","BCL7C","DPF1","BCL11A","ACTB","DPF2","BRD7","ARID2","DPF3","ARID1A","ACTL6A","SMARCA2","BCL7B","SS18L1","ACTL6B","BCL7A")

exclude_list_ARID1A <- c("Q8WUB8","Q86U86","Q9NPI1","Q9H8M2","Q9H165","Q68CP9","Q8NFD5","P60709", "O94805","Q9C0K0")#"O14497")

labels_A1A <- c("SMARCD3","SMARCA4","SMARCC2","SMARCB1","SMARCD1","SMARCC1","SMARCE1","SS18","BCL11B","SMARCD2","BCL7C","DPF1","DPF2","DPF3","ARID1A","ACTL6A","SMARCA2","BCL7B","SS18L1","ACTL6B","BCL7A")

rna_proteins <- c('P05387','P62917','Q9P015','P62888','Q14684-1','P62241','P18077','P62906','P27635','P39023','Q02543','P18124','Q6DKI1-1','P30050-1','Q9Y676','P42677','P62899-1','P62913-1','P49207','P15880','P18621-1','Q9Y3U8','P46777','P05388-1','Q07020-1','Q9Y4W2-1','O60287','P62424','P62277','P61254','P62829','P62244','P46776','P40429','P32969','Q02878','Q9Y3B7-1','P83731','P62750','P46781','P84098','P0CW22,','P56182','P62280','O43159','O76021-1','P61247','P50914','P62854','P62701','P46778','P46779-1','P61353','P25398','P35268','Q9Y3A4','P61513','P46782','P62851','Q92552-1','P62263','Q96EU6-1','P62753','P61313-1','P62249','P82650-1','P36578','Q9BYN8','P62847-4','P23396-1','P63173','Q9Y3D9','P83881','P42766','P62269','P62273-1','P62910','P62979','P47914','P60866-1','P62857','P62266','P39019','P82933','P26373-1','P08865','P46783','P62081')

rna_proteins <- c("P07305", "Q8IZA3", "Q92522", "P0C5Y9", "P0C5Z0", "H0YFX9", "Q9BTM1", "A8MQC5", "C9J0D1", "C9J386", "E5RJU1", 
		  "Q71UI9", "P16104", "B4DJC3", "D6RCF2", "O75367", "Q5SQT3", "Q9P0M6", "P0C0S5", "P0C1H6", "A9UJN3", "P57053",
                  "Q7Z2G1", "B4DEB1", "P84243", "B2R4P9", "K7EMV3", "K7ES00", "K7EK07", "K7EP01", "Q6NXT2", "Q02539", "P16401",
                  "P16403", "P16402", "Q4VB24", "P10412", "A3R0T8", "A1L407", "P22492", "Q96QV6", "P04908", "Q08AJ9", "Q93077",
                  "P20671", "P0C0S8", "A3KPC7", "Q96KK5", "Q99878", "A4FTV9", "Q92646", "Q96A08", "P33778", "P62807", "P58876",
                  "B2R4S9", "Q93079", "P06899", "O60814", "Q99880", "I6L9F7", "Q99879", "Q99877", "P23527", "P68431", "P62805",
                  "Q99525", "Q0VAS5", "B2R4R0", "Q6FI13", "Q8IUE6", "Q16777", "Q16778", "B4DR52", "Q5QNW6", "Q71DI3", "Q5TEC6",
                  "Q7L7L0", "Q8N257", "Q16695", "Q6TXQ4", "Q14463", "B4E0B3", "B2R5B6", "A2RUA4", "B2R5B3", "Q9HA11", "A8K9J7",
                  "B2R6Y1", "B4E380", "A8K4Y7", "Q6B823", "Q6LBZ2", "A3R0T7")

getAllRatios <- function(replicate="A",ratioData=BRG1.rr,abundanceData=BRG1.ab, exclude_list=c(), capAt2=T){
    pattern <- paste("Abundance.Ratios.by.Bio.Rep.*",replicate,"$",sep="")
    columns <- grep(pattern, colnames(ratioData),perl=T)
    columns2 <- grep("Abundance.F", colnames(abundanceData),perl=T)

    F1 <- ratioData[which(!(ratioData[,"Accession"] %in% exclude_list)),columns]
    M <- dim(F1)[1]
    nC <- dim(F1)[2]
    for(i in 1:M){
            if (capAt2) F1[i,which(F1[i,]==0.01)]=NA
            if (capAt2) F1[i,which(F1[i,]>2)]=2.0
    }
    for(i in 1:nC){
            colnames(F1)[i] <- substring(colnames(F1)[i],29,nchar(colnames(F1)[i])-7)
    }

    F2 <- abundanceData[which(!(ratioData[,"Accession"] %in% exclude_list)),columns2]

    cbind(F1,F2)
}

getAllRatios2 <- function(replicate="A",ratioData=BRG1.rr,abundanceData=BRG1.ab, exclude_list=c(), capAt2=T){
    columns2 <- grep("Abundance.F", colnames(abundanceData),perl=T)

    F1 <- getRatios2(replicate=replicate,ratioData=ratioData, exclude_list=exclude_list, full=TRUE, capAt2 = capAt2)

    columns2 <- grep("Abundance.F", colnames(abundanceData),perl=T)
    F2 <- abundanceData[which(!(ratioData[,"Accession"] %in% exclude_list)),columns2]

    cbind(F1,F2)
}


testMe <- function(data.rr, data.ab, thresh = 5e07, cor_threshold= 0.5, name, validV = c(1:7,10:11,15:30), maxMissing=7, minMissing=-1){
    rARID1A.A <- getAllRatios2(replicate="A",ratioData=data.rr,abundanceData=data.ab,capAt2=F)
    rARID1A.B <- getAllRatios2(replicate="B",ratioData=data.rr,abundanceData=data.ab,capAt2=F)
    N <- dim(rARID1A.A)[1]
    
    val <- rep(0,N)
    for (i in 1:N){
    	#val[i] <- tryCatch(0.01 + sqrt( mean(((as.numeric(rARID1A.A[i,validV]) - as.numeric(rARID1A.B[i,validV]))^2),na.rm=T )) , error=function(e) {NA})
	val[i] <- 0.1 + 3* tryCatch(1.0-cor(as.numeric(rARID1A.A[i,validV]),as.numeric(rARID1A.B[i,validV]),use="complete.obs"), error=function(e) {NA})
      if (sum(is.na(as.numeric(rARID1A.A[i,validV]) + as.numeric(rARID1A.B[i,validV])))>maxMissing && data.rr[i,"Marked.as"]!="BAF_complex") val[i] = NA
      #if (sum(is.na(as.numeric(rARID1A.A[i,validV]) + as.numeric(rARID1A.B[i,validV])))<=minMissing && data.rr[i,"Marked.as"]!="BAF_complex") val[i] = NA
    }

    rARID1A.A <- getAllRatios2(replicate="A",ratioData=data.rr,abundanceData=data.ab,capAt2=T)
    rARID1A.B <- getAllRatios2(replicate="B",ratioData=data.rr,abundanceData=data.ab,capAt2=T)
    N <- dim(rARID1A.A)[1]

    val2 <- rep(0,N)
    for (i in 1:N){
        val2[i] <- tryCatch( max(0.1, sqrt( mean(((as.numeric(rARID1A.A[i,validV]) - as.numeric(rARID1A.B[i,validV]))^2),na.rm=T) ) ) , error=function(e) {NA})
        #val[i] <- tryCatch(1-cor(as.numeric(rARID1A.A[i,validV]),as.numeric(rARID1A.B[i,validV]),use="complete.obs"), error=function(e) {NA})
      if (sum(is.na(as.numeric(rARID1A.A[i,validV]) + as.numeric(rARID1A.B[i,validV])))>maxMissing && data.rr[i,"Marked.as"]!="BAF_complex") val2[i] = NA
      #if (sum(is.na(as.numeric(rARID1A.A[i,validV]) + as.numeric(rARID1A.B[i,validV])))<=minMissing && data.rr[i,"Marked.as"]!="BAF_complex") val2[i] = NA
      val[i] <- min(val[i],val2[i],na.rm=T)
      if (is.infinite(val[i])) val[i] = NA
      #val[i] <- min(1.3,val[i])
      #if(i==5484) print(c(rARID1A.A[i,validV],rARID1A.B[i,validV]))
    }


    labels <- rep("",N)
    for (i in 1:N){
      descr <- data.rr[i,"Description"]
      res <- regexpr("GN=[^ ]*",descr)
      labels[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
      if (labels[i]=="BRD9") print(c(i,mm[i],val[i]))
    }

    size <- rep(0.3,N)
    for (i in 1:N){
	ws <- which(rabbitControl[,1]==data.rr[i,"Accession"])
    	if (length(ws>0)){
	  if (name=="BRG1" && rabbitControl[ws[1],2]+rabbitControl[ws[1],3]<(rabbitControl[ws[1],4]+rabbitControl[ws[1],5])/5){
		size[i] = 1;
	  }else if (name!="BRG1" && rabbitControl[ws[1],2]+rabbitControl[ws[1],3]<(rabbitControl[ws[1],6]+rabbitControl[ws[1],7])/5)
		  size[i] = 1;
	}

    }

    mm <- rep(0,N)
    for(i in 1:N){
	    #print(i)
	    theoret.peptides <- length(observable.peptides(seq=gsub("X","A",data.rr[i,"Sequence"]),nmc=1)[[1]]) 
	    mm[i] <- tryCatch( mean(as.numeric(rARID1A.A[i,31:36]),na.rm=T),  error=function(e) {NA})/theoret.peptides
    }

    w <- which(data.rr[,"Master"]=="IsMasterProtein" & data.rr[,"Protein.FDR.Confidence.Combined"]=="High" & labels!="")

    #plot(val[w],mm[w],log="xy",col="grey",pch=20,main=name,ylab="Abundance",xlab="difference between bio-replicates")
    #return()

    w2 <- which((size == 1 & mm>thresh & val<cor_threshold & 
		data.rr[,"Master"]=="IsMasterProtein" & 
		data.rr[,"Protein.FDR.Confidence.Combined"]=="High" & labels!="") | data.rr[,"Marked.as"]=="BAF_complex" )

    #text(val[w2],mm[w2],labels[w2],cex=0.8,col="black")
    
    color <- rep("not observed in rabbitIGG control",N)
    color[data.rr[,"Marked.as"]=="BAF_complex"]="known BAFcomplex members"

    ss <- labelPlotBackground(val[w],mm[w],val[w2],mm[w2],size[w2],color[w2] ,labels[w2],title=name, logX=F, logY=T)

    # add red BRG1 known complex members
    #w3 <- which(data.rr[,"Marked.as"]=="BAF_complex")
    #dd <- cbind(val[w3], mm[w3]) 
    #colnames(dd) <- c("x","y");    dd <- as.data.frame(dd)
    #ss <- ss + geom_point(data = dd, aes(x = x, y = y), colour = "red", size = 1) +
    #	  geom_text_repel(data = dd, aes(x = x, y = y), label=labels[w3], colour = "red") 
    
    #lines(val[which(data.rr[,"Marked.as"]=="BAF_complex")],mm[which(data.rr[,"Marked.as"]=="BAF_complex")],col="red",type="p")
    print(ss)
    
    #print(data.rr[w2,"Description"])
    return(labels[w2])
}


labelPlot <- function(x,y,labels,title="",myMin=0){
  #x[which(x>3)]=3
  #y[which(y>3)]=3
  x[which(x==0)]=0.1
  y[which(y==0)]=0.1
  mm <- rep(0,length(x))
  for(i in 1:length(x)) mm[i] <- max(x[i],y[i])
  MM <- max(mm,na.rm=T)
  dd <- cbind(x,y); colnames(dd) <- c("replicate1","replicate2")
  dd <- as.data.frame(dd)
  ggplot(dd,aes(x = replicate1, y = replicate2)) + geom_point(colour="black") +
  geom_text_repel(data = dd, aes(label=labels)) + 
  ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) + 
  #coord_cartesian(xlim = c(0., MM), ylim = c(myMin, MM))  
  #labs(x = "Abundance in ARID1A pulldown", y= "Ratio of abundance in ARID1A pulldown vs. BRG1 pulldown") +
  #scale_y_log10(limits=c(0.005,200)) + scale_x_log10(limits=c(1e5,1e10)) #+ geom_abline(intercept = 0.03887, slope = 0.98279)
  labs(x = "Abundance in ARID1A pulldown", y= "Abundance in BRG1 pulldown") +
  scale_y_log10(limits=c(1e3,2e8)) + scale_x_log10(limits=c(1e3,2e8)) + geom_abline(intercept = 0.03887, slope = 0.98279)
  

  #geom_point(colour="grey") +
  #geom_point(data=dd, aes(x = replicate1, y = replicate2)) +
}

labelPlotBackground <- function(xb,yb,x,y,size, color,labels,title="", logX=F, logY=F){
dd <- cbind(xb,yb); colnames(dd) <- c("interaction_signal","mean_MS1_abundance")
  dd <- as.data.frame(dd)

  dd2 <- cbind(x,y,size); colnames(dd2) <- c("interaction_signal","mean_MS1_abundance","size")
  dd2 <- as.data.frame(dd2)
  dd2$Protein_type <- color

  ss <- ggplot(dd,aes(x = interaction_signal, y = mean_MS1_abundance, colour= color)) + geom_point(colour="grey") +
  geom_point(data = dd2, aes(x = interaction_signal, y = mean_MS1_abundance, colour = color)) +
  #scale_size(range=c(1,1)) +
  geom_text_repel(data = dd2, aes(label=labels, colour = color)) +
  ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
  { if (logX) scale_x_log10() } + { if (logX==F) scale_x_continuous() } +
  { if (logY) scale_y_log10(limits=c(1e4,1e9)) } + { if (logY==F) scale_y_continuous() }
  return(ss)
}

plot_BRG1_peptides_all <- function(filename="BRG1_peptides_details.pdf"){
    pdf(filename,width=8,height=38)
    par(mfrow=c(34,1),mar=c(2, 1, 1.5, 1) +0.1)
    for(i in 66:99)
	    plot_BRG1_peptides(column=i)
    dev.off()
}



plot_BRG1_peptides <- function(protein="PBRM1",KO="ARID2A",column=NA,threshold=100){
   if (is.na(column))
	column <- 65 + grep(KO,colnames(BRG.pep)[66:99])[1] 
   id <- BRG1.rr[which(BRG1.lab==protein)[1],4]
   seq <- ARID1A.rr[which(ARID1A.rr[,6]==id),8]
   L <- nchar(seq)
   w <- which(BRG.pep[,11]==id)
   N <- length(w)
   
   dd <- matrix(NA,nrow=N,ncol=4)
   for(i in 1:N){
      pos <- BRG.pep[w[i],12]
      start <- regexpr("\\[",pos)[1]
      end <- regexpr("-",pos)[1]
      end2 <- regexpr("\\]",pos)[1]
      #print(pos)
      s1 <- strtoi(substring(pos,start+1,end-1), base = 0L)
      e1 <- strtoi(substring(pos,end+1,end2-1), base = 0L)
      #print(c(s1,e1))
      #lines(c(s1,e1),c(i,i))
      if(regexpr("Phospho",BRG.pep[w[i],4])[1]!=-2 && !is.na(BRG.pep[w[i],192]) && BRG.pep[w[i],192]>threshold) #&& !is.na(BRG.pep[w[i],column]))
      {
        dd[i,1]=s1
        dd[i,2]=e1
	dd[i,3]=round(min(BRG.pep[w[i],column],2)/2.0 * 99) + 1
	if(regexpr("Phospho",BRG.pep[w[i],4])[1]!=-1)
		dd[i,4]=30
	else
		dd[i,4]=NA
	if(is.na(dd[i,3])) dd[i,3]=1#101
      }
   }

   maxpos=-1
   oo <- order(dd[,1])
   #print(dd[1:10,3])
   oldE = -1
   oldS = -1
   for(i in 1:N){
      if(!is.na(dd[oo[i],1])){
      #print(i)
      start <- dd[oo[i],1]
      end <- dd[oo[i],2]
      if(oldE<start)
	   pos = 1
      else
	   pos = pos + 1
      #pos = i
      #rect(dd[oo[i],1],pos,dd[oo[i],2],pos+1,col=color3[dd[oo[i],3]],lw=2,density=dd[oo[i],4],border="black")
      if(pos>maxpos) maxpos=pos
      oldS <- start
      oldE <- end
      }
   }

   name <- substring(colnames(BRG.pep)[column],29,nchar(colnames(BRG.pep)[column])-7)
   plot(0,0,xlim=c(1,L),ylim=c(1,maxpos+1),col="white",main=paste(name," ratios of ",protein," peptides.",sep=""),yaxt="n")
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="grey")

   oldE = -1
   oldS = -1
   for(i in 1:N){
      if(!is.na(dd[oo[i],1])){
      #print(i)
      start <- dd[oo[i],1]
      end <- dd[oo[i],2]
      if(oldE<start)
           pos = 1
      else
           pos = pos + 1
      #pos = i
      rect(dd[oo[i],1],pos,dd[oo[i],2],pos+1,col=color3[dd[oo[i],3]],lw=2.,density=dd[oo[i],4],border="black")
      oldS <- start
      oldE <- end
      }
   }

}

otherProteins <- function(){	
       
       columns = c(1:7,10:11,15:30)
       pdf("otherProteins-ARID1A_1.pdf",width=17,height=14)
         proteins <- testMe(ARID1A.rr2,ARID1A.ab2,thresh=1e0,cor_threshold=0.5,name="ARID1A", validV=columns) 
       dev.off()
       otherProteinsDetails(proteins,filename="otherProteins-ARID1A-list.pdf",data=ARID1A.rr2,columns=columns,name="ARID1A")
       otherProteinsDetails(NA,filename="otherProteins-ARID1A-SandraList.pdf",data=ARID1A.rr2,columns=columns,name="ARID1A")
       
       
       columns = c(1:4,7:11,15:30)
       pdf("otherProteins-BRG1_1.pdf",width=17,height=14)
         proteins <- testMe(BRG1.rr2,BRG1.ab2,thresh=1e0, name = "BRG1", cor_threshold=0.5, validV=columns)
       dev.off()
       return(1)
       otherProteinsDetails(proteins,filename="otherProteins-BRG1-list.pdf",data=BRG1.rr2,columns=columns,name="BRG1")
       otherProteinsDetails(NA,filename="otherProteins-BRG1-SandraList.pdf",data=BRG1.rr2,columns=columns,name="BRG1")
}

otherProteinsDetails <- function(proteins=NA,filename="otherProteins-BRG1-2.pdf", data=BRG1.rr2,columns=c(1:4,7:11,15:30),name="BRG1")
{
  pdf(filename,width=7,height=7)

    rBRG1.A <- getRatios2(replicate="A",ratioData=data,capAt2=F, full=T, exclude_list=c())
    rBRG1.B <- getRatios2(replicate="B",ratioData=data,capAt2=F, full=T, exclude_list=c())
    if (is.na(proteins)){
    	proteins <- c("BRD2","BRD3","BRD4","BRD8","BICRA","NSF","REST","BICRAL","CHD4","GRN","CDK2AP2","FOSL1","SPICE1","TBPL2","NUP50","KEAP1","KDM4A","SS18L2","LDB2","CBX1","MBD3","NUP50","TP53BP1","GATAD2B","TAF6","TAF5","TAF4","TAF10","KPNA4","ADNP","SSBP3","SSBP4","KLF5","KLF4","EBF3","AK7","MED14","CECR2","NUP153","SALL3","JUNB","MTA2","MTA3","C1QBP","SCAI","GRB2","KPNA6","CTNNB1","HDAC2","HDAC1","CRTAP","MYH9","MYL6","NES","NEXN","LHX2","ZBTB33")
    }
    N <- length(proteins)

    M <- dim(data)[1]
    labels <- rep("",M)
    for (i in 1:M){
      descr <- data[i,"Description"]
      res <- regexpr("GN=[^ ]*",descr)
      labels[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
    }

    for(j in 1:N){
      i <- which(labels==proteins[j])[1];
      #rBRG1.A[i,which(is.na(rBRG1.A[i,]))]=0.0
      #rBRG1.B[i,which(is.na(rBRG1.B[i,]))]=0.0
      print(labelPlot(as.numeric(rBRG1.A[i,columns]),as.numeric(rBRG1.B[i,columns]),colnames(rBRG1.A)[columns],paste(labels[i],"in",name,"pulldowns")))
    }
  dev.off()
}

# returns data about WT and WT.R relative abundances (for given columns)
getWT_WTR <- function(i1,i2,i3,i4, minPep=1){
  w3 <- which(BRG1.rr2[,3]=="IsMasterProtein")
  acs <- BRG1.rr2[w3,4]
  N <- length(acs)
  dd <- matrix(NA,nrow=N,ncol=5)
  i <- 1
  #i1 <- 331
  #i2 <- 353
  #i3 <- 332
  #i4 <- 354
  for(aa in acs){
    w2 <- which(BRG.pep2[,"Marked.as"]=="BAF_complex" & BRG.pep2[,"Master.Protein.Accessions"]==aa & BRG.pep2[,"Quan.Info"]=="" & !is.na(BRG.pep2[,i1]) & !is.na(BRG.pep2[,i2]) & !is.na(BRG.pep2[,i3]) & !is.na(BRG.pep2[,i4]))
    if(length(w2)>=minPep){
      if(aa=="P60709"){   # actin
	 print(BRG.pep2[w2,i1])
    	 print(BRG.pep2[w2,i2])
	 print(BRG.pep2[w2,i3])
	 print(BRG.pep2[w2,i4])
      }
      else{
        dd[i,1] <- sum(BRG.pep2[w2,i1],na.rm=T)
        dd[i,2] <- sum(BRG.pep2[w2,i2],na.rm=T)
        dd[i,3] <- sum(BRG.pep2[w2,i3],na.rm=T)
        dd[i,4] <- sum(BRG.pep2[w2,i4],na.rm=T)
        dd[i,5] <- length(w2)
	#print(aa)
      }
    }
    i <- i + 1
  }
  dd
}


get_peptide_WT_WTR <- function(i1,i2,i3,i4){
  w3 <- which(BRG1.rr2[,3]=="IsMasterProtein")
  acs <- BRG1.rr2[w3,4]
  N <- length(acs)
  dd <- matrix(NA,nrow=N,ncol=5)
  i <- 1
  #i1 <- 331
  #i2 <- 353
  #i3 <- 332
  #i4 <- 354
  xx <- c()
  for(aa in acs){
      if(aa!="P60709"){  # no actin
        w2 <- which(BRG.pep2[,"Marked.as"]=="BAF_complex" & BRG.pep2[,"Master.Protein.Accessions"]==aa & BRG.pep2[,"Quan.Info"]=="" & !is.na(BRG.pep2[,i1]) & !is.na(BRG.pep2[,i2]) & !is.na(BRG.pep2[,i3]) & !is.na(BRG.pep2[,i4]))
        L <- length(w2)
	if (L>0){
	  for(j in 1:L)
            xx <- c(xx,BRG.pep2[w2[j],c(i1,i2,i3,i4)])
	    #break;
	}
      }
    i <- i + 1
  }
  dd <- matrix(xx,ncol=4,byrow=T)
  dd
  #xx
}



getPeptideRatios <- function(protein,replicate="A",ratioData=BRG.pep, proteinLabels=BRG1.lab, protein.accesions=BRG1.rr[,"Accession"]){
    pattern <- paste("Abundance.Ratios.by.Bio.Rep.*",replicate,"$",sep="")
    columns <- grep(pattern, colnames(ratioData),perl=T)

    id <- protein.accesions[which(proteinLabels==protein)[1]]
    seq <- ARID1A.rr[which(ARID1A.rr[,"Accession"]==id),"Sequence"]
    L <- nchar(seq)
    w <- which(ratioData[,"Master.Protein.Accessions"]==id)
    M <- length(w)


    F1 <- ratioData[w,columns]
    nC <- dim(F1)[2]

    for(i in 1:M){
            #F1[i,which(F1[i,]==0.01)]=NA

            F1[i,which(F1[i,]>2)]=2.0
    }

    for(i in 1:nC){
            colnames(F1)[i] <- substring(colnames(F1)[i],29,nchar(colnames(F1)[i])-7)
    }

    cc <- rep(0,M)
    for(i in 1:M){
      pos <- ratioData[w[i],"Positions.in.Master.Proteins"]
      start <- regexpr("\\[",pos)[1]
      end <- regexpr("-",pos)[1]
      end2 <- regexpr("\\]",pos)[1]
      s1 <- strtoi(substring(pos,start+1,end-1), base = 0L)
      e1 <- strtoi(substring(pos,end+1,end2-1), base = 0L)
      cc[i] <- s1
      name <- paste(s1,"-",e1,":",ratioData[w[i],"Modifications"],sep="")#":",ratioData[w[i],5],sep="")
      name <- gsub("Carbamidomethyl", "",name)
      name <- gsub("Oxidation", "",name)
      name <- gsub("\\[", "", name)
      name <- gsub("\\]", "", name)
      name <- gsub("N-Term", "Nterm", name)
      name <- gsub("TMT6plex", "TMT", name)
      rownames(F1)[i] <- name
    }

    oo <- order(cc)
	
    F1[oo,]
}


getPeptideRatios2 <- function(protein,replicate="A",ratioData=BRG.pep3, proteinLabels=BRG1.rr2[,"Description"], protein.accesions=BRG1.rr2[,"Accession"]){
    
    # preprocessing  ========================================================================================

    labels <- rep("",length(proteinLabels))
    for (i in 1:N){
      descr <- proteinLabels[i]
      res <- regexpr("GN=[^ ]*",descr)
      labels[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
    }
    
	
    id <- protein.accesions[which(labels==protein)[1]]
    myW <- which(ratioData[,"Master.Protein.Accessions"]==id & ratioData[,"Modifications.in.Master.Proteins"]!="")

    #print(c(protein,length(myW)))
    if (length(myW)==0) return (NULL);
    #=======================================================================================================
    # procesing body
	
    WR <- grep(paste("Abundances.Normalized.*\\.",replicate,"\\.WT\\.R",sep=""), colnames(ratioData),perl=T)
    wrN <- length(WR)

    # get fileID of each WT.R replicate (eg. "F1")
    WRfileID <- rep("", wrN)
    for (j in 1:wrN){
        WRfileID[j] <- substring(colnames(ratioData)[WR[j]],23,24) # eg WRfileID[j]="F1"
    }

    pattern <- paste("Abundances.Normalized.*\\.",replicate,"\\.",sep="")
    columns <- grep(pattern, colnames(ratioData),perl=T)
    columnsDivide <- columns
    nC <- length(columns)

    # get file ID of each replicate and divide by corresponding WR replicate
    for (i in 1:nC){
          fileID <- substring(colnames(ratioData)[columns[i]],23,24)
          columnsDivide[i] <- WR[which(WRfileID==fileID)]
    }

    #print(columns)
    #print(columnsDivide)

    F1 <- ratioData[myW,columns]
    M <- dim(F1)[1]
    nC <- dim(F1)[2]

    #print(dim(F1))
    for(i in 1:nC){
            F1[,i] <- F1[,i] / ratioData[myW,columnsDivide[i]]
            if(substring(colnames(F1)[i],39,39)==".") { start=40; } else { start=39; }
            name <- substring(colnames(F1)[i],start,nchar(colnames(F1)[i])-2)
            if (name=="WT")
                    name <- substring(colnames(F1)[i],start,nchar(colnames(F1)[i]))
            if (substring(name,2,2)==".")
                    name <- substring(name,3,nchar(name))
            colnames(F1)[i] <- paste(name,replicate,sep="-")
    }

    # postprocessing =================================================================

    for(i in 1:M){
	    F1[i,which(F1[i,]<=0.015)] = 0.02
            #F1[i,which(is.na(F1[i,]))]=0.015
	    F1[i,which(is.na(F1[i,]))]=0.0
            F1[i,which(F1[i,]>2)]=2.0
    }

    #print(dim(F1))
    cc <- rep(0,M)
    for(i in 1:M){
      pos <- ratioData[myW[i],"Positions.in.Master.Proteins"]
      start <- regexpr("\\[",pos)[1]
      end <- regexpr("-",pos)[1]
      end2 <- regexpr("\\]",pos)[1]
      s1 <- strtoi(substring(pos,start+1,end-1), base = 0L)
      e1 <- strtoi(substring(pos,end+1,end2-1), base = 0L)
      cc[i] <- s1
      name <- paste(s1,"-",e1,":",ratioData[myW[i],"Modifications"],sep="")#":",ratioData[myW[i],5],sep="")
      name <- gsub("Carbamidomethyl", "",name)
      name <- gsub("Oxidation", "",name)
      name <- gsub("\\[", "", name)
      name <- gsub("\\]", "", name)
      name <- gsub("N-Term", "Nterm", name)
      name <- gsub("TMT6plex", "TMT", name)
      rownames(F1)[i] <- name
    }

    oo <- order(cc)

    F1[oo,]
}

getRatios <- function(replicate="A",ratioData=BRG1.rr, exclude_list){
    pattern <- paste("Abundance.Ratios.by.Bio.Rep.*",replicate,"$",sep="")
    columns <- grep(pattern, colnames(ratioData),perl=T)

    w <- which(ratioData[,"Marked.as"]=="BAF_complex" & !(ratioData[,"Accession"] %in% exclude_list))
    F1 <- ratioData[w,columns]
    M <- dim(F1)[1]
    nC <- dim(F1)[2]
    for(i in 1:M){
	    #F1[i,which(F1[i,]==0.01)]=NA
	    F1[i,which(F1[i,]>2)]=2.0
    }

    for(i in 1:nC){
	    colnames(F1)[i] <- substring(colnames(F1)[i],29,nchar(colnames(F1)[i])-7)
    }
    
    descriptionID <-  which(colnames(ratioData)=="Description")
    for(i in 1: M){
       descr <- ratioData[w[i],descriptionID]
       res <- regexpr("GN=[^ ]*",descr)
       rownames(F1)[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
       #if (rownames(F1)[i]=="NA") rownames(F1)[i] <- descr
    }
    
    F1
}

getRatios2 <- function(replicate="A",ratioData=BRG1.rr, exclude_list, full=FALSE, capAt2 = T){
    WR <- grep(paste("Abundances.Normalized.*\\.",replicate,"\\.WT\\.R",sep=""), colnames(ratioData),perl=T)
    wrN <- length(WR)

    # get fileID of each WT.R replicate (eg. "F1")
    WRfileID <- rep("", wrN)
    for (j in 1:wrN){
	WRfileID[j] <- substring(colnames(ratioData)[WR[j]],23,24)
    }

    pattern <- paste("Abundances.Normalized.*\\.",replicate,"\\.",sep="")
    columns <- grep(pattern, colnames(ratioData),perl=T)
    columnsDivide <- columns
    nC <- length(columns)
    
    # get file ID of each replicate and divide by corresponding WR replicate
    for (i in 1:nC){
	  fileID <- substring(colnames(ratioData)[columns[i]],23,24)
	  columnsDivide[i] <- WR[which(WRfileID==fileID)]
    }

    #print(columns)
    #print(columnsDivide)

    w <- which(ratioData[,"Marked.as"]=="BAF_complex" & !(ratioData[,"Accession"] %in% exclude_list))
    if (full==T) 
	    w <- 1:nrow(ratioData)

    F1 <- ratioData[w,columns]
    M <- dim(F1)[1]

    for(i in 1:nC){
	    F1[,i] <- F1[,i] / ratioData[w,columnsDivide[i]]
	    if(substring(colnames(F1)[i],39,39)==".") { start=40; } else { start=39; } 
            name <- substring(colnames(F1)[i],start,nchar(colnames(F1)[i])-2)
	    if (name=="WT") 
		    name <- substring(colnames(F1)[i],start,nchar(colnames(F1)[i]))
	    if (substring(name,2,2)==".") 
		    name <- substring(name,3,nchar(name))
	    colnames(F1)[i] <- name
    }

    for(i in 1:M){
            #F1[i,which(F1[i,]==0.01)]=NA
            if (capAt2) F1[i,which(F1[i,]>2)]=2.0
    }

    for(i in 1: M){
       descr <- ratioData[w[i],"Description"]
       res <- regexpr("GN=[^ ]*",descr)
       if (full)
	       rownames(F1)[i] <- paste(i,substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1))
       else
	       rownames(F1)[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
       #if (rownames(F1)[i]=="NA") rownames(F1)[i] <- descr
    }

    F1
}

getGroupedRatios <- function(ratioData=BRG1.rr2, exclude_list){
  pattern <- "Abundance.Ratio."
  columns <- grep(pattern, colnames(ratioData),perl=T)[1:25]

  w <- which(ratioData[,"Marked.as"]=="BAF_complex" & !(ratioData[,"Accession"] %in% exclude_list))
  F1 <- ratioData[w,columns]

  M <- dim(F1)[1]
    nC <- dim(F1)[2]
    for(i in 1:M){
            #F1[i,which(F1[i,]==0.01)]=NA
            F1[i,which(F1[i,]>2)]=2.0
    }

    for(i in 1:nC){
       colnames(F1)[i] <- substring(colnames(F1)[i],17,nchar(colnames(F1)[i])-6)
    }

    descriptionID <-  which(colnames(ratioData)=="Description")
    for(i in 1: M){
       descr <- ratioData[w[i],descriptionID]
       res <- regexpr("GN=[^ ]*",descr)
       rownames(F1)[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
    }

    F1
}

getGroupedRatios2 <- function(ratioData=BRG1.rr2, exclude_list){
   rBRG1.A <- getRatios2(ratioData=ratioData,replicate="A",exclude_list=exclude_list)
   rBRG1.B <- getRatios2(ratioData=ratioData,replicate="B",exclude_list=exclude_list)
   return (sqrt(rBRG1.A * rBRG1.B))
}

getRatiosPvalue <- function(ratioData=BRG1.rr2, exclude_list){
  pattern <- "Abundance.Ratio.P.Value"
  columns <- grep(pattern, colnames(ratioData),perl=T)

  w <- which(ratioData[,"Marked.as"]=="BAF_complex" & !(ratioData[,"Accession"] %in% exclude_list))
  F1 <- ratioData[w,columns]

  M <- dim(F1)[1]
    nC <- dim(F1)[2]
    for(i in 1:M){
            #F1[i,which(F1[i,]==0.01)]=NA
            F1[i,which(F1[i,]>2)]=2.0
    }

    for(i in 1:nC){
       colnames(F1)[i] <- substring(colnames(F1)[i],25,nchar(colnames(F1)[i])-6)
    }

    descriptionID <-  which(colnames(ratioData)=="Description")
    for(i in 1: M){
       descr <- ratioData[w[i],descriptionID]
       res <- regexpr("GN=[^ ]*",descr)
       rownames(F1)[i] <- substring(descr,res[1]+3,res[1]+attr(res,"match.length")-1)
    }

    F1
}

getPvalue <- function(ratio1,ratio2,mean=1.0,sd=0.3936362){
    # oposite directions, no significant
    if (log(ratio1/mean)*log(ratio2/mean)<0) return(1.0)
    

    #else
    P1 <- pnorm(log(ratio1),mean=log(mean), sd=sd,lower.tail=F)
    P2 <- pnorm(log(ratio1),mean=log(mean), sd=sd,lower.tail=T)
    PA <- min(P1,P2)

    P1 <- pnorm(log(ratio2),mean=log(mean), sd=sd,lower.tail=F)
    P2 <- pnorm(log(ratio2),mean=log(mean), sd=sd,lower.tail=T)
    PB <- min(P1,P2)

    return(pchisq( -2*sum(log(c(PA,PB))), 4, lower.tail=FALSE))
}


getRatiosPvalue2 <- function(ratioData=BRG1.rr2, exclude_list, columns){
   rBRG1.A <- getRatios2(ratioData=ratioData,replicate="A",exclude_list=exclude_list, full=T)
   rBRG1.B <- getRatios2(ratioData=ratioData,replicate="B",exclude_list=exclude_list, full=T)

   #WT_variability <- 0.1034426
   WT_variability <- 0.09546746

   P1 <- matrix(NA,nrow=nrow(rBRG1.A),ncol=ncol(rBRG1.A))
   P2 <- P1;   PA <- P1
   for(i in 1:dim(P1)[1]){
   	P1[i,] <- pnorm(as.numeric(log(rBRG1.A[i,])),mean = 0, sd = WT_variability, lower.tail=F)
        P2[i,] <- pnorm(as.numeric(log(rBRG1.A[i,])),mean = 0, sd = WT_variability, lower.tail=T)
	for(j in 1:dim(P1)[2])
		PA[i,j] <- min(P1[i,j], P2[i,j])	
   }


   PB <- P1
   for(i in 1:dim(P1)[1]){
	P1[i,] <- pnorm(as.numeric(log(rBRG1.B[i,])),mean = 0, sd = WT_variability, lower.tail=F)
   	P2[i,] <- pnorm(as.numeric(log(rBRG1.B[i,])),mean = 0, sd = WT_variability, lower.tail=T)
        for(j in 1:dim(P1)[2])
                PB[i,j] <- min(P1[i,j], P2[i,j])
   }

   P <- PA
   for(i in 1:dim(PA)[1])
        for(j in 1:dim(PA)[2])
		if (is.na(rBRG1.A[i,j]) || is.na(rBRG1.B[i,j]))
			P[i,j] = min(PA[i,j],PB[i,j],na.rm=T)	
		else if	(log(rBRG1.A[i,j])*log(rBRG1.B[i,j])<0)
			P[i,j] = 1.0
   		else
			P[i,j] = pchisq( -2*sum(log(c(PA[i,j],PB[i,j]))), 4, lower.tail=FALSE)

   P <- P[,columns]
   P2 <- P
   oo <- order(P2)
   N <- nrow(P2)*ncol(P2)
   for(i in 1:N) P2[oo[i]] <- min(P[oo[i]] * N / i, 1.0)
   P2
}


getF1 <- function(replicate, ratioData = BRG1.rr, absoluteData = BRG1.ab,suffix="")
{
  print("reading data")
  name <- paste("Abundance.",replicate,sep="")
  f1_cols <- grep(name, colnames(ratioData), perl=TRUE, value=FALSE)
  f1_ab_col <- grep(paste(name,".*.WT",suffix,"$",sep=""), colnames(absoluteData), perl=TRUE, value=FALSE)
  #print(f1_ab_col)
  accesionID <- which(colnames(ratioData)=="Accession")
  markedAsID <- which(colnames(absoluteData)=="Marked.as")
  F1 <- cbind(absoluteData[,c(4,markedAsID,5,16,f1_ab_col)],ratioData[,c(accesionID,f1_cols)])
  #print(F1[1,])

  exclude_list <- c()#"Q8WUB8","Q86U86","Q9NPI1","Q9H8M2","Q9H165","Q68CP9","Q8NFD5","P60709")

  baf <- which(F1[,2]=="BAF_complex" & !(F1[,6] %in% exclude_list))
  N <- length(baf)

  F1_data <- matrix(NA,nrow=N,ncol=11)
  #print(dim(F1))
  for(i in 1:N) { 
    F1_data[i,1] = 0 #F1[baf[i],4]; 
    #ss = sum(F1[baf[i],7:16]); 
    for(j in 1:10) { 
      F1_data[i,j+1] = F1[baf[i],6+j] #* F1[baf[i],5]/ss; 
    } 
  }

  #myRNA <- which(F1[,1] %in% rna_proteins)
  #print(myRNA)

  #for(j in 1:10) {
  #	F1_data[N+1,j+1] = median(F1[myRNA,6+j],na.rm=T)
  #}
  #print(F1_data[N,2:11])
  #print(F1_data[N+1,2:11])

  relF1_data <- F1_data[,2:11]
  #for(i in 1:(N)) 
  #  relF1_data[i,] = relF1_data[i,] / max(relF1_data[i,])

  WR <- grep(paste("WT.R",suffix,sep=""), colnames(ratioData)[f1_cols], perl=TRUE, value=FALSE)
  #WR <- grep(paste(".WT",suffix,"$",sep=""), colnames(ratioData)[f1_cols], perl=TRUE, value=FALSE)
  colWR <- relF1_data[,WR]
  for(j in 1:10)
    relF1_data[,j] = relF1_data[,j]/colWR

  relF1_data
}

production <- function()
{

   relF1_data <- getF1("F1",BRG1.rr,BRG1.rr)
   relF2_data <- getF1("F2",BRG1.rr,BRG1.rr)
   relF3_data <- getF1("F3",BRG1.rr,BRG1.rr)
   relF4_data <- getF1("F4",BRG1.rr,BRG1.rr)   


   pdf("BRG1_pulldowns.pdf",width=14,height=7)
   barPlot(relF1_data, relF1_data, ymax=1.6,logV="",start=2,legends=labels_f1, legendPos=1.0)
   barPlot(relF2_data, relF2_data, ymax=1.6,logV="",start=2,legends=labels_f1, legendPos=1.2)
   barPlot(relF3_data, relF3_data, ymax=1.6,logV="",start=3,legends=labels_f3, legendPos=1.3)
   barPlot(relF4_data, relF4_data, ymax=1.6,logV="",start=3,legends=labels_f3, legendPos=1.3)
   dev.off()

   relF1_data_A1A <- getF1("F1",ARID1A.rr,ARID1A.rr,suffix=".2")
   relF2_data_A1A <- getF1("F2",ARID1A.rr,ARID1A.rr,suffix=".2")
   relF3_data_A1A <- getF1("F3",ARID1A.rr,ARID1A.rr,suffix=".1")
   relF4_data_A1A <- getF1("F4",ARID1A.rr,ARID1A.rr,suffix=".1")

   pdf("ARID1A_pulldowns.pdf",width=14,height=7)
   barPlot(relF1_data_A1A, relF1_data_A1A, ymax=1.6,logV="",start=2,legends=labels_f1_A1A, legendPos=1.3, labels = labels_A1A)
   barPlot(relF2_data_A1A, relF2_data_A1A, ymax=1.6,logV="",start=2,legends=labels_f1_A1A, legendPos=1.6, labels = labels_A1A)
   barPlot(relF3_data_A1A, relF3_data_A1A, ymax=1.6,logV="",start=3,legends=labels_f3_A1A, legendPos=0.9, labels = labels_A1A)
   barPlot(relF4_data_A1A, relF4_data_A1A, ymax=1.6,logV="",start=3,legends=labels_f3_A1A, legendPos=1.3, labels = labels_A1A)
   dev.off()
}

plotHeatmaps <- function(filename="heatmaps.pdf"){

    relF1_data <- scale(getF1("F1",BRG1.rr,BRG1.rr))
    relF2_data <- scale(getF1("F2",BRG1.rr,BRG1.rr))
    relF3_data <- scale(getF1("F3",BRG1.rr,BRG1.rr))
    relF4_data <- scale(getF1("F4",BRG1.rr,BRG1.rr))

    pdf(filename,height=7,width=14)
    myLabels=c(paste(labels_f1,"A",sep="-"),paste(labels_f1,"B",sep="-"),paste(labels_f3,"A",sep="-"),paste(labels_f3,"B",sep="-"))
    pheatmap(2*cbind(relF1_data,relF2_data,relF3_data,relF4_data),labels_col=myLabels,labels_row=labels_BRG1_full,scale='none',col=color,main="BRG1 pulldowns",legend_breaks=c(0.023,0.5,1,1.5,2),legend_labels=c("0","0.5xmedian","median","1.5xmedian","2xmedian"))

   relF1_data_A1A <- scale(getF1("F1",ARID1A.rr,ARID1A.rr,suffix=".2"))
   relF2_data_A1A <- scale(getF1("F2",ARID1A.rr,ARID1A.rr,suffix=".2"))
   relF3_data_A1A <- scale(getF1("F3",ARID1A.rr,ARID1A.rr,suffix=".1"))
   relF4_data_A1A <- scale(getF1("F4",ARID1A.rr,ARID1A.rr,suffix=".1"))

    myLabels=c(paste(labels_f1_A1A,"A",sep="-"),paste(labels_f1_A1A,"B",sep="-"),paste(labels_f3_A1A,"A",sep="-"),paste(labels_f3_A1A,"B",sep="-"))
    pheatmap(2*cbind(relF1_data_A1A,relF2_data_A1A,relF3_data_A1A,relF4_data_A1A)[c(1:12,14:27,29),],labels_col=myLabels,labels_row=labels_A1A_full[c(1:12,14:27,29)],scale='none',col=color, main="ARID1A pulldowns",legend_breaks=c(0.023,0.5,1,1.5,2),legend_labels=c("0","0.5xmedian","median","1.5xmedian","2xmedian"))
   #pheatmap(2*cbind(relF1_data_A1A,relF2_data_A1A,relF3_data_A1A,relF4_data_A1A)[c(1:8,10:19,21),],labels_col=myLabels,labels_row=labels_A1A_full[c(1:8,10:19,21)],scale='none',col=color, main="ARID1A pulldowns",legend_breaks=c(0.023,0.5,1,1.5,2),legend_labels=c("0","0.5xmedian","median","1.5xmedian","2xmedian"))
    dev.off()
}

plotRatioHeatmaps <- function(filename="heatmaps2.pdf"){

	pdf(filename,height=7,width=9)
	rBRG1.A <- getRatios2(ratioData=BRG1.rr2,replicate="A",exclude_list=exclude_list_BRG1,capAt2=F)
	rBRG1.B <- getRatios2(ratioData=BRG1.rr2,replicate="B",exclude_list=exclude_list_BRG1,capAt2=F)
        cols <- which(colnames(rBRG1.A)!="WT.R")

	colnames(rBRG1.A) <- paste(colnames(rBRG1.A),"-A",sep="")
        colnames(rBRG1.B) <- paste(colnames(rBRG1.B),"-B",sep="")
	data <- cbind(rBRG1.A[,cols],rBRG1.B[,cols])
	oC <- order(colnames(data))
        oR <- order(rownames(data))
	#pheatmap(data[oR,oC],col=c(color2),breaks=seq(0,2,length.out=101),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main="BRG1 pulldowns",cluster_cols =F, cluster_rows=F)
	for(i in 1:dim(data)[1]){
	   data[i,which(data[i,]>4)]=4
	   data[i,which(data[i,]<1/4)]=1/4
	}
        pheatmap(log(data[oR,oC],base=2),col=c(color2),breaks=seq(-2,2,length.out=101),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main="BRG1 pulldowns",cluster_cols =F, cluster_rows=F)

	rARID1A.B <- getRatios2(ratioData=ARID1A.rr2,replicate="B",exclude_list=exclude_list_ARID1A,capAt2=F)
	rARID1A.A <- getRatios2(ratioData=ARID1A.rr2,replicate="A",exclude_list=exclude_list_ARID1A,capAt2=F)
	cols <- which(colnames(rARID1A.A)!="WT.R")

	data <- cbind(rARID1A.A[,cols],rARID1A.B[,cols])
	validrows <- c()
	Nrow <- dim(data)[1]
	Ncol <- dim(data)[2]
	for(i in 1:Nrow)
		if(sum(is.na(data[i,]))!=Ncol)
			validrows <- c(validrows,i)
	#print(validrows)
	colnames(rARID1A.A) <- paste(colnames(rARID1A.A),"-A",sep="")
	colnames(rARID1A.B) <- paste(colnames(rARID1A.B),"-B",sep="")
	data <- cbind(rARID1A.A[,cols],rARID1A.B[,cols])[validrows,]
	oC <- order(colnames(data))
        oR <- order(rownames(data))
	#pheatmap(data[oR,oC],col=c(color2),breaks=seq(0,2,length.out=101),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main="ARID1A pulldowns",cluster_cols =F, cluster_rows=F)
	for(i in 1:dim(data)[1]){
           data[i,which(data[i,]>4)]=4
           data[i,which(data[i,]<1/4)]=1/4
        }
	pheatmap(log(data[oR,oC],base=2),col=c(color2),breaks=seq(-2,2,length.out=101),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main="ARID1A pulldowns",cluster_cols =F, cluster_rows=F)
	dev.off()
}

plotGroupedRatioHeatmaps <- function(filename="heatmapGrouped.pdf",data,exclude_list,name, columns = c(4:11,15:30), pval=NA){
	pdf(filename,height=7,width=9)
        rBRG1.A <- getGroupedRatios2(ratioData=data,exclude_list=exclude_list)
	w <- which(data[,"Marked.as"]=="BAF_complex" & !(data[,"Accession"] %in% exclude_list))

	if (is.na(pval)){
	   pval <- getRatiosPvalue2(ratioData=data[w,],exclude_list=exclude_list, columns=columns)
	}

        #pval <- pval[w,]
	rval <- pval
	N1 <- dim(pval)[1]
	N2 <- dim(pval)[2]
	for(i in 1:N1)
	   for(j in 1:N2){
		if(!is.na(pval[i,j])) rval[i,j] <- round(pval[i,j],digits=2)
		if(is.na(pval[i,j]) || pval[i,j]>0.01) rval[i,j] = ""
		if(!is.na(pval[i,j]) && pval[i,j]<=0.01) rval[i,j] = "*"
		if(!is.na(pval[i,j]) && pval[i,j]<=0.001) rval[i,j] = "**"
		if(!is.na(pval[i,j]) && pval[i,j]<=0.0001) rval[i,j] = "***"
		#if((rBRG1.A[,columns])[i,j]==1.0) rval[i,j]=""
	   }

	oC <- order(colnames(rBRG1.A[,columns]))
	oR <- order(rownames(rBRG1.A[,columns]))
	pheatmap(rBRG1.A[,columns][oR,oC],col=c(color2),breaks=seq(0,2,length.out=101),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main=name,display_numbers = rval[oR,oC],fontsize_number=15,cluster_cols =T, cluster_rows=T)
	write.table(rBRG1.A[,columns][oR,oC],file="data.csv",sep="\t",quote=F)
	write.table(pval[oR,oC],file="pval.csv",sep="\t",quote=F)
	dev.off()
}

plotPeptideHeatmaps <- function(filename="heatmaps_peptides_BRG1.pdf",threshold=0.5){
     pdf(filename,height=18,width=9)
     
     for (protein in labels_BRG1_full){
     pBRG.A <- getPeptideRatios2(protein, replicate="A", ratioData=BRG.pep3, proteinLabels=BRG1.rr2[,"Description"], protein.accesions=BRG1.rr2[,"Accession"])
     pBRG.B <- getPeptideRatios2(protein, replicate="B", ratioData=BRG.pep3, proteinLabels=BRG1.rr2[,"Description"], protein.accesions=BRG1.rr2[,"Accession"])
     if (!is.null(pBRG.A) && !is.null(pBRG.B)){
     N <- dim(pBRG.A)[1]
     print(N) 
     cc <-rep(NA,N) 
     for(i in 1:N) 
        cc[i] <- tryCatch(cor(as.numeric(pBRG.A[i,c(1:4,7:30)]),as.numeric(pBRG.B[i,c(1:4,7:30)]),use="complete.obs"), error= function(e) {NA})

     w <- which(cc>threshold)
     if (length(w>1)){
     data <- cbind(pBRG.A[w,],pBRG.B[w,])
     N <- dim(data)[1]
     M <- dim(data)[2]
     dd <- c()
     for(i in 1:N){ 
	#data[i,which(is.na(data[i,]))]=2.1
	if(sum((data[i,])==0.01,na.rm=T)!=M && sum((data[i,])==1,na.rm=T)!=M) dd <- c(dd,i)
     }
     data <- data[dd,]

     cellheight <- min(1000/(dim(data)[1]),20)
     print(c(cellheight, dim(data), N))
     if(dim(data)[1]>=2)
     pheatmap(data[,c(7,37,8,38,9,39,10,40,11,41,12,42,13,43,14,44,15,45,16,46,17,47,18,48,19,49,20,50,21,51,22,52,23,53,26,56,27,57,28,58,29,59,30,60)],col=c("grey",color2),breaks=seq(0,2,length.out=102),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main=paste("Peptides of ",protein," in BRG1 pulldowns.",sep=""), fontsize_row=7,cellheight=cellheight,cluster_cols = F)

     }
     }
     }

     dev.off()
}

plotPeptideHeatmaps_BRG1 <- function(proteinInput=BRG1.rr2, peptideInput=BRG.pep3, filename="heatmaps_peptides_BRG1_signif.pdf", threshold=0.5, relevantColums=c(7:23,26:30), exclude_list = exclude_list_BRG1, label="BRG1",proteinList = labels_BRG1_full ){
    plotPeptideHeatmaps_ARID1A(proteinInput, peptideInput, filename, threshold=0.5, relevantColums, exclude_list, label, proteinList)
}

plotPeptideHeatmaps_ARID1A <- function(proteinInput=ARID1A.rr2, peptideInput=ARID1A.pep3, filename="heatmaps_peptides_ARID1A_signif.pdf", threshold=0.5, relevantColums=c(4:7,10:30), exclude_list = exclude_list_ARID1A, label="ARID1A", proteinList = labels_A1A_full){
     pdf(filename,height=18,width=9)

     rrA <- getRatios2(ratioData=proteinInput,replicate="A",exclude_list=exclude_list)
     rrB <- getRatios2(ratioData=proteinInput,replicate="B",exclude_list=exclude_list)

     for (protein in proteinList){
     print(protein)
     p.A <- getPeptideRatios2(protein,replicate="A",ratioData=peptideInput, proteinLabels=proteinInput[,"Description"], protein.accesions=proteinInput[,"Accession"])
     p.B <- getPeptideRatios2(protein,replicate="B",ratioData=peptideInput, proteinLabels=proteinInput[,"Description"], protein.accesions=proteinInput[,"Accession"])
     #print(dim(p.A))
     if (!is.null(p.A) && !is.null(p.B)){
     N <- dim(p.A)[1]
     
     cc <-rep(NA,N)
     for(i in 1:N) 
     cc[i] <- tryCatch(cor(as.numeric(p.A[i,relevantColums]),as.numeric(p.B[i,relevantColums]),use="complete.obs"), error= function(e) {NA})

     w <- which(cc>=threshold)
     if (length(w>1)){
     
     data <- cbind(p.A[w,],p.B[w,])
     N <- dim(data)[1]
     M <- dim(data)[2]
     ddt <- c()
     for(i in 1:N){
        #data[i,which(is.na(data[i,]))]=2.1
        if(sum((data[i,])==0.01,na.rm=T)!=M && sum((data[i,])==1,na.rm=T)!=M) ddt <- c(ddt,i)
     }
     data <- data[ddt,]
     #print(sort(rownames(data)))

     #if(dim(data)[1]>=2){
	dd <- data[,relevantColums]
        dd2 <- data[,30+relevantColums]
	data2 <- sqrt(dd * dd2)
	#print ("HERE3")
	#print (dim(data2))	
	N1 <- dim(data2)[1]
        N2 <- dim(data2)[2]
        pval <- matrix(1,N1,N2)
	rProt <- pval;

	for(i in 1:N1)
	   for(j in 1:N2){
		col <- substring(colnames(data2)[j],1,nchar(colnames(data2)[j])-2)
		rProtein <- sqrt( rrA[protein,col] * rrB[protein,col]  )
		if (!grepl("WT",col) && !is.na(dd[i,j]) && !is.na(dd2[i,j]) && !is.na(rProtein) && dd[i,j]!=0 && dd2[i,j]!=0){
		  #print(c(col,dd[i,j],dd2[i,j],rProtein))
		  pval[i,j] <- getPvalue(dd[i,j],dd2[i,j],rProtein,sd=0.3936362)
		  rProt[i,j] <- rProtein
		}
	   }
	#oo <- order(pval)
        #for(i in 1:(N1*N2)) 
	#  pval[oo[i]] <- min(pval[oo[i]] * N1*N2 / i, 1.0)
	
	rval <- pval
	Rows2Show <- rep(0,N1)
	#print("HERE")
        for(i in 1:N1)
           for(j in 1:N2){
		if (!is.na(dd[i,j]) && !is.na(dd2[i,j]) &&  dd[i,j]!=0.015 && dd2[i,j]!=0.015){
		  if(pval[i,j]<=0.01) { 
			print (paste(rownames(data2)[i],colnames(data2)[j], data2[i,j], pval[i,j], rProt[i,j],sep=" "))
			Rows2Show[i]=i
			peptide <- rownames(data2)[i]
		        s1 <- regexpr(":",peptide)[1]	
			s2 <- regexpr("-",peptide)[1] 
			start <- strtoi(substring(peptide,1,s2-1))
			end <- strtoi(substring(peptide,s2+1,s1-1))
			for(k in 1:N1){
			    pp <- rownames(data2)[k]
			    s1 <- regexpr(":",pp)[1]
                            s2 <- regexpr("-",pp)[1]
                            start2 <- strtoi(substring(pp,1,s2-1))
                            end2 <- strtoi(substring(pp,s2+1,s1-1))
			    if (abs(start-start2)<=1 && abs(end-end2)<=1) Rows2Show[k]=i

			}
		  }
		#if(pval[i,j]<0.01) print (paste(rownames(data2)[i],colnames(data2)[j], data2[i,j], pval[i,j], rProt[i,j],sep=" "))
                #if(!is.na(pval[i,j])) rval[i,j] <- round(pval[i,j],digits=2)
                  if(is.na(pval[i,j]) || pval[i,j]>0.01) rval[i,j] = ""
                  if(!is.na(pval[i,j]) && pval[i,j]<=0.01) rval[i,j] = "*"
                  if(!is.na(pval[i,j]) && pval[i,j]<=0.001) rval[i,j] = "**"
                  if(!is.na(pval[i,j]) && pval[i,j]<=0.0001) rval[i,j] = "***"
		}
		else rval[i,j]=""
                #if((rBRG1.A[,columns])[i,j]==1.0) rval[i,j]=""
           }
	#print("HERE2")
	if(length(which(Rows2Show!=0))>=1){

	  print("... heatmap ...")

	  data2 <- data2[which(Rows2Show!=0),]
	  rval <-   rval[which(Rows2Show!=0),]
	  Rows2Show <- Rows2Show[which(Rows2Show!=0)]

	  # heatmap cannot show just a single row... adding a dummy row
	  if(length(which(Rows2Show!=0))==1){
	    data2 <- rbind(data2,data2)
	    rval <- rbind(rval,rval)
	    Rows2Show <- 1:2
	  }

	  oo <- order(Rows2Show)

	  cellheight <- max(min(1000/length(oo),20),5)
	  pheatmap(data2[oo,],col=c("grey",color2),breaks=seq(0,2,length.out=102),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main=paste("Peptides of ",protein," in ",label," pulldowns.",sep=""), fontsize_row=7,cellheight=cellheight, cluster_cols = F, cluster_rows=F, display_numbers= rval[oo,], number_color = "green")
	}
	else{
	   print(paste("Nothing significant for",protein,length(which(Rows2Show!=0))))
	   cellheight <- max(min(1000/N1,20),5)
	   #print(dim(data2))
	   #print(cellheight)
	   #pheatmap(data2,col=c("grey",color2),breaks=seq(0,2,length.out=102),norm="none",clustering_distance_cols="euclidean",clustering_distance_rows="euclidean",main=paste("Nothing significant for ",protein," in ",label," pulldowns.",sep=""), fontsize_row=7,cellheight=cellheight, cluster_cols = F, cluster_rows=F, display_numbers= rval, number_color = "red")
	}

     #}
     #else print(paste("not enought observations for",protein))
     }
     }
     }
     dev.off()
}

pulledAbundance <- function(){
   #pdf("ARID1A.pdf",height=7,width=10)
   #	arid1APlot(legendPos=-38000)
   pdf("ARID1A-3.pdf",height=7,width=10)
	arid1APlot(legendPos=-38000,data=ARID1A.rr2,labelprefix=27)
   dev.off()

   #pdf("BRG1.pdf",height=7,width=10)
   #     arid1APlot(ID="P51532",data=BRG1.rr,labelprefix=29,legendPos=-21000,name="SMARCA4")
   pdf("BRG1-3.pdf",height=7,width=10)
   	arid1APlot(ID="P51532",data=BRG1.rr2,labelprefix=27,legendPos=-21000,name="SMARCA4")
   dev.off()

   pdf("BAF_complex.pdf",height=7,width=10)
   	barPlot3()
   dev.off()
}

scale <- function(data){
   scaled_data <- data
   columns <- dim(data)[2]
   rows <- dim(data)[1]
   for(j in 1:columns){
     maxV <- max(1.0,max(2*median(data[,j],na.rm=T),min(3.0,max(data[,j],na.rm=T))))
     #maxV <- max(0.5,2*median(data[,j],na.rm=T))
     #print(maxV)
     scaled_data[,j] <- data[,j]/maxV
     for(i in 1:rows) scaled_data[i,j] <- min(scaled_data[i,j],1.0)
   }
   scaled_data
}



