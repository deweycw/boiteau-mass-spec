#Rene Boiteau, 12/21/20 - File for Batch processing of mzXML files. 

library(xcms)

########################   Input settings (Change settings here)

timerange=c(0,30*60) #time range in seconds
background<-3000 #minimum signal intensity
peakwidth<-30 #expected peak width in seconds (larger number is less stringent for discarding false positives)
slope<-c(0,2) #minimum normalized slope cutoff.(Light/Heavy)
correlation<-0.7 #minimum correlation cutoff.
eicwidth<-0.003 #m/z width of extracted ion chromatogram (+/- Da)

filetype<-".mzML" #define file tag for mass spec files. 

######################## Files are loaded here

#Set list of Mass Spectrometry files for analysis
mzdatafiles <- list.files(recursive=FALSE, full.names=TRUE, pattern=paste("\\",filetype,sep=""))
#Set list of Isotope pattern files for analysis
isotopes <- list.files(recursive=FALSE, full.names=TRUE, pattern="\\isotope.csv")

######################## Functions are defined here


###########Function to find isotope matches

##Note, also uses 'background'
isotopehunter8<-function(mzxcms,isotopefile,scantime){
  
  starttime<-Sys.time()
  
  if(scantime[1]=="all"){
    startscan<-1
    endscan<-length(mzxcms@scantime)
  }else{
    startscan<-which(mzxcms@scantime>scantime[1])[1]
    endscan<-tail(which(mzxcms@scantime<scantime[2]),1)
  }
  
  readpattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  pattern<-readpattern[(readpattern[,6]=='Y'),]
  
  uppermass<-pattern[-1,2]+pattern[-1,4]-pattern[1,2]
  lowermass<-pattern[-1,2]-pattern[-1,4]-pattern[1,2]
  upperratio<-pattern[-1,3]*pattern[-1,5]/pattern[1,3]
  lowerratio<-pattern[-1,3]/pattern[-1,5]/pattern[1,3]
  nisotope<-length(uppermass)
  
  columnnames<-c(sapply(pattern[,1],FUN=function(x) c(paste(x,'mass'),paste(x,'intensity'))))
  
  results<-matrix(0,ncol=(4+nisotope*2),nrow=1E6)
  r<-1
  colnames(results) = c('scan','time', columnnames)
  
  
  mzint<-20
  mzbuffer<-uppermass[length(uppermass)]
  
  for(i in startscan:endscan) {
    tscan<-getScan(mzxcms,i)
    tscan<-tscan[which(tscan[,2]>background),]
    
    #print(paste(i,'of',endscan))
    
    first<-tscan[1,1]
    last<-first+mzint
    
    while(last<tscan[nrow(tscan),1]){
      
      scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
      
      if(length(scan)>2){
        final<-length(which(scan[,1]<(last-mzbuffer)))
        for(j in 1:final){
          x<-scan[j,1]
          y<-scan[j,2]
          isotopes<-matrix(NA,ncol=nisotope)
          k<-1
          while(k<nisotope+1){
            isotopes[k]<-(which(
              scan[,1]>(x+lowermass[k]) & 
                scan[,1]<(x+uppermass[k]) & 
                scan[,2]>(y*lowerratio[k]) & 
                scan[,2]<(y*upperratio[k])
            )[1])
            if(is.na(isotopes[k])){k<-(nisotope+2)
            } else {k<-k+1
            }
          }
          if(k==(nisotope+1)){
            results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
            r<-r+1
          }
          
        }
      }
      first<-last-mzbuffer
      last<-last+mzint
    }
    scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
    if(length(scan)>2){
      for(j in 1:nrow(scan)){
        x<-scan[j,1]
        y<-scan[j,2]
        isotopes<-matrix(NA,ncol=nisotope)
        k<-1
        while(k<nisotope+1){
          isotopes[k]<-(which(
            scan[,1]>(x+lowermass[k]) & 
              scan[,1]<(x+uppermass[k]) & 
              scan[,2]>(y*lowerratio[k]) & 
              scan[,2]<(y*upperratio[k])
          )[1])
          if(is.na(isotopes[k])){k<-(nisotope+2)
          } else {k<-k+1
          }
        }
        if(k==(nisotope+1)){
          results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
          r<-r+1
        }
      }
    }
  }
  
  print(Sys.time()-starttime)
  
  return(results[1:(r-1),])
}

########Simplify list:

featurebinning<-function(results,mzxcms,background,peakwidth){
  
  masslist2<-c()

      #create index row for results
      results2<-results[order(results[,4],decreasing=TRUE),]
      
      #Remove elements that don't appear at least 2x
      results2<-results2[!isUnique(round(results2[,5],digits=3)),]
      
      #Remove duplicate masses, selecting max intensity scan
      results2<-results2[!duplicated(round(results2[,5],digits=3)),]
      
      #Calculate correlation coefficient for each element
      slope<-0
      rsquare<-0  
      
      if(class(results2)=='matrix'){
        if(nrow(results2)>1){
          for(m in 1:nrow(results2)){
            EIC1<-rawEIC(mzxcms,mzrange=c(-eicwidth,eicwidth)+results2[m,3],rtrange=c(-peakwidth,peakwidth)+results2[m,2])
            EIC2<-rawEIC(mzxcms,mzrange=c(-eicwidth,eicwidth)+results2[m,5],rtrange=c(-peakwidth,peakwidth)+results2[m,2])
            peakcorr<-lm(EIC1$intensity ~ EIC2$intensity)
            slope[m]=coef(peakcorr)[2]
            rsquare[m]=summary(peakcorr)$r.square
          }
        }
      }
      results3<-cbind(results2,slope,rsquare)
      
      return(data.frame(results3))
}

################# Plotting functions

#Isotope plotter for final report. Makes one plot with EIC's of every isotope found, scaled to the 'same' intensity.
isotopeplotter<-function(mzxcms,pattern,lowmass,timerange,maxscan){
  
  plotrange<-pattern[,2]-pattern[1,2]+lowmass
  isotoperatio<-pattern[,3]/pattern[1,3]
  namerange<-as.vector(pattern[,1])
  print(namerange)
  print(isotoperatio)
  colors<-c('blue4','darkorange2','green',  'burlywood4','black','red')#
  usedcolors<-colors[1:ncol(pattern)]
  print(usedcolors[2])
  maxy<-0
  times<-mzxcms@scantime
  
  for(i in 1:nrow(pattern)){
    mzrange<-c(-eicwidth,eicwidth)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    EIC1<-unlist(EIC[2])/isotoperatio[i]    #EIC1<-unlist(EIC[2])/isotoperatio[i] <-- original code; change made 2/6/21 CWD
    newmax<-max(EIC1[which(times>timerange[1]&times<timerange[2])])
    maxy<-max(c(maxy,newmax))
  }
  
  print(pattern)
  mzrange<-c(-eicwidth,eicwidth)+lowmass
  EIC<-rawEIC(mzxcms,mzrange)
  #print(EIC[[2]])
  plot(times,EIC[[2]],
       type='l',
       lwd=1,
       xlim=timerange,
       ylim=c(0,maxy*1.1),
       ylab='Scaled Intensity',
       xlab='Retention Time (min)',
       col=usedcolors[1],
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  
  for(i in 2:ncol(pattern)){
    mzrange<-c(-eicwidth,eicwidth)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    #print(i)
    #print(isotoperatio[i])
    #print(usedcolors[i])
    lines(times,EIC[[2]],col=usedcolors[i],lwd=1)  #lines(times,EIC[[2]]/isotoperatio[i],col=usedcolors[i],lwd=1)  <-- original code; change made 2/6/21  CWD
  }
  
  abline(v=maxscan, lty=3,col='gray48')
  title(paste('EIC = ',toString(round(plotrange,digits=4))),line=0,adj=0,cex.main=1)
  legend('topright',bty="n",namerange,lwd=2,col=usedcolors,cex=0.7,horiz=TRUE)
  print(namerange[2])
  print(usedcolors[2])
}

#EIC plotter for final report.
EICplot<-function(mzxcms,mass,title,timerange,maxscan){
  mzrange<-c(-eicwidth,eicwidth)+mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime
  plot(times,EIC[[2]],
       type='l',
       lwd=1,
       xlim=timerange,
       ylab='Intensity',
       xlab='Retention Time (min)',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  title(paste(title,'EIC = ',round(mass,digits=4)),line=1,adj=0,cex.main=1)
  
  abline(v=maxscan, lty=3,col='gray48')
}

#Mass spectra plotter for final report.
MSplotter<-function(mzxcms,results,i,pattern){
  
  #determine correct isotope pattern
  theorymass<-pattern[,2]-pattern[1,2]+results[i,3]
  theoryratio<-pattern[,3]/pattern[1,3]*results[i,4]
  
  isotopes<-matrix(results[i,-1:-2],byrow=TRUE,ncol=2)
  massrange<-c(mean(isotopes[,1])-10,mean(isotopes[,1])+10)
  scanrange<-getScan(mzxcms,results[i,1],mzrange=massrange)
  
  plot(scanrange,type='h',ylim=c(0,max(scanrange[,2])*1.1),bty='n',ylab='Intensity',xlab='m/z')
  lines(theorymass,theoryratio,type='h',col='cornsilk3',lwd=3)
  lines(isotopes,type='h',col='red')
  title(paste('MS: ',toString(round(results[i,2]/60,digits=2)),' min'),line=1,adj=0,cex.main=1)
}



#Final plot function for generating report for isotope pattern. Includes Apo form MS spectra as last panel
Reportplotter<-function(mzxcms,mass,results,isotopepattern,timerange,i){
  
  mzxcms@scantime<-mzxcms@scantime/60
  timerange_m<-timerange/60
  readpattern<-read.csv(isotopepattern, header = TRUE, sep = ",")
  pattern<-readpattern[which(readpattern[,6]=='Y'|readpattern[,6]=='N'),]
  apo<-readpattern[which(readpattern[,6]=='O'),]
  
  nplot<-2+nrow(apo)
  par(mfrow=c(nplot,1))
  par(mar=c(4,4,3,3))
  
  isotopeplotter(mzxcms,pattern,mass,timerange_m,(results[i,2]/60))
  title(savefile)
  #title(savefile,outer=TRUE)
  title(paste(' R-squared: ',toString(round(results$rsquare[i],digits=3))),line=-1,adj=0.0,cex.main=1)
  title(paste(' Slope: ',toString(round(results$slope[i],digits=3))),line=-2,adj=0.0,cex.main=1)
  title(paste(' d m/z: ',toString(round((results[i,5]-results[i,3]),digits=4))),line=-3,adj=0.0,cex.main=1)
  
  if(nrow(apo)>0){
    for(m in 1:nrow(apo)){
      EICplot(mzxcms,mass+apo[m,2],apo[m,1],timerange_m,(results[i,2]/60))
    }
  }
  #print(results)
  MSplotter(mzxcms,as.matrix(results[,1:(ncol(results)-2)]),i,pattern) ##
}

# Initialize Pipeline for algorithm search, correlation analysis, and report generation.

Isotope_pipeline <- function(isotopepattern,mzxcms,timerange,background,peakwidth,savefile){
  
  pattern<-read.csv(isotopepattern, header = TRUE, sep = ",")
  ratio<-pattern[2,3]/pattern[1,3]
  
  #Search for all matches to the defined isotope pattern
  results<-isotopehunter8(mzxcms,isotopepattern,timerange)
  print(paste('Results:',nrow(results)))
  
  
  if(class(results)=='matrix'){
  
    #Bin results by mass and remove patterns that only appear 1x:
    binnedresults<-featurebinning(results,mzxcms,background,peakwidth)
    binnedresults$slope<- binnedresults$slope*ratio
    
    #Filter clean results based on R-squared value (should be close to 1 for best hits)
    cleanresults<-binnedresults[which(binnedresults$rsquare>correlation),]
    
    #Filter clean results based on normalized slope (should be 1 if peak ratio really matches isotope pattern)
    cleanresults<-cleanresults[which(cleanresults$slope<slope[2] &
                                         cleanresults$slope>slope[1]),]
    
    print(paste("Cleaned results:",nrow(cleanresults)))
    
    #Generate report: Fix this cleanresults part at some point...Right now, wont generate report if just one hit...
  
    if(nrow(cleanresults)>0){
      #Will generate pdf w/ title savefile.
      pdf(paste(savefile,'.pdf'))
      
      #First, creates QC page:
      par(mfrow=c(3,1))
      #Histogram of Rsquared distribution
      hist(binnedresults$rsquare,breaks=seq(0,1,0.05),xlab='R-squared',main='R-squared Distribution')
      abline(v=correlation,col='red')
      #Histogram of Slope distribution
      hist(cleanresults$slope,breaks=seq(slope[1],slope[2],0.05),xlab='Slope (light/heavy)',main='Normalized Slope Distribution')
      abline(v=slope,col='red')
      
      ## Plot m/z and intensity ratio of all results vs. clean results..
      plot(results[,5]-results[,3],results[,4]/results[,6],col='lightgray',xlab='Mass difference',ylab='Intensity ratio (light/heavy)')
      points(binnedresults[,5]-binnedresults[,3],binnedresults[,4]/binnedresults[,6],col='black')
      points(cleanresults[,5]-cleanresults[,3],cleanresults[,4]/cleanresults[,6],col='red')
      legend('top',legend=c('All','Binned','Cleaned'),col=c('lightgray','black','red'),pch=1, xpd=TRUE, horiz=TRUE, inset=c(0,-.2),bty='n')
  
      
      masslist<-cleanresults[,3]
      
      
      for(k in 1:length(masslist)){
        Reportplotter(mzxcms,masslist[k],cleanresults,isotopepattern,timerange,k)
      }
      
      dev.off()
      
      #write.table(cleanresults,paste(savefile,'.txt'),sep="\t")
    }
    return(cleanresults)
  }
}

######################### To run batch mode:

#Loop for each MS data file:
for(i in 1:length(mzdatafiles)){
  
  #Load MS data file as an xcmsRaw object:
  mzxcms <- xcmsRaw(mzdatafiles[i],profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
  
  #Loop for each isotope pattern file:
  for(j in 1:length(isotopes)){
    
    isotopepattern<-isotopes[j]
    
    #Creates name of save file    
    savefile<-paste(gsub(paste("^./(.*)",filetype,sep=""),"\\1",mzdatafiles[i]),"_",gsub("^./(.*).csv","\\1",isotopes[j]),sep="")
    print(savefile)
    
    #Run isotope pipeline
    Isotope_pipeline(isotopepattern,mzxcms,timerange,background,peakwidth,savefile)
  }
}


