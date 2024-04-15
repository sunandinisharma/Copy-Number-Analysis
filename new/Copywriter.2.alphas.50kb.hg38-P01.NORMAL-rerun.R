.libPaths("/home/javeediqbal/abouska/R/x86_64-pc-linux-gnu-library/3.5")
library(CopywriteR)


args <- commandArgs(TRUE)
print(args)
pfx<-args[1]
norm<-args[2]
print(pfx)
print(norm)

samples<-c(paste(pfx,"/Alignment/",pfx,".recalibrated.final.bam",sep=""),paste("/lustre/work/javeediqbal/shared/weiwei/P01/",norm,"/Alignment/",norm,".recalibrated.final.bam",sep=""))
controls<-c(paste("/lustre/work/javeediqbal/shared/weiwei/P01/",norm,"/Alignment/",norm,".recalibrated.final.bam",sep=""),paste("/lustre/work/javeediqbal/shared/weiwei/P01/",norm,"/Alignment/",norm,".recalibrated.final.bam",sep=""))
sample.control <- data.frame(samples, controls)
sample.control<-unique(sample.control)

print(sample.control)

bp.param <- SnowParam(workers =1, type = "SOCK")
outPath<-paste("/lustre/work/javeediqbal/shared/MCL_Sunandini/MCL_WES/",pfx,"/", sep="")

refPath<-"/lustre/work/javeediqbal/abouska/CopywriteR/R3.2/"

CopywriteR(sample.control=sample.control,destination.folder=file.path(outPath),reference.folder=file.path(refPath, "hg38_50kb_chr"),bp.param=bp.param,keep.intermediary.files=T)

plotCNA(destination.folder=file.path(outPath))

file.rename(from=paste(outPath,"/CNAprofiles/segment.Rdata",sep=""),to=paste(outPath,"/CNAprofiles/segment.default.Rdata",sep=""))

file.rename(from=paste(outPath,"/CNAprofiles/plots",sep=""),to=paste(outPath,"/CNAprofiles/plots.default",sep=""))

plotCNA(destination.folder=file.path(outPath),alpha=0.1)

library(DNAcopy)

load(paste(outPath,"CNAprofiles/segment.default.Rdata",sep=""))
segs.p <- segments.p(segment.CNA.object)
samples<-unique(segs.p[,1])

samples<-cbind(samples,samples)
samples[,2]<-gsub("log2.","",samples[,2])
samples[,2]<-gsub(".recalibrated.final.bam","",samples[,2])
samples


for (i in 1:nrow(samples))
{
x<-segs.p[segs.p[,1]==as.character(samples[i,1]),]
colnames(x)[1]<-"sample"
x[,1]<-samples[i,2]
write.table(x,file=paste(outPath,"/CNAprofiles/",as.character(samples[i,2]),".seg.txt",sep=""),row.names=F,sep="\t")
}


load(paste(outPath,"CNAprofiles/segment.Rdata",sep=""))
segs.p <- segments.p(segment.CNA.object)
samples<-unique(segs.p[,1])

samples<-cbind(samples,samples)
samples[,2]<-gsub("log2.","",samples[,2])
samples[,2]<-gsub(".recalibrated.final.bam","",samples[,2])
samples


for (i in 1:nrow(samples))
{
x<-segs.p[segs.p[,1]==as.character(samples[i,1]),]
colnames(x)[1]<-"sample"
x[,1]<-samples[i,2]
write.table(x,file=paste(outPath,"/CNAprofiles/",as.character(samples[i,2]),".alpha.0.1.seg.txt",sep=""),row.names=F,sep="\t")
}



################ CHANGE FOLDER NAME  ######################


file.rename(from=paste(outPath,"/CNAprofiles",sep=""),to=paste(outPath,"/CNAprofiles_vsTGEN-PBMC-13_50kb_rerun",sep=""))

#####################CallSegments##########################

setwd(paste(outPath,"/CNAprofiles_vsTGEN-PBMC-13_50kb", sep=""))
 library(DNAcopy)
library(reshape2)


getindex <- function(x) (1:length(x))[x]
##
getmean <- function(tcons, copy = 3, center = 0) log2(tcons * (copy - 2) + 2) + center -1
##
getcopy <- function(mn, tcons, center = 0) (2^(mn - center + 1) - 2)/tcons + 2
##

lsqfun <- function(param,weight,mnval,guess) sum(weight*(getmean(param[2],guess,param[1])-mnval)^2)

################
conserv.dat <- function(probSegments)
## Gets initial estimates of center and tumor content and then generates a new dataset for least-squares fitting
{
   sampleID=probSegments[1,1]
   param <- estimate.param(probSegments)
   dat2 <-cbind(Sample=sampleID,SegmentID=1:nrow(probSegments),probSegments[,2:6])
   dat2[,8] <- round(getcopy(dat2$seg.mean,param[2],param[1]),0)
   colnames(dat2)[8]<- "guess"
   dat2[dat2[,8]< 1 | dat2[,8]>3,8]<-0
   dat3<-dat2[dat2$chrom<23,]
   pdf(paste(sampleID,".ConsSeg.pdf", sep= ""))
   plot(dat3[,7], dat3[,6], 
	xlab="Segment Log2 Mean", ylab = "Number of SNPs",main = sampleID)
   points(dat3[dat3[,8]==1,7], y=dat3[dat3[,8]==1,6], col="blue")
   points(dat3[dat3[,8]==2,7], y=dat3[dat3[,8]==2,6], col="red")
   points(dat3[dat3[,8]==3,7], y=dat3[dat3[,8]==3,6], col="green")
   dev.off()
   o<- order(dat3[,7])
   dat3 <- dat3[o, ]
   write.table(dat3, file= paste(sampleID,".GuessTable.txt", sep= ""),row.names = FALSE, sep = "\t")
   list(center=param[1], tumor=param[2],is.working=param[4],data=dat3)
}

######################

estimate.param <- function(dat)
{
   autSegments<-dat[dat$chrom<23,]
   o<-order(autSegments[,6])
   sortSegments <- autSegments[o,]
   total.mark<-sum(sortSegments$num.mark)
      x<-0; i<-0
   for (j in 1:nrow(sortSegments)) 
      {i <- j
      x <- x+ sortSegments$num.mark[j]
      if (x>total.mark/2) break
      }
   drift <- sortSegments$seg.mean[i] - (x-total.mark/2)*(sortSegments$seg.mean[i]-sortSegments$seg.mean[i-1])/sortSegments$num.mark[i]
   sortSegments$seg.mean <- sortSegments$seg.mean - drift
   ### calculate p to find local minimum
   p <- 0.5
   err<-abs(sortSegments$seg.mean)
   oldErr<-1; change<-0.05
   for (j  in 1:100)
      {
      p<- p + change
      for (i in 1:nrow(sortSegments)) err[i] <- min(abs(sortSegments$seg.mean[i]), abs(sortSegments$seg.mean[i] - log2(1+p/2)), abs(sortSegments$seg.mean[i] - log2(1-p/2)))
      av.err <- sum(err[] * sortSegments$num.mark)/total.mark
      if (av.err>oldErr) change<--0.3*change
      oldErr<-av.err
      }
   p.localmin <- p
   ### calculate p to find global minimum   change=0.008
   min.p <- 0.55
   min.err <- 1000
   for (k in 1:3)
      {
      p<- min.p - 50 * change
      for (j  in 1:100)
         {
         p<- p + change
         err <- sortSegments$seg.mean
         for (i in 1:nrow(sortSegments)) err[i] <- min(abs(sortSegments$seg.mean[i]), abs(sortSegments$seg.mean[i] - log2(1+p/2)), abs(sortSegments$seg.mean[i] - log2(1-p/2)))
         av.err <- sum(err[] * sortSegments$num.mark)/total.mark
         if (av.err<min.err) 
           {
           min.err<-av.err
           min.p<-p
           }
        }
      p<-min.p
      change=change*0.03
      param <- c(drift, p.localmin, p)
      }
   param[4] <-TRUE
   if (abs(param[2]-param[3])>0.01 | param[3] < 0.1 | param[3]>0.85) param[4] <- FALSE
   param
}

##############
Estimate.copy <- function(dat,probes.col=6,seg.mean.col=7,guess.col=8,paramstart=c(0,0.3),minconc=0.10)
# estimates center, concentration, and copy number based on data set of initial guesses
# dat is the segment data
# probes.col,segmean.col and guess.col are the indices of the Probes, segment mean, and copy guess columns
# paramstart is a starting guess for the center and concentration estimates
# minimum allowed concentration
{
   dat <- dat[,c(probes.col,seg.mean.col,guess.col)]
   dat1 <- dat[dat[,3]!=0,]
   try <- nlminb(paramstart,lsqfun,weight=dat1[,1],mnval=dat1[,2],guess=dat1[,3],lower=c(-Inf,minconc))
   center <- try$par[1]
   conc <- try$par[2]
   rawcopy <- getcopy(dat[,2],conc,center)
   copy <- round(rawcopy)
   copy[copy<0] <- 0
   copy[copy>4] <- 4
   #cat("\n",conc)
   list(cons=try$par[2],center=try$par[1],rawcopy=rawcopy,copy=copy)   
}

#############

collapseCGH <-function(dat,chrom=3,strt=4,end=5,probes=6,mn=7, raw=8,copy=9)
{   #combines consecutive segments with the same copy number
   # dat is the data set, the rest are indices of the columns that have that statistic
   
   dat <- dat[!is.na(dat[,copy]),]
   dat <- dat[dat[,copy]>=0,] 
   ndat <- dim(dat)[1]
   chg <- 100*dat[,chrom]+dat[,copy] # trick to compare chr and copy together
   chg <- chg[-1]!=chg[-ndat] # are neighbors different in chr or copy?
   chg <- c(TRUE,chg)
   setx <- rev(getindex(!chg)) #reverse set of row numbers for segments without change
   for(i in setx)
   {    # cat(i,)
      dat[i-1,end] <- dat[i,end] #copies the last "end" of a group of unchanged to the changed above
      nt <- dat[i,probes]+dat[i-1,probes]

      for(j in c(mn)) #This allows use of a different column if mn is not column 7
      {      dat[i-1,j] <- (dat[i-1,j]*dat[i-1,probes]+dat[i,j]*dat[i,probes])/nt
      }
      dat[i-1,probes] <- nt
   }
   dat[chg,] #includes only changed segments
}


#########
dropsmall <-function(dat,chrom=3,strt=4,end=5,probes=6,mn=7,raw=8,copy=9,minnum=4)
#eliminates segments of length less than minnum
{	set <- getindex(dat[,probes]<minnum)
	while(length(set)>0)
	{	#cat("\n nshort",length(set),)
		if(is.element(1,set))
		{	dat[1,copy] <- dat[2,copy]
		}
		n <- dim(dat)[1]
		if(is.element(n,set))
		{	dat[n,copy] <- dat[n-1,copy]
		}
1
		
		c1 <- abs(dat[set,mn]-dat[set-1,mn])+100*abs(dat[set,chrom]-dat[set-1,chrom])	
		c2 <- abs(dat[set,mn]-dat[set+1,mn])+100*abs(dat[set,chrom]-dat[set+1,chrom])

		set1 <- set[c1<c2]
		set2 <- set[c2<=c1]
		set1 <- set1[!is.element(set1-1,set2)]
		dat[set1,copy] <- dat[set1-1,copy]
		dat[set2,copy] <- dat[set2+1,copy]
		dat <- collapseCGH(dat,chrom,strt,end,probes,mn,raw,copy)
		#dat <- dat[,-dim(dat)[2]]
		set <- getindex(dat[,probes]<minnum)
	}
	dat
}	





###############
ConsSegfiles <- list.files(pattern = ".seg.txt")
ConsSegfiles <-ConsSegfiles [!grepl("alpha",ConsSegfiles)]
#ConsSegfiles <-ConsSegfiles [grepl("none",ConsSegfiles)]

c<-cbind(ConsSegfiles,ConsSegfiles)
c[,2]<-gsub(".seg.txt","",c[,2])
c<-cbind(c,NA,NA)
a<-colsplit(c[,2],".vs.",c("f1","f2"))
c<-cbind(c,a)
c<-c[c$f1!=c$f2,]
ConsSegfiles<-as.character(c[,1])


#qiang Method:for 20% tumor

combo<-NULL
for(i in 1:length(ConsSegfiles))
{
name=ConsSegfiles[i]
name<-gsub(".txt","",name)
dat <- read.table (ConsSegfiles[i], header = TRUE,sep="\t",stringsAsFactors=F)
dat<-cbind(dat,2)
colnames(dat)[11]<-"copy"
dat[which(dat$seg.mean >= log2(2.2 / 2)),'copy'] <- '3'
dat[which(dat$seg.mean <= log2(1.8 / 2)),'copy'] <- '1'
dat[which(dat$seg.mean >= log2(3.6 / 2)),'copy'] <- '4'
dat[which(dat$seg.mean <= log2(0.35 / 2)),'copy'] <- '0'

dat[,11]<-as.numeric(dat[,11])
dat<-collapseCGH(dat,chrom=2,strt=3,end=4,probes=5,mn=6,raw=11,copy=11)
dat <- dropsmall(dat,chrom=2,strt=3,end=4,probes=5,mn=6,raw=11,copy=11,minnum=4)
dat<-cbind(dat[,c(1:6,11)],dat$copy-2)
colnames(dat)[8]<- "ivg.copy"
write.table(dat,file=paste(name,".simple.thresh.qiang.setting.seg",sep=""),row.names=F,sep="\t",quote=F)
}

combo<-NULL
for(i in 1:length(ConsSegfiles))
{
name=ConsSegfiles[i]
name<-gsub(".txt","",name)
dat <- read.table (ConsSegfiles[i], header = TRUE,sep="\t",stringsAsFactors=F)
dat<-cbind(dat,2)
colnames(dat)[11]<-"copy"
dat[which(dat$seg.mean >= log2(2.4 / 2)),'copy'] <- '3'
dat[which(dat$seg.mean <= log2(1.6/ 2)),'copy'] <- '1'
dat[which(dat$seg.mean >= log2(3.6 / 2)),'copy'] <- '4'
dat[which(dat$seg.mean <= log2(0.35 / 2)),'copy'] <- '0'
dat[,11]<-as.numeric(dat[,11])
dat<-collapseCGH(dat,chrom=2,strt=3,end=4,probes=5,mn=6,raw=11,copy=11)
dat <- dropsmall(dat,chrom=2,strt=3,end=4,probes=5,mn=6,raw=11,copy=11,minnum=4)
dat<-cbind(dat[,c(1:6,11)],dat$copy-2)
colnames(dat)[8]<- "ivg.copy"
write.table(dat,file=paste(name,".2.4-1.6.thresh.seg",sep=""),row.names=F,sep="\t",quote=F)
}


combo<-NULL
for(i in 1:length(ConsSegfiles))
{
name=ConsSegfiles[i]
name<-gsub(".txt","",name)
dat <- read.table (ConsSegfiles[i], header = TRUE,sep="\t",stringsAsFactors=F)
dat<-cbind(dat,2)
colnames(dat)[11]<-"copy"
dat[which(dat$seg.mean >= log2(2.6 / 2)),'copy'] <- '3'
dat[which(dat$seg.mean <= log2(1.4/ 2)),'copy'] <- '1'
dat[which(dat$seg.mean >= log2(3.6 / 2)),'copy'] <- '4'
dat[which(dat$seg.mean <= log2(0.35 / 2)),'copy'] <- '0'
dat[,11]<-as.numeric(dat[,11])
dat<-collapseCGH(dat,chrom=2,strt=3,end=4,probes=5,mn=6,raw=11,copy=11)
dat <- dropsmall(dat,chrom=2,strt=3,end=4,probes=5,mn=6,raw=11,copy=11,minnum=4)

dat<-cbind(dat[,c(1:6,11)],dat$copy-2)
colnames(dat)[8]<- "ivg.copy"
write.table(dat,file=paste(name,".2.6-1.4.thresh.seg",sep=""),row.names=F,sep="\t",quote=F)
}


param.list <- NULL

######MCK METHOD###########

 ConsSegfiles <- list.files(pattern = ".seg.txt")
ConsSegfiles <-ConsSegfiles [!grepl("alpha",ConsSegfiles)]
#ConsSegfiles <-ConsSegfiles [grepl("none",ConsSegfiles)]

c<-cbind(ConsSegfiles,ConsSegfiles)
c[,2]<-gsub(".seg.txt","",c[,2])
c<-cbind(c,NA,NA)
a<-colsplit(c[,2],".vs.",c("f1","f2"))
c<-cbind(c,a)
c<-c[c$f1!=c$f2,]
ConsSegfiles<-as.character(c[,1])

param.list <- NULL

#

i=1

print(i)
dat <- read.table (ConsSegfiles[i], header = TRUE)
    sampleID <- as.character(ConsSegfiles[i])
    sampleID <- sub(".seg.txt","",sampleID)
sampleID <- gsub("_trimmed","",sampleID)
   LiberalName <- as.character(ConsSegfiles[i])
    LiberalName <- sub(".seg.txt",".alpha.0.1.seg.txt", LiberalName)
   dat1 <- conserv.dat(dat)
## Get better parameters
   param <- Estimate.copy(dat1$data)
   param.list<- rbind(param.list,c("sampleID"=sampleID, "1st guess fraction"=dat1$tumor, "1st guess center"=dat1$center, "tumor fraction" = param$cons, "center"= param$center))
   write.table(param.list, file= paste("ParameterList",Sys.Date(), "txt", sep="."), row.names=FALSE, sep="\t")
   dat<- read.table(LiberalName, header = TRUE)
dat$chrom<-gsub("chrX","chr23",dat$chrom)
dat$chrom<-gsub("chrY","chr24",dat$chrom)
dat$chrom<-as.numeric(gsub("chr","",dat$chrom))
dat[,1]<-sampleID
      dat <-cbind(Sample=sampleID,SegmentID=1:nrow(dat),dat[,2:6])
   dat[,8] <- getcopy(dat$seg.mean,param$cons,param$center)
   dat[,9]<-2
   colnames(dat)[8:9] <- c("raw.copy","copy")
     dat[dat[,8]< 1.4,9]<-1
   dat[dat[,8]> 2.6,9]<-3
   dat[dat[,8]< 0.35,9]<-0
   dat[dat[,8]> 3.6,9]<-4
   dat <-collapseCGH(dat)
   dat <- dropsmall(dat)
   dat[,2]<- 1:dim(dat)[1]
   write.table(dat, file= paste(sampleID,".CombinedSeg.txt", sep= ""), row.names=FALSE, sep="\t")
   pdf(paste(sampleID, ".CombinedSeg.pdf", sep= ""))
   plot(dat[,7], dat[,6], 
	xlab="Segment Log2 Mean", ylab = "Number of SNPs",main = sampleID)
   points(dat[dat[,9]==1,7], y=dat[dat[,9]==1,6], col="blue")
   points(dat[dat[,9]==2,7], y=dat[dat[,9]==2,6], col="red")
   points(dat[dat[,9]==3,7], y=dat[dat[,9]==3,6], col="green")
   dev.off()

i=2

print(i)
dat <- read.table (ConsSegfiles[i], header = TRUE)
    sampleID <- as.character(ConsSegfiles[i])
    sampleID <- sub(".seg.txt","",sampleID)
sampleID <- gsub("_trimmed","",sampleID)
   LiberalName <- as.character(ConsSegfiles[i])
    LiberalName <- sub(".seg.txt",".alpha.0.1.seg.txt", LiberalName)
   dat1 <- conserv.dat(dat)
## Get better parameters
   param <- Estimate.copy(dat1$data)
   param.list<- rbind(param.list,c("sampleID"=sampleID, "1st guess fraction"=dat1$tumor, "1st guess center"=dat1$center, "tumor fraction" = param$cons, "center"= param$center))
   write.table(param.list, file= paste("ParameterList",Sys.Date(), "txt", sep="."), row.names=FALSE, sep="\t")
   dat<- read.table(LiberalName, header = TRUE)
dat$chrom<-gsub("chrX","chr23",dat$chrom)
dat$chrom<-gsub("chrY","chr24",dat$chrom)
dat$chrom<-as.numeric(gsub("chr","",dat$chrom))
dat[,1]<-sampleID
      dat <-cbind(Sample=sampleID,SegmentID=1:nrow(dat),dat[,2:6])
   dat[,8] <- getcopy(dat$seg.mean,param$cons,param$center)
   dat[,9]<-2
   colnames(dat)[8:9] <- c("raw.copy","copy")
     dat[dat[,8]< 1.4,9]<-1
   dat[dat[,8]> 2.6,9]<-3
   dat[dat[,8]< 0.35,9]<-0
   dat[dat[,8]> 3.6,9]<-4
   dat <-collapseCGH(dat)
   dat <- dropsmall(dat)
   dat[,2]<- 1:dim(dat)[1]
   write.table(dat, file= paste(sampleID,".CombinedSeg.txt", sep= ""), row.names=FALSE, sep="\t")
   pdf(paste(sampleID, ".CombinedSeg.pdf", sep= ""))
   plot(dat[,7], dat[,6], 
	xlab="Segment Log2 Mean", ylab = "Number of SNPs",main = sampleID)
   points(dat[dat[,9]==1,7], y=dat[dat[,9]==1,6], col="blue")
   points(dat[dat[,9]==2,7], y=dat[dat[,9]==2,6], col="red")
   points(dat[dat[,9]==3,7], y=dat[dat[,9]==3,6], col="green")
   dev.off()


########
i=3

print(i)
dat <- read.table (ConsSegfiles[i], header = TRUE)
    sampleID <- as.character(ConsSegfiles[i])
    sampleID <- sub(".seg.txt","",sampleID)
sampleID <- gsub("_trimmed","",sampleID)
   LiberalName <- as.character(ConsSegfiles[i])
    LiberalName <- sub(".seg.txt",".alpha.0.1.seg.txt", LiberalName)
   dat1 <- conserv.dat(dat)
## Get better parameters
   param <- Estimate.copy(dat1$data)
   param.list<- rbind(param.list,c("sampleID"=sampleID, "1st guess fraction"=dat1$tumor, "1st guess center"=dat1$center, "tumor fraction" = param$cons, "center"= param$center))
   write.table(param.list, file= paste("ParameterList",Sys.Date(), "txt", sep="."), row.names=FALSE, sep="\t")
   dat<- read.table(LiberalName, header = TRUE)
dat$chrom<-gsub("chrX","chr23",dat$chrom)
dat$chrom<-gsub("chrY","chr24",dat$chrom)
dat$chrom<-as.numeric(gsub("chr","",dat$chrom))
dat[,1]<-sampleID
      dat <-cbind(Sample=sampleID,SegmentID=1:nrow(dat),dat[,2:6])
   dat[,8] <- getcopy(dat$seg.mean,param$cons,param$center)
   dat[,9]<-2
   colnames(dat)[8:9] <- c("raw.copy","copy")
     dat[dat[,8]< 1.4,9]<-1
   dat[dat[,8]> 2.6,9]<-3
   dat[dat[,8]< 0.35,9]<-0
   dat[dat[,8]> 3.6,9]<-4
   dat <-collapseCGH(dat)
   dat <- dropsmall(dat)
   dat[,2]<- 1:dim(dat)[1]
   write.table(dat, file= paste(sampleID,".CombinedSeg.txt", sep= ""), row.names=FALSE, sep="\t")
   pdf(paste(sampleID, ".CombinedSeg.pdf", sep= ""))
   plot(dat[,7], dat[,6], 
	xlab="Segment Log2 Mean", ylab = "Number of SNPs",main = sampleID)
   points(dat[dat[,9]==1,7], y=dat[dat[,9]==1,6], col="blue")
   points(dat[dat[,9]==2,7], y=dat[dat[,9]==2,6], col="red")
   points(dat[dat[,9]==3,7], y=dat[dat[,9]==3,6], col="green")
   dev.off()

i=4

print(i)
dat <- read.table (ConsSegfiles[i], header = TRUE)
    sampleID <- as.character(ConsSegfiles[i])
    sampleID <- sub(".seg.txt","",sampleID)
sampleID <- gsub("_trimmed","",sampleID)
   LiberalName <- as.character(ConsSegfiles[i])
    LiberalName <- sub(".seg.txt",".alpha.0.1.seg.txt", LiberalName)
   dat1 <- conserv.dat(dat)
## Get better parameters
   param <- Estimate.copy(dat1$data)
   param.list<- rbind(param.list,c("sampleID"=sampleID, "1st guess fraction"=dat1$tumor, "1st guess center"=dat1$center, "tumor fraction" = param$cons, "center"= param$center))
   write.table(param.list, file= paste("ParameterList",Sys.Date(), "txt", sep="."), row.names=FALSE, sep="\t")
   dat<- read.table(LiberalName, header = TRUE)
dat$chrom<-gsub("chrX","chr23",dat$chrom)
dat$chrom<-gsub("chrY","chr24",dat$chrom)
dat$chrom<-as.numeric(gsub("chr","",dat$chrom))
dat[,1]<-sampleID
      dat <-cbind(Sample=sampleID,SegmentID=1:nrow(dat),dat[,2:6])
   dat[,8] <- getcopy(dat$seg.mean,param$cons,param$center)
   dat[,9]<-2
   colnames(dat)[8:9] <- c("raw.copy","copy")
     dat[dat[,8]< 1.4,9]<-1
   dat[dat[,8]> 2.6,9]<-3
   dat[dat[,8]< 0.35,9]<-0
   dat[dat[,8]> 3.6,9]<-4
   dat <-collapseCGH(dat)
   dat <- dropsmall(dat)
   dat[,2]<- 1:dim(dat)[1]
   write.table(dat, file= paste(sampleID,".CombinedSeg.txt", sep= ""), row.names=FALSE, sep="\t")
   pdf(paste(sampleID, ".CombinedSeg.pdf", sep= ""))
   plot(dat[,7], dat[,6], 
	xlab="Segment Log2 Mean", ylab = "Number of SNPs",main = sampleID)
   points(dat[dat[,9]==1,7], y=dat[dat[,9]==1,6], col="blue")
   points(dat[dat[,9]==2,7], y=dat[dat[,9]==2,6], col="red")
   points(dat[dat[,9]==3,7], y=dat[dat[,9]==3,6], col="green")
   dev.off()

