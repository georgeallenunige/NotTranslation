


#load R packages
library(riboWaltz)
library(GenomicAlignments)
library(GenomicFeatures)
library(RiboProfiling)
library(Biostrings)
library(rtracklayer)
library(openxlsx)



#Get reads and perform the P-site offset...
deletion="NOT5"
bamfolder="/path/to/bams/"
setwd(bamfolder)

bams=c("WT_1.bam",
"WT_2.bam",
"Deletion_1.bam",
"Deletion_2.bam")


annotation_file=read.table("./annotation_file_RiboWaltz.xls",sep="\t",header=T,quote="")

reads_list <- bamtolist(bamfolder = bamfolder, annotation = annotation_file)
reads_list=reads_list[match(gsub(".bam","",bams),names(reads_list))]

#Filter by read length
filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom", length_filter_vector = 28:31)

#Perform P-site offset
psite_offset <- psite(filtered_list, extremity = "auto")
reads_psite_listO <- psite_info(filtered_list, psite_offset,granges=TRUE)


reads_psite_list=reads_psite_listO
for(i in 1:(length(bams))){
start(reads_psite_list[[i]])=reads_psite_list[[i]]$psite
end(reads_psite_list[[i]])=reads_psite_list[[i]]$psite
}


#pscod is the total P-site counts in the coding sequence
pscod=c()
for(i in 1:(length(bams))){
co=countOverlaps(grc,reads_psite_list[[i]])
pscod=c(pscod,sum(co))
}



numerize=function(tab){
    outTab=c()
    for(i in 1:ncol(tab)){
    outTab=cbind(outTab,as.numeric(as.character(tab[,i])))
    }
    colnames(outTab)=colnames(tab)
    rownames(outTab)=rownames(tab)
    return(as.matrix(outTab))
}


#read in rpkms
rpkm=numerize(read.table("./RPKMs_WT_NOT4_Not5_MinRPKM1_NoDubiousOrfs.txt",sep="\t",header=T))
rpkm=cbind(rpkm[,1:2],rpkm[,grep(deletion,colnames(rpkm))])


#Get the list of RP Genes
rpgs=read.table("./rpgs.xls",sep="\t",header=T)[,4]
idgrc=gsub(".+_","",as.character(rownames(rpkm)))
rpIDs=rownames(rpkm)[match(rpgs,idgrc)]


####read in the CDS locations and create a granges object with them
cds=read.table("./CDSLocsWithWidths.xls",sep="\t",header=T,quote="")
grco <- GRanges(seqnames = cds[,1], ranges = IRanges(start = cds[,2], width = cds[,3]-cds[,2]+1))


grc=grco
#retain all transcripts with mean log2 rpkm > 5
fiveUp=rpkm[which(log(rowMeans(rpkm),2)>5),]
grc=grc[which(as.character(seqnames(grc))%in%rownames(fiveUp)==T)]

#use this if only RP transcripts are to be used
#grc=grc[which(as.character(seqnames(grc))%in%rpIDs==F)]


#retain all transcripts with at least 100 bases
subG=grc[which(width(grc)>=100)]
outPDF=paste0(deletion,"MetageneNormGitHubTest_ScaledCDSPlot.pdf")
subG=subG[which(duplicated(as.character(seqnames(subG)))==F)]
grct=subG

#split all CDS ranges into 1 base windows
grct1=unlist(slidingWindows(grct, width=1L, step=1L))

#For each 1 base window establish which integer quantile of the CDS it falls into (i.e. which of the 100 equal bins it falls into)
starts=start(grct)[match(as.character(seqnames(grct1)),as.character(seqnames(grct)))]
widths=width(grct)[match(as.character(seqnames(grct1)),as.character(seqnames(grct)))]
relLoc=start(grct1)-starts
quant=floor(100*relLoc/widths)+1

#assign integer quantiles to each base
widl=width(grct)
values(grct1)=quant


#Get CDS counts to use in normalisation and take the per-base mean
coTabN=c()
for(i in 1:length(bams)){
co=countOverlaps(subG,reads_psite_list[[i]])
coTabN=cbind(coTabN,co)
}

rownames(coTabN)=as.character(seqnames(subG))
normTab=sweep(coTabN,1,width(subG),FUN="/")


colnames(normTab)=bams

wids=width(subG)
normTab=cbind(normTab,wids)


plotTab=c()
for(i in 1:length(bams)){

#get the normalising factor for each base (same for all within a CDS)
norms=normTab[match(as.character(seqnames(grct1)),rownames(normTab)),]
#get counts for each bin and scale by the normalising factor
co=countOverlaps(grct1,reads_psite_list[[i]])/norms[,i]

bylis=grct1$X
#get the per-base mean for each CDS quantile - averaged over all CDSs
pdep=aggregate(co,by=list(bylis),FUN=mean)
plotTab=cbind(plotTab,pdep[,2])

}

maxy=max(plotTab[,1:4])


nots=c("NOT4","NOT5")
cols=c("chartreuse3","firebrick1")
delCol=cols[match(deletion,nots)]

#Plot the mean of WT and deletion (normalised)
pdf(outPDF)
scaler=1.2
scalerl=1.2
par(cex.axis=1*scaler*1.2)
par(cex.lab=1*scalerl*1.2)
plot(0:99,(plotTab[,1]+plotTab[,2])/2,col="black",type="l",ylim=c(0,maxy),xlab="Scaled position, whole CDS (100 equal bins)",ylab="Normalised RPF Counts ",lwd=4,xaxt="n",bty="n")
axis(1, at=0:5*20, labels=c("START",as.character(1:4*20),"STOP"))

abline(v=0,lwd=7,col="lightgrey")
abline(v=99,lwd=7,col="lightgrey")

points(0:99,(plotTab[,1]+plotTab[,2])/2,col="black",type="l",ylim=c(0,maxy),lwd=4)
points(0:99,(plotTab[,3]+plotTab[,4])/2,col=delCol,type="l",lwd=4)

if(deletion=="NOT4"){
legend("bottom",col=c("black",delCol),legend=c("WT",expression(paste("",italic("not4"),Delta,""))),lwd=10,cex=1.5,bty="n")
}
if(deletion=="NOT5"){
legend("bottom",col=c("black",delCol),legend=c("WT",expression(paste("",italic("not5"),Delta,""))),lwd=10,cex=1.5,bty="n")
}
dev.off()


#########



plotTab=c()
for(i in 1:length(bams)){

#Get counts overlapping bins and normalise by total counts over all CDSs (millions)
co=countOverlaps(grct1,reads_psite_list[[i]])/(pscod[i]/1000000)

bylis=grct1$X
pdep=aggregate(co,by=list(bylis),FUN=mean)

plotTab=cbind(plotTab,pdep[,2])

}

maxy=max(plotTab[,1:4])


#Plot the mean of WT and deletion expression weighted
outPDF2=gsub("Norm","NoNorm",outPDF)
pdf(outPDF2)
scaler=1.2
scalerl=1.2
par(cex.axis=1*scaler*1.2)
par(cex.lab=1*scalerl*1.2)
plot(0:99,(plotTab[,1]+plotTab[,2])/2,col="black",type="l",ylim=c(0,maxy),xlab="Scaled position, whole CDS (100 equal bins)",ylab="RPF Counts - expression weighted ",lwd=4,xaxt="n",bty="n")
axis(1, at=0:5*20, labels=c("START",as.character(1:4*20),"STOP"))
abline(v=0,lwd=7,col="lightgrey")
abline(v=99,lwd=7,col="lightgrey")
points(0:99,(plotTab[,1]+plotTab[,2])/2,col="black",type="l",ylim=c(0,maxy),lwd=3)

points(0:99,(plotTab[,3]+plotTab[,4])/2,col="firebrick1",type="l",lwd=4)
legend("bottom",col=c("black","firebrick1"),legend=c("WT",expression(paste("",italic("not5"),Delta,""))),lwd=10,cex=1.5,bty="n")
dev.off()



