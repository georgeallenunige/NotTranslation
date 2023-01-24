
# This plots metagenes scaled into 100 equal bins split by expression and lengths


5pSeqDepths="path_to_5'PSeqDepths"
wtRPKMs="path_to_WTrpkms"
base="path_to_outputfolder"

ct=c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot")
# colz= is colours associated with each strain

# 5'P-Seq depths is a list for each strain/replicate and within each of these lists is a list of lists of transcripts 
# giving 5'P-Seq counts for each position on their CDS (only kept if CDS length > 100)
cds=readRDS(5pSeqDepths)

# pick out the strains you want
grabs=c()
for(cac in 1:length(ct)){
grabs=c(grabs,grep(paste0("^",ct[cac]),names(cl)))
}

kp=cl[grabs]
cat("\n\n")
print(names(kp))
kpo=kp

# function to separate each CDS into 100 bins
hun=function(vv){
return(ceiling(100*(1:length(vv))/length(vv)))
}

# get wt rpkms
wtAI=read.table(wtRPKMs)
wtAI=wtAI[which(wtAI>0)]
wtAI=wtAI[order(wtAI)]

# get CDS lengths
len=unlist(lapply(cl[[1]],length))
names(len)=names(cl[[1]])
len=log(1+len,2)
len=len[order(len)]

# split expression and length in to low, medium, high
expL=list(names(wtAI[1:floor(length(wtAI)/3)]),names(wtAI[(ceiling(length(wtAI)/3)):(floor(2*length(wtAI)/3))]),names(wtAI[(ceiling(2*length(wtAI)/3)):length(wtAI)]))
lenL=list(names(len[1:floor(length(len)/3)]),names(len[(ceiling(length(len)/3)):(floor(2*length(len)/3))]),names(len[(ceiling(2*length(len)/3)):length(len)]))
cc=c("Low","Med","High")



pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"SplitByExp_Len_SCALED.pdf"),width=14)
par(mfrow=c(3,3))


meanLisAdd=list()
for(ii in 1:3){
for(jj in 1:3){


# select depths for the expression/length group in question 
picks=intersect(expL[[ii]],lenL[[jj]])

# cumulate counts for each of 100 bins over CDS
scal=list()
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
nums=unlist(lapply(cds,hun))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[1]]=ag

for(k in 2:length(kp)){
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag
}

names(scal)=names(kp)

# normalise to library size
adf=as.data.frame(scal)
adf=adf[,2*(1:(ncol(adf)/2))]
af=t(sweep(adf,STATS=colSums(adf)/1000000,MARGIN=2,FUN="/"))

# average for each strain
meanLisS=list()
for(state in 1:length(ct)){
meanLisS[[state]]=colMeans(af[grep(paste0("^",ct[state]),rownames(af)),])
}

meanLisAdd[[jj]]=meanLisS

miny=min(unlist(meanLisS))
maxy=max(unlist(meanLisS))


for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(1:100,meanLisS[[ci]],ylim=c(miny,maxy),col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P Depth",bty="n",paste0(cc[ii]," Expr, ",cc[jj]," Length,"))
} else {
print(ci)
points(1:100,meanLisS[[ci]],col=colz[ci],type="l",lwd=2)
}
}

if(ii*jj==9){
legend("topleft",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1)
abline(v=0,lty=2)
abline(v=100,lty=2)
}


}}

dev.off()





















library(openxlsx)


setwd("~/Desktop/Not5_5P/DESeq2/xlsx/")
lf=list.files()
for(l in 1:length(lf)){
xls=read.table(lf[l],sep="\t",header=TRUE)

xlsx=cbind(gsub("_.+","",rownames(xls)),gsub(".+_","",rownames(xls)),rownames(xls),xls)
colnames(xlsx)[1:3]=c("Symbol","ID","Symbol_ID")

write.xlsx(xlsx,file=paste0(lf[l],"x"))

}


library(riboWaltz)
library(GenomicAlignments)
library(GenomicFeatures)
library(RiboProfiling)
library(Biostrings)

library(rtracklayer)

heads=readLines("~/Desktop/Not5/Saccharomyces_cerevisiae.R64-1-1.cds.all.header.txt")
geneId=gsub(" .+","",gsub(".+gene:","",heads))
dubious=rep("good",length(heads))
#remove dubious orfs and some putative proteins



bads=unique(c(grep("YNL042W-B|YOR072W|YPR108W-A|YOL013W-A",heads),grep("dubious|Dubious|YLR162W-A|YLR154W-C|YLR162W|YLR154C-G",heads)))
dubious[bads]="dubious"

dubTab=cbind(geneId,dubious)
dubs=dubTab[which(dubTab[,2]=="dubious"),1]


####
cds=read.table("~/Desktop/Not5/CDSLocsWithWidths.xls",sep="\t",header=T,quote="")
transOut <- readDNAStringSet('~/Desktop/Not5/CDSPlusUTRs.fa')
wids=width(transOut)
nams=names(transOut)

grc <- GRanges(seqnames = cds[,1], ranges = IRanges(start = cds[,2], width = cds[,3]-cds[,2]+1))
idgrc=gsub(".+_","",as.character(seqnames(grc)))
grc=grc[which(idgrc%in%dubs==FALSE)]


fa = FaFile('~/Desktop/Not5/CDSPlusUTRs.fa')
idx = scanFaIndex(fa)

cseq = getSeq(fa, unlist(grc))

grcoo=grc[grep("^ATG",as.character(cseq),invert=F)]








#####






#write.table(allTab,paste0(bam,".ScaledPSiteDepthPlotCDS.xls"),sep="\t",quote=F,row.names=F)

library(GenomicAlignments)
library(GenomicFeatures)
library(RiboProfiling)
library(Biostrings)




trans=readDNAStringSet('~/Desktop/Not5/CDSPlusUTRs.fa')
cds=read.table("~/Desktop/Not5/CDSLocs.xls",header=T)
#trans=trans[match(cds$names,


grc <- GRanges(seqnames = cds[,1], ranges = IRanges(start = cds[,2], width = cds[,3]-cds[,2]+1))

cdstr=trans[grc]
dna=as.character(cdstr)
codons=as.character(translate(cdstr, genetic.code=GENETIC_CODE, if.fuzzy.codon="error"))

cdstr=cdstr[grep("^M",codons)]
dna=dna[grep("^M",codons)]
grc=grc[grep("^M",codons)]
cds=cds[grep("^M",codons),]
codons=codons[grep("^M",codons)]


grc1=unlist(slidingWindows(grc, width=1L, step=1L))
grc3=unlist(slidingWindows(grc, width=3L, step=3L))



###


dnaString=strsplit(paste(dna,collapse=""),"")[[1]]
codString=strsplit(paste(codons,collapse=""),"")[[1]]

f0=(1:(length(dnaString)/3))*3-2
aaString=paste0(dnaString[f0],dnaString[f0+1],dnaString[f0+2])

aaOut=rep("ERROR",length(aaString)*3)
aaOut[f0]=aaString
aaOut[f0+1]=aaString
aaOut[f0+2]=aaString

codOut=rep("ERROR",length(codString)*3)
codOut[f0]=codString
codOut[f0+1]=codString
codOut[f0+2]=codString

tabC=cbind(dnaString,aaOut,codOut)
colnames(tabC)=c("DNA","Codon","AA")


###


AA_cod=paste(tabC[,3],tabC[,2],sep="_")
AA_cod=AA_cod[(1:(length(AA_cod)/3))*3-2]



grc3=unlist(slidingWindows(grc, width=3L, step=3L))
values(grc3) <- DataFrame(AA_cod = AA_cod)


grc=grc[which(as.character(seqnames(grc))%in%as.character(seqnames(grcoo))==T)]
grc3=grc3[which(as.character(seqnames(grc3))%in%as.character(seqnames(grcoo))==T)]


baseFreq=table(grc3$AA_cod)





rps=as.character(read.csv("~/Desktop/UNIGE_TOOLS/GRCh38.84.TranMatchCDS/RPS.csv")[,2])
rpl=as.character(read.csv("~/Desktop/UNIGE_TOOLS/GRCh38.84.TranMatchCDS/RPL.csv")[,2])

rpgo=c(rps,rpl)

idgrc=gsub("_.+","",as.character(seqnames(grc)))
rpg=grc[na.omit(match(rpgo,idgrc))]

rp=grc[as.integer(match(as.character(seqnames(rpg)),as.character(seqnames(grc))))]




grc3o=grc3




#####First 20%


#~/Downloads/fivepseq_results02/fivepseq_frg_tot/fivepseq_counts/frg_WTtot3_09_0423_GTGATC/protein_coding/data_summary.txt
#
#sum(unlist(lapply(cl$frg_WTtot3_09_0423_GTGATC,sum)))

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))
#clo=cl

for(wido in c(100,200,300)){
#for(wido in c(100)){
cl=clo
utr=0


lens=unlist(lapply(cl[[1]],length))
if(utr==0){
lens=lens-200
}
for(ll in 1:length(cl)){

cl[[ll]]=cl[[ll]][which(lens>(wido+100))]
}


#hunF=function(v){
#return(sum(v[101:200]))
#}
#
#hunF_L=function(lis){
#return(lapply(lis,hunF))
#}
#
#cl100F=lapply(cl,hunF_L)


######
######


hunF=function(v){
return(sum(v[101:(wido+100)-utr]))
}

hunF_L=function(lis){
return(lapply(lis,hunF))
}

cl100F=lapply(cl,hunF_L)


hunL=function(v){
en=(length(v)-100+utr)
return(sum(v[(en-(wido-1)):en]))
}

hunL_L=function(lis){
return(lapply(lis,hunL))
}

cl100L=lapply(cl,hunL_L)


df100F=c()
for(d in 1:length(cl100F)){
df100F=cbind(df100F,cl100F[[d]])
}

colnames(df100F)=names(cl100F)

df100L=c()
for(d in 1:length(cl100L)){
df100L=cbind(df100L,cl100L[[d]])
}

colnames(df100L)=names(cl100L)


numerize=function(tab){
    outTab=c()
    for(i in 1:ncol(tab)){
        outTab=cbind(outTab,as.numeric(as.character(tab[,i])))
    }
    colnames(outTab)=colnames(tab)
    rownames(outTab)=rownames(tab)
    return(as.matrix(outTab))
}

df100F=numerize(df100F)
df100L=numerize(df100L)

summary(df100F)

mf=apply(df100F,1,max)
ml=apply(df100L,1,max)

mm=apply(cbind(mf,ml),1,max)



df100F2=df100F[which(mm>9),]
df100L2=df100L[which(mm>9),]

dim(df100F2)
dim(df100L2)

fcFL=log((1+df100F2)/(1+df100L2),2)

fiveP=fcFL[,grep("frg_",colnames(fcFL),invert=TRUE)]
fiveP=cbind(fiveP[,grep("_AI_",colnames(fiveP),invert=TRUE)],fiveP[,grep("_AI_",colnames(fiveP))])


ctLis=list(c("WTsol","WTtot","not5D_sol","not5D_tot"),c("WTsol","WTtot","not4D_sol","not4D_tot"),
c("Wt_AI_sol","Wt_AI_tot","Not1p_AI_sol","Not1p_AI_tot"),c("Wt_AI_sol","Wt_AI_tot","Not4p_AI_sol","Not4p_AI_tot"),c("Wt_AI_sol","Wt_AI_tot","Not5p_AI_sol","Not5p_AI_tot"))


pdf(paste0("~/Desktop/Not5_5P/boxplotFirst",wido,"OverLast",wido,"with",utr,"nt_ofUTRs_23DEC2020.pdf"))
boxplot(fiveP,outline=FALSE,ylim=c(-8,6),names=rep("",ncol(fiveP)),ylab=paste0("logFC First ",wido," bases / Last Including ",utr," bp up/downstream of CDS"))
cn=colnames(fiveP)
cn=gsub("_0.+","",cn)
text(1:ncol(fiveP)+1,-5,cn,col="darkcyan",srt=90,pos=2)
dev.off()



for(cl in 1:length(ctLis)){

cn=colnames(fiveP)
cto=ctLis[[cl]]
locs=c()
for(st in 1:length(cto)){
locs=c(locs,grep(cto[st],cn))
}

outNam=cn[locs]
outLis=fiveP[,locs]


pdf(paste0("~/Desktop/Not5_5P/boxplotFirst",paste(cto,collapse="_"),wido,"OverLast",wido,"with",utr,"nt_ofUTRs_23DEC2020.pdf"))
boxplot(outLis,outline=FALSE,ylim=c(-8,6),names=rep("",ncol(outLis)),ylab=paste0("logFC First ",wido," bases / Last Including ",utr," bp up/downstream of CDS"))
#cn=colnames(fiveP)
cn=gsub("_0.+","",outNam)
text(1:ncol(outLis)+0.25,-5,cn,col="darkcyan",srt=90,pos=2)
dev.off()

}

}








head(cl100[[1]])

cg=cl
into=intersect(names(cl[[1]]),as.character(seqnames(grc)))

for(gg in 1:length(cg)){
cg[[gg]]=cl[[gg]][match(into,names(cl[[gg]]))]
cg[[gg]]=lapply(cg[[gg]],keep)
}


w1=unlist(lapply(cg[[1]],length))

grg=grc[match(into,seqnames(grc))]
wr=width(grg)
names(wr)=seqnames(grg)

into=names(wr[which(wr==w1)])

####check that the lengths are right then filter out the 5 that don't match


for(gg in 1:length(cg)){
cg[[gg]]=cg[[gg]][match(into,names(cg[[gg]]))]
}

grg=grc[match(into,seqnames(grc))]


cga=lapply(cg,unlist)

len=length(cga[[11]])/3

cgm=list()
for(gg in 1:length(cga)){

cgm[[gg]]=cga[[gg]][(1:len)*3-2]+cga[[gg]][(1:len)*3-1]+cga[[gg]][(1:len)*3]
}

names(cgm)=names(cga)

g3g=grc3[which(seqnames(grc3)%in%seqnames(grg)==TRUE)]

codo=width(grg)/3

tar=table(seqnames(grc3))
ta3=tar[match(into,names(tar))]


setwd("~/Downloads/fivepseq_results02/")
lf=list.files()
lf=grep("fivepseq_log",lf,invert=TRUE,value=TRUE)
#lf=grep("not5|WT",lf,value=TRUE)[c(3,4,1,2)]


tra=read.table("~/Downloads/fivepseq_results02/fivepseq_WTsol/fivepseq_counts/WTsol1_04_0318_CGTGCA/protein_coding/transcript_assembly.txt",header=TRUE)

rct=gsub("transcript:|_mRNA","",tra[,1])
names=rtracklayer::import("~/Desktop/Not5/Saccharomyces_cerevisiae.R64-1-1.84.gtf")

nam=names$gene_name[match(rct,names$gene_id)]
nam[which(is.na(nam))]=rct[which(is.na(nam))]

pn=paste(nam,rct,sep="_")


namz=c()
coLis=list()
bas=0

twen=function(v){
return(as.numeric(v[1:floor(length(v)/2)]))
#return(as.numeric(v[(length(v)-floor(length(v)/2)+1):length(v)]))
}

for(l in 1:length(lf)){
#for(l in 1:length(1:2)){

setwd(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts"))
lf2=list.files()

for(l2 in 1:length(lf2)){

fl=readLines(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts/",lf2[l2],"/protein_coding/counts_FULL_LENGTH.txt"))

namz=c(namz,lf2[l2])

st=strsplit(fl,split="\t")
names(st)=pn

bas=bas+1


coLis[[bas]]=lapply(st,twen)

}}

names(coLis)=namz


saveRDS(coLis,file="~/Desktop/Not5_5P/CDSCountsFirst50Percent.rds")


######





cl20=readRDS(file="~/Desktop/Not5_5P/CDSCountsFirst50Percent.rds")
names(cl20)=gsub("�","D_",names(cl20))


cl80=readRDS(file="~/Desktop/Not5_5P/CDSCountsLast50Percent.rds")
names(cl80)=gsub("�","D_",names(cl80))


ctab20=c()
for(lit in 1:length(cl20)){
ctab20=cbind(ctab20,unlist(lapply(cl20[[lit]],sum)))
}
colnames(ctab20)=names(cl20)


ctab80=c()
for(lit in 1:length(cl80)){
ctab80=cbind(ctab80,unlist(lapply(cl80[[lit]],sum)))
}
colnames(ctab80)=names(cl80)


cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))


ctab=c()
for(lit in 1:length(cl)){
ctab=cbind(ctab,unlist(lapply(cl[[lit]],sum)))
}
colnames(ctab)=names(cl)


write.table(ctab20,file="~/Desktop/Not5/FivePCountsFirst50Percent.xls",sep="\t")
write.table(ctab80,file="~/Desktop/Not5/FivePCountsLast50Percent.xls",sep="\t")
write.table(ctab,file="~/Desktop/Not5/FivePCountsAll100Percent.xls",sep="\t")


ctaba=read.table("~/Desktop/Not5/FivePCountsAll100Percent.xls",sep="\t",header=TRUE)


setwd("~/Downloads/align_dedup_SC_SP_mapped/")

bms=list.files(recursive=TRUE)

bm=grep("bam$",bms,value=TRUE)
bm=bm[-1]

system("rm ~/Desktop/Not5_5P/mapped_SC_SP.txt")
system("touch ~/Desktop/Not5_5P/mapped_SC_SP.txt")

for(b in 1:length(bm)){

cmd=paste0("samtools flagstat ~/Downloads/align_dedup_SC_SP_mapped/",bm[b]," | grep mapped | grep 100 >> ~/Desktop/Not5_5P/mapped_SC_SP.txt")
system(cmd)

}

tex=read.csv("~/Desktop/Not5_5P/mapped_SC_SP.txt",sep=" ",header=FALSE)
rownames(tex)=gsub(".fastqAligned.primary.out_dedup.bam","",gsub("_SC_SP_","_",gsub("�","D",gsub(".+/","",bm))))
rownames(tex)=toupper(gsub("03","3",gsub("02","2",gsub("01","1",rownames(tex)))))

baser=gsub(".fastqAligned.primary.out_dedup.bam","",gsub("_SC_SP_.+","",gsub("�","D",gsub(".+/","",bm))))
baser=toupper(gsub("03","3",gsub("02","2",gsub("01","1",baser))))
baser=gsub("WT_TOT","WTTOT",baser)
baser=gsub("WT_SOL","WTSOL",baser)

tex=cbind(baser,tex)
tex=tex[,1:2]

#cst=colSums(ctaba)
#names(cst)=toupper(names(cst))
#names(cst)=gsub("NOTP_AI_TOT","NOT1P_AI_TOT",names(cst))

lisf=function(liz){
return(sum(unlist(lapply(liz,sum))))
}

cst=unlist(lapply(cl,lisf))
names(cst)=toupper(names(cst))
names(cst)=gsub("NOTP_AI_TOT","NOT1P_AI_TOT",names(cst))

#528190

~/Downloads/align_dedup_SC_SP_mapped/align_dedup_SC_SP_bam/01_0615/Not1p_AI_tot01_SC_SP_TAATGT.fastqAligned.primary.out_dedup.bam


library(GenomicAlignments)
b1=readGAlignments("~/Downloads/align_dedup_SC_SP_mapped/align_dedup_SC_SP_bam/01_0615/Not1p_AI_tot01_SC_SP_TAATGT.fastqAligned.primary.out_dedup.bam")

names=rtracklayer::import("~/Downloads/align_dedup_SC_SP_mapped/SC_Spombe/gff/SC_SP.gff")
cds=names[which(names$type=="CDS")]








library(GenomicAlignments)
b1=readGAlignments("~/Downloads/align_dedup_SC_SP_mapped/align_dedup_SC_SP_bam/01_0615/Not1p_AI_tot01_SC_SP_TAATGT.fastqAligned.primary.out_dedup.bam")


co=sum(countOverlaps(cds,b1))


nam=names$gene_name[match(rct,names$gene_id)]
nam[which(is.na(nam))]=rct[which(is.na(nam))]






nam=names$gene_name[match(rownames(coTab),names$gene_id)]
nam[which(is.na(nam))]=rownames(coTab)[which(is.na(nam))]



#FRG_WTTOT3

sum(unlist(lapply(cl$frg_WTtot3_09_0423_GTGATC,sum)))

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))



#tab=cbind(names(cst),gsub("SOL1


num=paste0(gsub("[0-9]_.+","",names(cst)),gsub("_.+","",gsub(".+SOL|.+TOT","",names(cst))))


grep("AI_TOT",gsub("[0-9]_.+","",names(cst)),value=TRUE)

cstab=cbind(names(cst),num,cst)


match(cstab[,2],tex[,1])
cstab[,1][is.na(match(cstab[,2],tex[,1]))]


binc=cbind(cstab,tex[match(cstab[,2],tex[,1]),])


plot(log(as.numeric(binc[,3]),2),log(as.numeric(binc[,5]),2),xlim=c(0,21),ylim=c(0,21),pch=16,cex=1,xlab="log CDS counts SC",ylab="log CDS counts SP - spike")

cti=cbind(ctab20[,1],ctab80[,1],ctab[,1])

colnames(cti)=c("20","80","all")


frg_Not1p_AI_sol1_10_0324_GACTAG


cd=read.table("~/Downloads/fivepseq_results02/fivepseq_frg_AI_sol/fivepseq_counts/frg_Not1p_AI_sol1_10_0324_GACTAG/protein_coding/count_distribution.txt")


sum(cd[,2])


library(GenomicAlignments)
b1=readGAlignments("~/Downloads/align_dedup_SC_SP_mapped/align_dedup_SC_SP_bam/04_0318/Not1p_AI_sol01_SC_SP_ACGCTG.fastqAligned.primary.out_dedup.bam")


basi="~/Downloads/align_dedup_SC_SP_mapped/align_dedup_SC_SP_bam/"
setwd(basi)




write.table(cbi,"~/Desktop/Not5_5P/POMBE_CDS_Counts.xls",sep="\t")

allz=cbind(binc,cbi[match(binc[,2],cbi[,2]),])

sum(unlist(lapply(cl$frg_WTtot3_09_0423_GTGATC,sum)))

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))

names(cl)=gsub("Notp_AI_tot","Not1p_AI_tot",names(cl))

cdf=function(v){
remov=c(101:(length(v)-100))
return(as.numeric(v[remov]))
}

cl2=cl
sumoz=c()
for(ca in 1:length(cl)){
cl2[[ca]]=lapply(cl[[ca]],cdf)
sumoz=c(sumoz,sum(unlist(lapply(cl2[[ca]],sum))))
}


trunk=toupper(gsub("3_.+","",(gsub("2_.+","",gsub("1_.+","",names(cl))))))
numiz=gsub("_.+","",gsub(".+tot","",gsub(".+sol","",names(cl))))

biz=paste0(trunk,numiz)


sc=cbind(names(cl),biz,sumoz)

write.table(sc,"~/Desktop/Not5_5P/CEREVISIAE_CDS_Counts.xls",sep="\t")


cbind(cbi,sc[match(cbi[,2],sc[,2]),])

sc_sp=cbind(sc,cbi[match(sc[,2],cbi[,2]),])

colnames(sc_sp)=c("Nam_SC","Nam2_SC","Counts_SC","Nam_SP","Nam2_SP","Counts_SP")
write.table(sc_sp,"~/Desktop/Not5_5P/CEREVISIAE_AND_POMBE_CDS_Counts.xls",sep="\t")


100-100*as.numeric(sc_sp[,3])/as.numeric(sc_sp[,6])

cbi[,2][is.na(match(cbi[,2],sc[,2]))]

######here



tra=read.table("~/Downloads/fivepseq_results02/fivepseq_WTsol/fivepseq_counts/WTsol1_04_0318_CGTGCA/protein_coding/transcript_assembly.txt",header=TRUE)

rct=gsub("transcript:|_mRNA","",tra[,1])
names=rtracklayer::import("~/Downloads/align_dedup_SC_SP_mapped/SC_Spombe/gff/SC_SP.gff")

cds=names[which(names$type=="CDS")]
ids=gsub("CDS:|_mRNA","",cds$ID)
cds=cds[which(ids%in%rct==TRUE)]

#####

setwd(basi)
coss=c()
ded=list.files(pattern="bam$",recursive=TRUE)
ded=ded[-1]

nt=names(tb)
tib=c()
for(b in 1:length(ded)){

library(GenomicAlignments)
b1=readGAlignments(paste0(basi,ded[b]))


tb=table(seqnames(b1))
tib=cbind(tib,tb[match(nt,names(tb))])
#co=countOverlaps(cds,b1)
#coss=c(coss,sum(co))

}



ded1=gsub(".fastqAligned.primary.out_dedup.bam","",gsub("_SC_SP_","_",gsub("�","D",gsub(".+/","",ded))))
ded1=toupper(gsub("03","3",gsub("02","2",gsub("01","1",ded1))))

ded2=gsub(".fastqAligned.primary.out_dedup.bam","",gsub("_SC_SP_.+","",gsub("�","D",gsub(".+/","",ded))))
ded2=toupper(gsub("03","3",gsub("02","2",gsub("01","1",ded2))))
ded2=gsub("WT_TOT","WTTOT",baser)
ded2=gsub("WT_SOL","WTSOL",baser)


cbi=cbind(ded1,ded2,coss)

write.table(cbi,"~/Desktop/Not5_5P/POMBE_CDS_Counts.xls",sep="\t")


cb3=t(tib)
rownames(cb3)=ded2

write.table(cb3,"~/Desktop/Not5_5P/POMBE_CDS_CountsMapByChr.xls",sep="\t")

cb3=read.table("~/Desktop/Not5_5P/POMBE_CDS_CountsMapByChr.xls",sep="\t",header=TRUE,row.names=3)

cb3=cb3[,-c(1:2)]
sumSP=rowSums(cb3[,1:4])
sumSC=rowSums(cb3[,5:ncol(cb3)])

tk=100*sumSP/(sumSP+sumSC)

to=tk[order(names(tk))]

so=round(sumSC[order(names(sumSC))]/1000000,2)

plot(to,so,cex=0.1,xlab="Spike-In %",ylab="Total SC Counts")
text(to,so,names(so),cex=0.3)

rp5a=read.table("~/Desktop/Not5_5P/AllDataTypes/Counts_RiboSeq_PolyASome_FivePSeq.xls",sep="\t",header=TRUE)[,1:6]

bina=bin[match(rownames(rp5a),rownames(bin)),]





riboStart=read.table("~/Desktop/Not5/CountsFirst50Percent.xls",sep="\t",header=TRUE)
riboEnd=read.table("~/Desktop/Not5/CountsLast50Percent.xls",sep="\t",header=TRUE)


fivePStart=read.table("~/Desktop/Not5/FivePCountsFirst50Percent.xls",sep="\t",header=TRUE)
fivePEnd=read.table("~/Desktop/Not5/FivePCountsLast50Percent.xls",sep="\t",header=TRUE)



#fivePStart=cbind(fivePStart[,grep("^WTsol",colnames(fivePStart))],fivePStart[,grep("^not5D_sol",colnames(fivePStart))])
#fivePEnd=cbind(fivePEnd[,grep("^WTsol",colnames(fivePEnd))],fivePEnd[,grep("^not5D_sol",colnames(fivePEnd))])

fivePStart=cbind(fivePStart[,grep("^WTsol",colnames(fivePStart))],fivePStart[,grep("^not5D_tot",colnames(fivePStart))])
fivePEnd=cbind(fivePEnd[,grep("^WTsol",colnames(fivePEnd))],fivePEnd[,grep("^not5D_tot",colnames(fivePEnd))])


riboStart=cbind(riboStart[,grep("WT",colnames(riboStart))],riboStart[,grep("NOT5",colnames(riboStart))])
riboEnd=cbind(riboEnd[,grep("WT",colnames(riboEnd))],riboEnd[,grep("NOT5",colnames(riboEnd))])

into=intersect(rownames(fivePStart),rownames(riboStart))

riboStart=riboStart[match(into,rownames(riboStart)),]
riboEnd=riboEnd[match(into,rownames(riboEnd)),]

fivePStart=fivePStart[match(into,rownames(fivePStart)),]
fivePEnd=fivePEnd[match(into,rownames(fivePEnd)),]

finTab=cbind(rowMeans(riboStart[,1:2]),rowMeans(riboEnd[,1:2]),rowMeans(riboStart[,1:2+2]),rowMeans(riboEnd[,1:2+2]),rowMeans(fivePStart[,1:3]),rowMeans(fivePEnd[,1:3]),rowMeans(fivePStart[,1:2+3]),rowMeans(fivePEnd[,1:2+3]))

colnames(finTab)=c("RiboWTStart","RiboWTEnd","RiboNOT5Start","RiboNOT5End","FivePWTStart","FivePWTEnd","FivePNOT5Start","FivePNOT5End")





cts=cbind(log(1+rowMeans(finTab[,1:4]),2),log(1+rowMeans(finTab[,1:4+4]),2))

pdf("~/Desktop/Not5_5P/RPF_vs_FivePCountsTot.pdf")
plot(cts[,1],cts[,2],pch=16,xlab="Log Mean RPF Counts",ylab="Log Mean FiveP Tot Counts")

text(10,2,paste0("R=",round(cor(cts[,1],cts[,2]),2)),cex=1.5)

dev.off()


exp=apply(finTab[,1:4+4],1,min)
hi=finTab[which(exp>median(exp)),]

fc=log((1+hi[,c(2,4,6,8)])/(1+hi[,c(2,4,6,8)-1]),2)
colnames(fc)=paste(colnames(fc),"_Over_Start")

dif=cbind(fc[,2]-fc[,1],fc[,4]-fc[,3])


pdf("~/Desktop/Not5_5P/RPF_vs_FivePTotCounts_End_Over_Start.pdf")

plot(dif[,1],dif[,2],pch=16,xlab="LogFC End/Start RPF Counts not5D-WT",ylab="LogFC End/Start fiveP Tot Counts not5D-WT")

abline(v=0,lty=2,col="red")
abline(h=0,lty=2,col="red")

text(1.5,-2,paste0("R=",round(cor(dif[,1],dif[,2]),2)),cex=1.5)

dev.off()


bin=na.omit(cbind(coti,rp5a[match(rownames(coti),rownames(rp5a)),1:6])













lcpm=log(1+rowMeans(cbind(coTab20[,1:2],coTab80[,1:2])),2)


n5k=cbind(coTab20[which(lcpm>median(lcpm)),c(1,2,6,7)],coTab80[which(lcpm>median(lcpm)),c(1,2,6,7)])
dif=log((1+rowMeans(n5k[,7:8]))/(1+rowMeans(n5k[,7:8-4])),2)-log((1+rowMeans(n5k[,7:8-2]))/(1+rowMeans(n5k[,7:8-4-2])),2)

#dif=na.omit(dif)


n5p=cbind(ctab20[,grep("^WTsol|^not5D_sol",colnames(ctab20))],ctab80[,grep("^WTsol|^not5D_sol",colnames(ctab80))])
n5p=n5p[match(names(dif),rownames(n5p)),]

dif5p=log((1+rowMeans(n5p[,10:12-3]))/(1+rowMeans(n5p[,10:12-3-6])),2)-log((1+rowMeans(n5p[,10:12]))/(1+rowMeans(n5p[,10:12-6])),2)



plot(log(1+rowMeans(n5p[,10:12-3]),2),log(1+rowMeans(n5k[,7:8]),2))



n5p=cbind(ctab20[,grep("^WTsol|^not5D_sol",colnames(ctab20))],ctab80[,grep("^WTsol|^not5D_sol",colnames(ctab80))])

n5pk=n5p[match(rownames(n5k),rownames(n5p)),]

r5=log(1+rowMeans(n5k),2)
r5p=log(1+rowMeans(n5pk),2)

plot(r5,r5p)












cd5p=ctab[,grep("^WTsol|^not5D_sol",colnames(ctab))]



coTab20=read.table("~/Desktop/Not5/CountsFirst20Percent.xls",sep="\t",header=TRUE)
coTab80=read.table("~/Desktop/Not5/CountsFirst80Percent.xls",sep="\t",header=TRUE)

c28=coTab20+coTab80
c28=c28[,grep("WT|NOT5",colnames(c28))]


p28=na.omit(cbind(log(1+rowMeans(c28),2),log(1+rowMeans(cd5p[match(rownames(c28),rownames(cd5p)),]),2)))

plot(p28[,1],p28[,2])




rp5=read.table("~/Desktop/Not5_5P/AllDataTypes/Counts_RiboSeq_PolyASome_FivePSeq.xls",sep="\t",header=TRUE)
cor(rp5,use="complete.obs")[,1]


rp5a=rp5[,16:18]

wtc=ctab[,grep("^WTsol",colnames(ctab))]
wtc=wtc[match(rownames(rp5),rownames(wtc)),]
#w20=ctab20[,grep("^WTsol",colnames(ctab20))]


cot5=c28[na.omit(match(rownames(rp5),rownames(c28))),]





coTab20=read.table("~/Desktop/Not5/CountsFirst20Percent.xls",sep="\t",header=TRUE)

wtc=ctab[,grep("^WTsol",colnames(ctab))]
wtc=wtc[match(rownames(coTab20),rownames(wtc)),]


cot=cbind(coTab20[,1],wtc[,1])

cor(cot,use="complete.obs")


setwd("~/Downloads/fivepseq_results02/")
lf=list.files()
lf=grep("fivepseq_log",lf,invert=TRUE,value=TRUE)
#lf=grep("not5|WT",lf,value=TRUE)[c(3,4,1,2)]


tra=read.table("~/Downloads/fivepseq_results02/fivepseq_WTsol/fivepseq_counts/WTsol1_04_0318_CGTGCA/protein_coding/transcript_assembly.txt",header=TRUE)

rct=gsub("transcript:|_mRNA","",tra[,1])
names=rtracklayer::import("~/Desktop/Not5/Saccharomyces_cerevisiae.R64-1-1.84.gtf")

nam=names$gene_name[match(rct,names$gene_id)]
nam[which(is.na(nam))]=rct[which(is.na(nam))]

pn=paste(nam,rct,sep="_")


namz=c()
coLis=list()
bas=0

for(l in 3:length(lf)){
#for(l in 1:length(1:2)){

setwd(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts"))
lf2=list.files()

for(l2 in 1:length(lf2)){

fl=readLines(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts/",lf2[l2],"/protein_coding/counts_FULL_LENGTH.txt"))

namz=c(namz,lf2[l2])

st=strsplit(fl,split="\t")
names(st)=pn

bas=bas+1

coLis[[bas]]=lapply(st,as.numeric)

}}

names(coLis)=namz


saveRDS(coLis,file="~/Desktop/Not5_5P/CDSCounts.rds")




cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))
#
#"Wt_AI_sol1_04_0318_GACTAG"        "Wt_AI_sol2_05_0219_ACGCTG"
#"Not1p_AI_sol1_04_0318_ACGCTG"     "Not1p_AI_sol2_05_0219_GTGATC"     "Not1p_AI_sol3_06_0120_ACTGGT"
#
#"Wt_AI_tot1_01_0615_CGTGCA"        "Wt_AI_tot2_02_0516_TAATGT"        "Wt_AI_tot3_03_0417_CTACAC"
#"Not1p_AI_tot1_01_0615_TAATGT"     "Not1p_AI_tot3_03_0417_GACTAG"     "Notp_AI_tot2_02_0516_CTACAC"
#


coTab=read.table("~/Desktop/Not5_5P/countTab_All_5PSeq_Samples.xls",sep="\t")
coTab=coTab[rev(order(rowMeans(coTab))),]
coTab=coTab[,grep("FRG",colnames(coTab),invert=TRUE)]

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))

names(cl)=toupper(names(cl))
#cl=cl[grep("WT_AI_SOL2_05_0219_ACGCTG",names(cl),invert=TRUE)]
cl=cl[match(colnames(coTab),names(cl))]


coTab[1,1]
wan=cl[[1]]
#wan[[grep("HSC82_YMR186W",names(wan))]]
sum(wan[[grep("HSC82_YMR186W",names(wan))]])


contos=colSums(coTab)
names(contos)=names(cl)


for(geni in c("DLD1"

gn=list()

for(lio in 1:length(cl)){

gn[[lio]]=as.numeric(unlist(cl[[lio]][grep(geni,names(cl[[lio]]))]))/((contos/(10^6))[lio])

}

names(gn)=names(cl)


tab=as.data.frame(gn)
colnames(tab)=gsub("1_.+|2_.+|3_.+","",colnames(tab))

uni=unique(colnames(tab))

tiby=c()
for(u in 1:length(uni)){
tiby=cbind(tiby,rowMeans(tab[,which(colnames(tab)==uni[u])]))
}

colnames(tiby)=uni



pdf(paste0("~/Desktop/Not5_5P/",geni,".pdf"),width=(10*nrow(tiby)/1000)/2,height=14/2)

par(mfrow=c(3,1))

for(mi in 1:14){

plot(tiby[,mi],type="l",ylab="Counts per million",main=colnames(tiby)[mi])

}
dev.off()

}





sum(unlist(cl[[1]]))

#ct=c("Wt_AI_sol","Not1p_AI_sol")
#ct=c("WTsol","not5D_sol")

ct=c("WTtot","not5D_tot","not4D_tot")

#ctLis=list(c("WTsol","not5D_sol","not4D_sol","Wt_AI_sol","Not1p_AI_sol"),c("WTtot","not5D_tot","not4D_tot","Wt_AI_tot","Not1p_AI_tot"))

ctLis=list(c("WTsol","not5D_sol","not4D_sol"),c("Wt_AI_sol","Not1p_AI_sol"),c("Wt_AI_sol","Not1p_AI_sol","Not4p_AI_sol","Not5p_AI_sol"),c("WTtot","not5D_tot","not4D_tot"),c("Wt_AI_tot","Not1p_AI_tot"),c("Wt_AI_tot","Not1p_AI_tot","Not4p_AI_tot","Not5p_AI_tot"))



ctLis=list(c("WTsol","not5D_sol","not4D_sol"),c("Wt_AI_sol","Not1p_AI_sol"),c("Wt_AI_sol","Not1p_AI_sol","Not4p_AI_sol","Not5p_AI_sol"),c("WTtot","not5D_tot","not4D_tot"),c("Wt_AI_tot","Not1p_AI_tot"),c("Wt_AI_tot","Not1p_AI_tot","Not4p_AI_tot","Not5p_AI_tot"))


#ctLis=list(c("WTsol","WTtot","not5D_sol","not5D_tot"),
#c("WTsol","WTtot","not4D_sol","not4D_tot"),
#c("Wt_AI_sol","Wt_AI_tot","Not1p_AI_sol","Not1p_AI_tot"),
#c("Wt_AI_sol","Wt_AI_tot","Not5p_AI_sol","Not5p_AI_tot"),
#c("Wt_AI_sol","Wt_AI_tot","Not4p_AI_sol","Not4p_AI_tot"))


#colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")
samp=c("WTsol","Wt_AI_sol","WTtot","Wt_AI_tot",
"not5D_sol","Not5p_AI_sol","not5D_tot","Not5p_AI_tot",
"not4D_sol","Not4p_AI_sol","not4D_tot","Not4p_AI_tot",
"Not1p_AI_sol","Not1p_AI_tot")

#colo=c("black","black","darkgrey","darkgrey",
#"firebrick3","firebrick3","darkred","darkred",
#"chartreuse3","chartreuse3","darkgreen","darkgreen",
#"dodgerblue","darkblue")

colo=c("black","black","black","black",
"firebrick3","firebrick3","firebrick3","firebrick3",
"chartreuse3","chartreuse3","chartreuse3","chartreuse3",
"dodgerblue","dodgerblue")

colz=colo[match(ct,samp)]

for (ccc in 1:length(ctLis)) {
#for (ccc in 1) {

ct=ctLis[[ccc]]

grabs=c()
for(cac in 1:length(ct)){
grabs=c(grabs,grep(paste0("^",ct[cac]),names(cl)))
}

kp=cl[grabs]

tots=function(lif){
return(sum(unlist(lif)))
}

cdsTot=lapply(kp,tots)

sta=function(v){
return(as.numeric(v[1:200]))
}

tLis=list()
for(i in 1:length(kp)){
tLis[[i]]=t(as.data.frame(lapply(kp[[i]],sta)))
}

names(tLis)=names(kp)

mnt=lapply(tLis,colSums)
weigh=unlist(cdsTot)/1000000

adf=t(as.data.frame(mnt))
sw=sweep(adf,MARGIN=1,STATS=weigh,FUN="/")

meanLis=list()
for(state in 1:length(ct)){
meanLis[[state]]=colMeans(sw[grep(paste0("^",ct[state]),rownames(sw)),])
}

names(meanLis)=ct

colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")


maxy=max(unlist(lapply(meanLis,max)))

pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"START_DEC23_2020.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(-100:99,meanLis[[ci]],col=colz[ci],ylim=c(0,maxy),type="l",lwd=2,xlab="Distance from START",ylab="Normalised 5P Depth",bty="n")
} else {
print(ci)
points(-100:99,meanLis[[ci]],col=colz[ci],type="l",lwd=2)
}
}

legend("topright",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)
abline(v=0,lty=2)

dev.off()


colz=colo[match(ct,samp)]


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"START_Shift17_DEC23_2020.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(-100:99+17,meanLis[[ci]],col=colz[ci],ylim=c(0,maxy),type="l",lwd=2,xlab="Distance from START - Shift 17",ylab="Normalised 5P Depth",bty="n")
} else {
print(ci)
points(-100:99+17,meanLis[[ci]],col=colz[ci],type="l",lwd=2)
}
}

abline(v=0,lty=2)
legend("topright",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)

dev.off()




#######Stops


sto=function(v){
return(as.numeric(v[(length(v)-199):length(v)]))
}


tLis=list()
for(i in 1:length(kp)){
tLis[[i]]=t(as.data.frame(lapply(kp[[i]],sto)))
}

names(tLis)=names(kp)

mnt=lapply(tLis,colSums)
weigh=unlist(cdsTot)/1000000

adf=t(as.data.frame(mnt))
sw=sweep(adf,MARGIN=1,STATS=weigh,FUN="/")

meanLis=list()
for(state in 1:length(ct)){
meanLis[[state]]=colMeans(sw[grep(paste0("^",ct[state]),rownames(sw)),])
}

names(meanLis)=ct

#colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")

maxy=max(unlist(lapply(meanLis,max)))


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"STOP_DEC23_2020.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(-100:99,meanLis[[ci]],col=colz[ci],ylim=c(0,maxy),type="l",lwd=2,xlab="Distance from STOP",ylab="Normalised 5P Depth",bty="n")
} else {
print(ci)
points(-100:99,meanLis[[ci]],col=colz[ci],type="l",lwd=2)
}
}

legend("topright",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)
abline(v=0,lty=2)

dev.off()


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"STOP_Shift17_DEC23_2020.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(-100:99+17,meanLis[[ci]],col=colz[ci],ylim=c(0,maxy),type="l",lwd=2,xlab="Distance from STOP - Shift 17",ylab="Normalised 5P Depth",bty="n")
} else {
print(ci)
points(-100:99+17,meanLis[[ci]],col=colz[ci],type="l",lwd=2)
}
}

abline(v=0,lty=2)
legend("topright",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)

dev.off()



}















#########
#########









library(openxslx)

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")

names(cl)=gsub("�","D_",names(cl))





at=read.table("~/Desktop/Not5_5P/ParCLIP_WithRIBO_PELECH_FirstLast100_LastOverCDS_RPKM_LogFC_AddSoluble_PerMRNA_BindingA_DEGRONSONLY_January14_2021_SPIKE_IN.xls",sep="\t",header=TRUE)

pcl=na.omit(at$frg_Not1p_AI_sol_Over_frg_Not1p_AI_tot_NormToWT)
pcl=pcl[!is.infinite(pcl)]

lopar=na.omit(at[which(at$frg_Not1p_AI_sol_Over_frg_Not1p_AI_tot_NormToWT<(-0.5)),])
hipar=na.omit(at[which(at$frg_Not1p_AI_sol_Over_frg_Not1p_AI_tot_NormToWT>(0.5)),])


upList=list(rownames(lopar)[which(lopar$frg_Not4p_AI_sol_Over_frg_Not4p_AI_tot_NormToWT<(-0.5))],
rownames(lopar)[which(lopar$frg_Not4p_AI_sol_Over_frg_Not4p_AI_tot_NormToWT>(0.5))],
rownames(hipar)[which(hipar$frg_Not4p_AI_sol_Over_frg_Not4p_AI_tot_NormToWT<(-0.5))],
rownames(hipar)[which(hipar$frg_Not4p_AI_sol_Over_frg_Not4p_AI_tot_NormToWT>(0.5))])

lapply(upList,length)

nots=c("LowNot1Sol_LowSOL_NOT4P_vs_WT_FRG_SPIKE_IN","LowNot1Sol_HighSOL_NOT4P_vs_WT_FRG_SPIKE_IN","HighNot1Sol_LowSOL_NOT4P_vs_WT_FRG_SPIKE_IN","HighNot1Sol_HighSOL_NOT4P_vs_WT_FRG_SPIKE_IN")



#at=read.table("~/Desktop/Not5_5P/ParCLIP_WithRIBO_PELECH_FirstLast100_LastOverCDS_RPKM_LogFC_AddSoluble_PerMRNA_BindingA_DEGRONS_January13_2021.xls",sep="\t",header=TRUE)
#
##edgeR_DiffExp_RPKM__FRG_NOT4P_AI_SOL_Over_FRG_NOT4P_AI_TOT
##edgeR_DiffExp_RPKM__FRG_NOT4P_AI_SOL_Over_FRG_WT_AI_SOL
#
#pcl=na.omit(at$logFC_NOT4D_Over_PolII_PARCLIP_CDS)
#pcl=pcl[!is.infinite(pcl)]
#
#lopar=na.omit(at[which(at$logFC_NOT4D_Over_PolII_PARCLIP_CDS<summary(pcl)[2]),])
#hipar=na.omit(at[which(at$logFC_NOT4D_Over_PolII_PARCLIP_CDS>summary(pcl)[5]),])
#
#upList=list(rownames(lopar)[which(lopar$edgeR_DiffExp_RPKM__FRG_NOT4P_AI_TOT_Over_FRG_WT_AI_TOT<(-0.5))],
#rownames(lopar)[which(lopar$edgeR_DiffExp_RPKM__FRG_NOT4P_AI_TOT_Over_FRG_WT_AI_TOT>(0.5))],
#rownames(hipar)[which(hipar$edgeR_DiffExp_RPKM__FRG_NOT4P_AI_TOT_Over_FRG_WT_AI_TOT<(-0.5))],
#rownames(hipar)[which(hipar$edgeR_DiffExp_RPKM__FRG_NOT4P_AI_TOT_Over_FRG_WT_AI_TOT>(0.5))])
#
#lapply(upList,length)
#
#nots=c("LowPARCLIP_Low_NOT4P_Over_WT_FRGTOT","LowPARCLIP_High_NOT4P_Over_WT_FRGTOT","HighPARCLIP_Low_NOT4P_Over_WT_FRGTOT","HighPARCLIP_High_NOT4P_Over_WT_FRGTOT")
#










#ct=c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol")
ct=c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot")
#ct=c("WTsol","not4D_sol","not5D_sol")
#ct=c("WTtot","not4D_tot","not5D_tot")
#ct=c("WT_AIsol","not4D_sol","not5D_sol")



lizo=list(c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot"),c("WTsol","not4D_sol","not5D_sol"),c("WTtot","not4D_tot","not5D_tot"),c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol"))


for(lo in 1:length(lizo)){
#for(lo in 1){

ct=lizo[[lo]]


grabs=c()
for(cac in 1:length(ct)){
grabs=c(grabs,grep(paste0("^",ct[cac]),names(cl)))
}

kp=cl[grabs]
cat("\n\n")
print(names(kp))
kpo=kp


cdf=function(v){
remov=c(1:100,(length(v)-99):length(v))
return(as.numeric(v[-remov]))
}

hun=function(vv){
return(ceiling(100*(1:length(vv))/length(vv)))
}


#8.5 mil CH, 5883
#10 mil, 7893
#100*57722/66650000


rpkm=read.table("~/Desktop/Not5_5P/rpkmTab_All_5PSeq_Samples.xls",sep="\t",header=TRUE)

match(rownames(rpkm),names(cds))


wtAI=log(1+rowMeans(rpkm[,grep("FRG_WT_AI_TOT",colnames(rpkm))]),2)
wtAI=wtAI[which(wtAI>0)]

len=unlist(lapply(cl[[1]],length))-200
names(len)=names(cl[[1]])
len=log(1+len,2)
len=len[match(names(wtAI),names(len))]


wtAI=wtAI[order(wtAI)]
len=len[order(len)]

expL=list(names(wtAI[1:floor(length(wtAI)/3)]),names(wtAI[(ceiling(length(wtAI)/3)):(floor(2*length(wtAI)/3))]),names(wtAI[(ceiling(2*length(wtAI)/3)):length(wtAI)]))

lenL=list(names(len[1:floor(length(len)/3)]),names(len[(ceiling(length(len)/3)):(floor(2*length(len)/3))]),names(len[(ceiling(2*length(len)/3)):length(len)]))


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"_",nots[1],"SCALE_Wide19Jan2021.pdf"),width=14)
par(mfrow=c(2,2))

meanLisAdd=list()
for(ii in 3){
for(jj in 1:length(upList)){

#picks=intersect(expL[[ii]],lenL[[jj]])
picks=upList[[jj]]

##ag=aggregate(cds[which(names(cds)%in%picks==TRUE)],by=list(nums),FUN=sum)
#
#print(match(picks,names(cds)))
#
#}}


#kp=kpo

scal=list()
k=1
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
c1=lapply(kp[[k]],cdf)
cds=unlist(c1)
nums=unlist(lapply(c1,hun))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag

for(k in 2:length(kp)){
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
cds=unlist(lapply(kp[[k]],cdf))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag
}

names(scal)=names(kp)

adf=as.data.frame(scal)
adf=adf[,2*(1:(ncol(adf)/2))]

af=t(sweep(adf,STATS=colSums(adf)/1000000,MARGIN=2,FUN="/"))

meanLisS=list()
for(state in 1:length(ct)){
meanLisS[[state]]=colMeans(af[grep(paste0("^",ct[state]),rownames(af)),])
}

meanLisAdd[[jj]]=meanLisS

colz=c("black","chartreuse3","firebrick3","dodgerblue","darkorange","purple")

miny=min(unlist(meanLisS))
maxy=max(unlist(meanLisS))


cc=c("Low","Med","High")
#pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"Exp",ii,"Len",jj,"SCALE_DEC2020.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(1:100,meanLisS[[ci]],ylim=c(miny,maxy),col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P Depth",bty="n",main=nots[jj])
#,main=paste(cc[ii],"Expr,",cc[jj],"Length"))
} else {
print(ci)
points(1:100,meanLisS[[ci]],col=colz[ci],type="l",lwd=2)
}
}

if(ii*jj==9){
legend("topleft",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1)
abline(v=0,lty=2)
abline(v=100,lty=2)
}


}}

dev.off()

}






pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"Exp","Len","SCALE_DEC2020_WideXX.pdf"),width=7)


#ylim=c(miny,maxy),

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(1:100,meanLisAdd[[1]][[ci]]+meanLisAdd[[2]][[ci]]+meanLisAdd[[3]][[ci]],col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P Depth",bty="n",main=paste(cc[ii],"Expr, Sum over 3 Length levels"))
} else {
print(ci)
points(1:100,meanLisAdd[[1]][[ci]]+meanLisAdd[[2]][[ci]]+meanLisAdd[[3]][[ci]],col=colz[ci],type="l",lwd=2)
}
}

if(ii*jj==9){
legend("bottomleft",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1)
abline(v=0,lty=2)
abline(v=100,lty=2)
}


dev.off()









########MetageneNPC
#####





library(openxslx)

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))


#ct=c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol")
ct=c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot")


inz=unlist(unlist(lizo))



samp=c("WTsol","Wt_AI_sol","WTtot","Wt_AI_tot",
"not5D_sol","Not5p_AI_sol","not5D_tot","Not5p_AI_tot",
"not4D_sol","Not4p_AI_sol","not4D_tot","Not4p_AI_tot",
"Not1p_AI_sol","Not1p_AI_tot")

colo=c("black","black","black","black",
"firebrick3","firebrick3","firebrick3","firebrick3",
"chartreuse3","chartreuse3","chartreuse3","chartreuse3",
"dodgerblue","dodgerblue")

colz=colo[match(ct,samp)]


# lizo=list(c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot"),c("WTsol","not4D_sol","not5D_sol"),c("WTtot","not4D_tot","not5D_tot"),c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol"))

lizo=list(c("Wt_AI_sol","Not1p_AI_sol","Not4p_AI_sol"),c("Wt_AI_tot","Not1p_AI_tot","Not4p_AI_tot"))


spike=FALSE

for(lo in 1:length(lizo)){
#for(lo in 2){

ct=lizo[[lo]]


grabs=c()
for(cac in 1:length(ct)){
grabs=c(grabs,grep(paste0("^",ct[cac]),names(cl)))
}

kp=cl[grabs]
cat("\n\n")
print(names(kp))
kpo=kp


cdf=function(v){
remov=c(1:100,(length(v)-99):length(v))
return(as.numeric(v[-remov]))
}

hun=function(vv){
return(ceiling(100*(1:length(vv))/length(vv)))
}

fram=function(v){
return(c(sum(v[which((1:length(v))%%3==1)]),sum(v[which((1:length(v))%%3==2)]),sum(v[which((1:length(v))%%3==0)])))
}

picks=solAss[,1]
picks=names(kpo[[1]])


pickLis=list(names(kpo[[1]]),solAss[,1],solEnc[,1],c(solAss[2,1],solAss[2,1]))
cc=c("All_Genes","NPC_Associated","NPC_Encoded","NUP1")


pdf(paste0("~/Desktop/Not5_5P/NPC/Frames_Nov30_2022_NPC_NUP1_Soluble.pdf"),width=12,height=4)
par(mfrow=c(1,5))

for(pp in 1:length(pickLis)){

picks=pickLis[[pp]]

scal=list()
k=1
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
c1=lapply(kp[[k]],cdf)
frams=colSums(t(as.data.frame(lapply(c1,fram))))

framTab=t(as.data.frame(frams))
for(k in 2:length(kpo)){
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
c1=lapply(kp[[k]],cdf)
frams=colSums(t(as.data.frame(lapply(c1,fram))))
framTab=rbind(framTab,frams)
}
rownames(framTab)=names(kpo)

per=100*sweep(framTab,1,rowSums(framTab),FUN="/")
perc=aggregate(per,by=list(gsub("1_0.+|2_0.+|3_0.+","",names(kpo))),FUN=mean)
rownames(perc)=perc[,1]
perc=perc[c(4,2,3)-1,2:4]
rownames(perc)=c("WT","not1d","not4d")
barplot(t(perc),col=c("firebrick3","chartreuse3","dodgerblue"),ylab="% depth in frame",main=cc[pp])

}

plot(1,1,col="white",xaxt='n',yaxt="n",ylab="",xlab="",frame=FALSE)
legend("left",legend=c("Soluble","Frame 0","Frame 1","Frame 2"),col=c("white","firebrick3","chartreuse3","dodgerblue"),lwd=10,bty="n",cex=1.3)

dev.off()



pickLis=list(solAss[,1],solEnc[,1],c(solAss[2,1],solAss[2,1]))
cc=c("NPC_Associated","NPC_Encoded","NUP1")


for(ii in 1:3){


picks=pickLis[[ii]]


scal=list()
k=1
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
c1=lapply(kp[[k]],cdf)





#8.5 mil CH, 5883
#10 mil, 7893
#100*57722/66650000


rpkm=read.table("~/Desktop/Not5_5P/rpkmTab_All_5PSeq_Samples.xls",sep="\t",header=TRUE)

#match(rownames(rpkm),names(cds))


wtAI=log(1+rowMeans(rpkm[,grep("FRG_WT_AI_TOT",colnames(rpkm))]),2)
wtAI=wtAI[which(wtAI>0)]

len=unlist(lapply(cl[[1]],length))-200
names(len)=names(cl[[1]])
len=log(1+len,2)
len=len[match(names(wtAI),names(len))]


wtAI=wtAI[order(wtAI)]
len=len[order(len)]

expL=list(names(wtAI[1:floor(length(wtAI)/3)]),names(wtAI[(ceiling(length(wtAI)/3)):(floor(2*length(wtAI)/3))]),names(wtAI[(ceiling(2*length(wtAI)/3)):length(wtAI)]))

lenL=list(names(len[1:floor(length(len)/3)]),names(len[(ceiling(length(len)/3)):(floor(2*length(len)/3))]),names(len[(ceiling(2*length(len)/3)):length(len)]))



solEnc=read.xlsx("~/Desktop/Not5_5P/NPC/Solubility_NPC_Encoded_30Nov_2022.xlsx")
solAss=read.xlsx("~/Desktop/Not5_5P/NPC/Solubility_NPC_Associated_30Nov_2022.xlsx")


pdf(paste0("~/Desktop/Not5_5P/NPC/Metagene_5P_",paste(ct,collapse="_"),"Exp","Len","SCALE_Spike=",spike,"_WideFC_Nov30_2022_NPC_NUP1.pdf"),width=10,height=4)
par(mfrow=c(1,3))


lili=list()

# pickLis=list(names(kpo[[1]]),solEnc[,1],solAss[,1])
# cc=c("All_Genes","NPC_Encoded","NPC_Associated")

pickLis=list(solAss[,1],solEnc[,1],c(solAss[2,1],solAss[2,1]))
cc=c("NPC_Associated","NPC_Encoded","NUP1")


for(ii in 1:3){


picks=pickLis[[ii]]


scal=list()
k=1
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
c1=lapply(kp[[k]],cdf)
cds=unlist(c1)
nums=unlist(lapply(c1,hun))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag

for(k in 2:length(kp)){
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
cds=unlist(lapply(kp[[k]],cdf))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag
}

names(scal)=names(kp)

adf=as.data.frame(scal)
adf=adf[,2*(1:(ncol(adf)/2))]

# af=t(sweep(adf,STATS=colSums(adf)/1000000,MARGIN=2,FUN="/"))
af=t(sweep(adf,STATS=colMeans(adf),MARGIN=2,FUN="/"))


setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_SpikeIn_PUF3Indicated/")
pi=list.files(pattern="*PUF3Indicated.xlsx")

xl=paste0(ct,"_Over_frg_",ct,"PUF3Indicated.xlsx")
mv=match(xl,pi)

meanLisS=list()
for(state in 1:length(ct)){
meanLisS[[state]]=colMeans(af[grep(paste0("^",ct[state]),rownames(af)),])
rx=2^mean(read.xlsx(pi[mv[state]])[,6])
rn=2^mean(read.xlsx(paste0("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_NoSpikeIn_PUF3Indicated/",gsub("PUF3Indicated.xlsx","NoSpikePUF3Indicated.xlsx",pi[mv[state]])))[,6])
print(state)
print(rx)
print(rn)
if(spike==TRUE){
meanLisS[[state]]=meanLisS[[state]]*(rx/rn)
}
}

#colz=c("black","chartreuse3","firebrick3","dodgerblue","darkorange","purple")
#colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")

colz=colo[match(ct,samp)]


miny=min(unlist(meanLisS))
maxy=max(unlist(meanLisS))


#pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"Exp",ii,"Len",jj,"SCALE_DEC2020.pdf"))

#for(ci in 1:length(ct)){
for(ci in 1:length(ct)){

#ylim=c(miny,maxy)
ylabi=""
if(spike==TRUE){
ylabi=" spike-in weighted to RNASeq"
}

if(ci == 1){
print("TRUE")
print(ci)
plot(1:100,meanLisS[[ci]],ylim=c(miny,maxy),col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P depth",bty="n",main=paste0(cc[ii]))
} else {
print(ci)
points(1:100,meanLisS[[ci]],col=colz[ci],type="l",lwd=2)
}
}

if(ii==2){
legend("topleft",legend=gsub("t_AI_","T ",gsub("p_AI_","d ",ct)),col=colz[1:length(ct)],lwd=5,bty="n",cex=1)
}
abline(v=0,lty=2)
abline(v=100,lty=2)




}

dev.off()

}














#aggregations of not proteins
#recruitment of chaperones to help folding?




library(openxslx)

cl=readRDS(file="~/Desktop/Not5_5P/CDSCounts.rds")
names(cl)=gsub("�","D_",names(cl))



#ct=c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol")
ct=c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot")
#ct=c("WTsol","not4D_sol","not5D_sol")
#ct=c("WTtot","not4D_tot","not5D_tot")
#ct=c("WT_AIsol","not4D_sol","not5D_sol")


#


inz=unlist(unlist(lizo))


#lizo=list(c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot"),c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol"))
#
#
#lizo=list(c("Wt_AI_tot"),c("Wt_AI_sol"),c("WTsol"),c("WTtot"))
#
#lizo=list(c("Wt_AI_tot"))
#
#lizo=list(c("WTsol","not4D_sol","not5D_sol"))

#
#setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_SpikeIn_PUF3Indicated/")
#
#library(openxlsx)
#
#pi=list.files(pattern="*PUF3Indicated.xlsx")
#
#
#genes=read.xlsx(pi[grep(paste0("^",inz[1]),pi)])[,4]
#
#for(pp in 2:length(inz)){
#
#genes=intersect(genes,read.xlsx(pi[grep(paste0("^",inz[1]),pi)])[,4])
#
#nam=pi[pp]
#top=gsub("_Over.+","",nam)
#bot=gsub("PUF3Indicated.xlsx","",gsub(".+_Over_","",nam))
#
#print(paste(top,bot))
#rt=read.xlsx(nam)
#print(mean(rt[,6]))
#
#diff=ag[match(c(top,bot),ag[,1]),][1,2]-ag[match(c(top,bot),ag[,1]),][2,2]
#print(diff)
#
#}


samp=c("WTsol","Wt_AI_sol","WTtot","Wt_AI_tot",
"not5D_sol","Not5p_AI_sol","not5D_tot","Not5p_AI_tot",
"not4D_sol","Not4p_AI_sol","not4D_tot","Not4p_AI_tot",
"Not1p_AI_sol","Not1p_AI_tot")

colo=c("black","black","black","black",
"firebrick3","firebrick3","firebrick3","firebrick3",
"chartreuse3","chartreuse3","chartreuse3","chartreuse3",
"dodgerblue","dodgerblue")

colz=colo[match(ct,samp)]


lizo=list(c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot"),c("WTsol","not4D_sol","not5D_sol"),c("WTtot","not4D_tot","not5D_tot"),c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol"))


spike=TRUE

for(lo in 1:length(lizo)){
#for(lo in 2){

ct=lizo[[lo]]


grabs=c()
for(cac in 1:length(ct)){
grabs=c(grabs,grep(paste0("^",ct[cac]),names(cl)))
}

kp=cl[grabs]
cat("\n\n")
print(names(kp))
kpo=kp


cdf=function(v){
remov=c(1:100,(length(v)-99):length(v))
return(as.numeric(v[-remov]))
}

hun=function(vv){
return(ceiling(100*(1:length(vv))/length(vv)))
}


#8.5 mil CH, 5883
#10 mil, 7893
#100*57722/66650000


rpkm=read.table("~/Desktop/Not5_5P/rpkmTab_All_5PSeq_Samples.xls",sep="\t",header=TRUE)

#match(rownames(rpkm),names(cds))


wtAI=log(1+rowMeans(rpkm[,grep("FRG_WT_AI_TOT",colnames(rpkm))]),2)
wtAI=wtAI[which(wtAI>0)]

len=unlist(lapply(cl[[1]],length))-200
names(len)=names(cl[[1]])
len=log(1+len,2)
len=len[match(names(wtAI),names(len))]


wtAI=wtAI[order(wtAI)]
len=len[order(len)]

expL=list(names(wtAI[1:floor(length(wtAI)/3)]),names(wtAI[(ceiling(length(wtAI)/3)):(floor(2*length(wtAI)/3))]),names(wtAI[(ceiling(2*length(wtAI)/3)):length(wtAI)]))

lenL=list(names(len[1:floor(length(len)/3)]),names(len[(ceiling(length(len)/3)):(floor(2*length(len)/3))]),names(len[(ceiling(2*length(len)/3)):length(len)]))


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"Exp","Len","SCALE_Spike=",spike,"_WideFC_Feb9_2021a.pdf"),width=14)
par(mfrow=c(3,3))

cc=c("Low","Med","High")

lili=list()

for(ii in 1:3){
lili[[ii]]=list(intersect(expL[[ii]],lenL[[1]]),intersect(expL[[ii]],lenL[[2]]),intersect(expL[[ii]],lenL[[3]]))
names(lili[[ii]])=paste(cc[ii],"Expr,",cc[1:3],"Length")
}

for(ii in 1:3){
for(jj in 1:3){

#for(ii in 3){
#for(jj in 3){

picks=intersect(expL[[ii]],lenL[[jj]])


##ag=aggregate(cds[which(names(cds)%in%picks==TRUE)],by=list(nums),FUN=sum)
#
#print(match(picks,names(cds)))
#
#}}


#kp=kpo

scal=list()
k=1
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
c1=lapply(kp[[k]],cdf)
cds=unlist(c1)
nums=unlist(lapply(c1,hun))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag

for(k in 2:length(kp)){
kp[[k]]=kpo[[k]][which(names(kpo[[k]])%in%picks==TRUE)]
cds=unlist(lapply(kp[[k]],cdf))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag
}

names(scal)=names(kp)

adf=as.data.frame(scal)
adf=adf[,2*(1:(ncol(adf)/2))]

af=t(sweep(adf,STATS=colSums(adf)/1000000,MARGIN=2,FUN="/"))


setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_SpikeIn_PUF3Indicated/")
pi=list.files(pattern="*PUF3Indicated.xlsx")

xl=paste0(ct,"_Over_frg_",ct,"PUF3Indicated.xlsx")
mv=match(xl,pi)

meanLisS=list()
for(state in 1:length(ct)){
meanLisS[[state]]=colMeans(af[grep(paste0("^",ct[state]),rownames(af)),])
rx=2^mean(read.xlsx(pi[mv[state]])[,6])
rn=2^mean(read.xlsx(paste0("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_NoSpikeIn_PUF3Indicated/",gsub("PUF3Indicated.xlsx","NoSpikePUF3Indicated.xlsx",pi[mv[state]])))[,6])
print(state)
print(rx)
print(rn)
if(spike==TRUE){
meanLisS[[state]]=meanLisS[[state]]*(rx/rn)
}
}

#colz=c("black","chartreuse3","firebrick3","dodgerblue","darkorange","purple")
#colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")

colz=colo[match(ct,samp)]


miny=min(unlist(meanLisS))
maxy=max(unlist(meanLisS))


#pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"Exp",ii,"Len",jj,"SCALE_DEC2020.pdf"))

#for(ci in 1:length(ct)){
for(ci in 1:length(ct)){

#ylim=c(miny,maxy)
ylabi=""
if(spike==TRUE){
ylabi=" spike-in weighted to RNASeq"
}

if(ci == 1){
print("TRUE")
print(ci)
plot(1:100,meanLisS[[ci]],ylim=c(miny,maxy),col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P depth",bty="n",main=paste0(cc[ii]," Expr, ",cc[jj]," Length,",ylabi))
} else {
print(ci)
points(1:100,meanLisS[[ci]],col=colz[ci],type="l",lwd=2)
}
}

if(ii*jj==9){
legend("topleft",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1)
abline(v=0,lty=2)
abline(v=100,lty=2)
}


#if(ci == 1){
#print("TRUE")
#print(ci)
#plot(1:100,meanLisS[[ci]],ylim=c(miny,maxy),col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P Depth",bty="n",main=paste(cc[ii],"Expr,",cc[jj],"Length"))
#} else {
#print(ci)
##points(1:100,meanLisS[[ci]],col=colz[ci],type="l",lwd=2)
#}
#}
#
#if(ii*jj==9){
#legend("topleft",legend=ct[1],col=colz[1],lwd=5,bty="n",cex=1)
#abline(v=0,lty=2)
#abline(v=100,lty=2)
#}


}}

dev.off()

}





setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_SpikeIn_PUF3Indicated/")

lf=list.files(pattern="*.totPUF3Indicated.xlsx")
lf=grep("^frg_",lf,value=TRUE)
lf=grep("sol_Over",lf,value=TRUE)

lf=lf[c(6,1,3,5,7,2,4)]

lili=list()

for(ii in 1:3){
lili[[ii]]=list(intersect(expL[[ii]],lenL[[1]]),intersect(expL[[ii]],lenL[[2]]),intersect(expL[[ii]],lenL[[3]]))
names(lili[[ii]])=paste(cc[ii],"Expr,",cc[1:3],"Length")
}



pdf(paste0("~/Desktop/Not5_5P/Boxplot_SolubilitySplitByExpressionAndLength_Feb2021.pdf"),width=7)

par(mfrow=c(3,3))
for(ii in 1:3){
for(jj in 1:3){

lf1=gsub("sol.+","",gsub("frg_","",gsub("_sol.+","",gsub("PUF3Indicated.xlsx","",lf))))
boxl=list()
for( ax in 1:length(lf) ){
rt=read.xlsx(lf[ax])
boxl[[ax]]=rt[na.omit(match(lili[[ii]][[jj]],rt[,4])),6]
}
names(boxl)=lf1

boxplot(boxl,main=paste(cc[ii],"Expr,",cc[jj],"Length"),pch=16,cex=0.5,ylab="Solubility",las=2,ylim=c(-4,5))

}}

dev.off()





setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_NoSpikeIn_PUF3Indicated/")

lf=list.files(pattern="*.totNoSpikePUF3Indicated.xlsx")
lf=grep("^frg_",lf,value=TRUE)
lf=grep("sol_Over",lf,value=TRUE)

lf=lf[c(6,1,3,5,7,2,4)]



pdf(paste0("~/Desktop/Not5_5P/Boxplot_SolubilitySplitByExpressionAndLength_Feb2021_NOSpikeN4MinN1_LIM.pdf"),width=7)

par(mfrow=c(3,3))
for(ii in 1:3){
for(jj in 1:3){

lf1=gsub("sol.+","",gsub("frg_","",gsub("_sol.+","",gsub("PUF3Indicated.xlsx","",lf))))
boxl=list()
for( ax in 1:length(lf) ){
rt=read.xlsx(lf[ax])
boxl[[ax]]=rt[na.omit(match(lili[[ii]][[jj]],rt[,4])),6]
names(boxl[[ax]])=rt[na.omit(match(lili[[ii]][[jj]],rt[,4])),4]
}
names(boxl)=lf1

intz=intersect(names(boxl[[2]]),names(boxl[[3]]))
N4mN1=boxl[[3]][match(intz,names(boxl[[3]]))]-boxl[[2]][match(intz,names(boxl[[2]]))]

boxl[[8]]=N4mN1
names(boxl)[8]="Not4-Not1_AI"

boxplot(boxl,main=paste(cc[ii],"Expr,",cc[jj],"Length"),pch=16,cex=0.5,ylab="Solubility No Spike",las=2,ylim=c(-3,3))
abline(h=0,col="red")
}}

dev.off()




######


setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_NoSpikeIn_PUF3Indicated/")

lf=list.files(pattern="*.totNoSpikePUF3Indicated.xlsx")
lf=grep("^frg_",lf,value=TRUE)
lf=grep("sol_Over",lf,value=TRUE)

lf=lf[c(6,1,3,5,7,2,4)]


pdf(paste0("~/Desktop/Not5_5P/Boxplot_SolubilityAllGenes_Feb2021_NOSpikeN4MinN1_LIM.pdf"),width=7)

lf1=gsub("sol.+","",gsub("frg_","",gsub("_sol.+","",gsub("PUF3Indicated.xlsx","",lf))))
boxl=list()
for( ax in 1:length(lf) ){
rt=read.xlsx(lf[ax])
boxl[[ax]]=rt[,6]
names(boxl[[ax]])=rt[,4]
}
names(boxl)=lf1

intz=intersect(names(boxl[[2]]),names(boxl[[3]]))
N4mN1=boxl[[3]][match(intz,names(boxl[[3]]))]-boxl[[2]][match(intz,names(boxl[[2]]))]

boxl[[8]]=N4mN1
names(boxl)[8]="Not4-Not1_AI"

boxplot(boxl,pch=16,cex=0.5,ylab="Solubility No Spike",las=2,ylim=c(-3,3))
abline(h=0,col="red")


dev.off()



#######




setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_SpikeIn_PUF3Indicated/")

lf=list.files(pattern="*.totPUF3Indicated.xlsx")
lf=grep("^frg_",lf,value=TRUE)
lf=grep("sol_Over",lf,value=TRUE)

lf=lf[c(6,1,3,5,7,2,4)]




pdf(paste0("~/Desktop/Not5_5P/Boxplot_SolubilityAllGenes_Feb2021_WITH_SpikeN4MinN1_LIM.pdf"),width=7)

lf1=gsub("sol.+","",gsub("frg_","",gsub("_sol.+","",gsub("PUF3Indicated.xlsx","",lf))))
boxl=list()
for( ax in 1:length(lf) ){
rt=read.xlsx(lf[ax])
boxl[[ax]]=rt[,6]
names(boxl[[ax]])=rt[,4]
}
names(boxl)=lf1

intz=intersect(names(boxl[[2]]),names(boxl[[3]]))
N4mN1=boxl[[3]][match(intz,names(boxl[[3]]))]-boxl[[2]][match(intz,names(boxl[[2]]))]

#boxl[[8]]=N4mN1
#names(boxl)[8]="Not4-Not1_AI"

boxplot(boxl,pch=16,cex=0.5,ylab="Solubility With Spike",las=2)
abline(h=0,col="red")


dev.off()





########



setwd("~/Desktop/Not5_5P/DESeq2/xlsx_SpikeIn_AndNotSpike_PUF3_Indicated/xlsx_SpikeIn_PUF3Indicated/")

lf=list.files(pattern="*.totPUF3Indicated.xlsx")
lf=grep("^frg_",lf,value=TRUE,invert=TRUE)
lf=grep("_Over_frg",lf,value=TRUE)
lfTot=lf[2:8]

lfTot=lfTot[c(6,1,3,5,7,2,4)]


lf=list.files(pattern="*.solPUF3Indicated.xlsx")
lf=grep("^frg_",lf,value=TRUE,invert=TRUE)
lf=grep("_Over_frg",lf,value=TRUE)
lfSol=lf

lfSol=lfSol[c(6,1,3,5,7,2,4)]

lf=c(lfSol,lfTot)

lf=grep("_AI",lf,value=TRUE)

pdf(paste0("~/Desktop/Not5_5P/Boxplot_AllGenes_Degradation_per_Transcript_WithSpikeN4MinN1_LIM_Nov23_2021.pdf"),width=7)

par(mar=c(7.1, 5.1, 1.1, 1.1))


lf1=gsub("tot.+","Tot",gsub("_tot.+","Tot",gsub("sol.+","Sol",gsub("frg_","",gsub("_sol.+","Sol",gsub("PUF3Indicated.xlsx","",lf))))))


boxl=list()
for( ax in 1:length(lf) ){
rt=read.xlsx(lf[ax])
boxl[[ax]]=rt[,6]
names(boxl[[ax]])=rt[,4]
}
names(boxl)=lf1

intz=intersect(names(boxl[[2]]),names(boxl[[3]]))
N4mN1=boxl[[3]][match(intz,names(boxl[[3]]))]-boxl[[2]][match(intz,names(boxl[[2]]))]

#boxl[[8]]=N4mN1
#names(boxl)[8]="Not4-Not1_AI"

boxplot(boxl,pch=16,cex=0.5,ylab="5'P Degradation per Transcript, With Spike",las=2)
abline(h=0,col="red")


dev.off()



for(ii in 1:4+4){
for(jj in 1:4+4){

if(jj>ii){
print(paste(ii,jj))
print(t.test(boxl[[ii]],boxl[[jj]])$p.value)

}}







####
####PUF3 part


#ct=c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot")
ct=c("WTsol","not4D_sol","not5D_sol")
#ct=c("WTtot","not4D_tot","not5D_tot")
#ct=c("WT_AIsol","not4D_sol","not5D_sol")

lizo=list(c("Wt_AI_tot","Not4p_AI_tot","Not5p_AI_tot","Not1p_AI_tot"),c("WTsol","not4D_sol","not5D_sol"),c("WTtot","not4D_tot","not5D_tot"),c("Wt_AI_sol","Not4p_AI_sol","Not5p_AI_sol","Not1p_AI_sol"))


#for(lo in 1:length(lizo)){
for(lo in 4){

ct=lizo[[lo]]

grabs=c()
for(cac in 1:length(ct)){
grabs=c(grabs,grep(paste0("^",ct[cac]),names(cl)))
}

kp=cl[grabs]

library(openxlsx)
c5=read.xlsx("~/Desktop/Not5/elife-40670-supp1-v2_PUF3TargetsPassmore_Core5.xlsx",2)
c5=c5[union(which(c5[,2]==5),which(c5[,3]==5)),]

fcTab=cbind(c5[,1],c5[,1])


p3=function(v){
return(na.omit(v[match(fcTab[,1],gsub(".+_","",names(v)))]))
}

kp=lapply(kp,p3)



cdf=function(v){
remov=c(1:100,(length(v)-99):length(v))
return(as.numeric(v[-remov]))
}

hun=function(vv){
return(ceiling(100*(1:length(vv))/length(vv)))
}


scal=list()
k=1
c1=lapply(kp[[k]],cdf)
cds=unlist(c1)
nums=unlist(lapply(c1,hun))
ag=aggregate(cds,by=list(nums),FUN=sum)

scal[[k]]=ag

for(k in 2:length(kp)){
cds=unlist(lapply(kp[[k]],cdf))
ag=aggregate(cds,by=list(nums),FUN=sum)
scal[[k]]=ag
}

names(scal)=names(kp)

adf=as.data.frame(scal)
adf=adf[,2*(1:(ncol(adf)/2))]

af=t(sweep(adf,STATS=colSums(adf)/1000000,MARGIN=2,FUN="/"))

meanLisS=list()
for(state in 1:length(ct)){
meanLisS[[state]]=colMeans(af[grep(paste0("^",ct[state]),rownames(af)),])
}


colz=c("black","chartreuse3","firebrick3","dodgerblue","darkorange","purple")

miny=min(unlist(meanLisS))
maxy=max(unlist(meanLisS))
pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"SCALE_DEC2020_PUF3TARGETS.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(1:100,meanLisS[[ci]],ylim=c(miny,maxy),col=colz[ci],type="l",lwd=2,xlab="Scaled CDS Position (100 equal bins)",ylab="Normalised 5P Depth",bty="n",main="PUF3 TARGETS ONLY")
} else {
print(ci)
points(1:100,meanLisS[[ci]],col=colz[ci],type="l",lwd=2)
}
}

legend("bottom",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)
abline(v=0,lty=2)
abline(v=100,lty=2)

dev.off()

}



#########
#########




sto=function(v){
return(as.numeric(v[(length(v)-699):length(v)]))
}

kp2=kp
lens=unlist(lapply(kp2[[1]],length))

tLis=list()
for(i in 1:length(kp)){
tLis[[i]]=t(as.data.frame(lapply(kp2[[i]][which(lens>700)],sto)))
}

names(tLis)=names(kp)



names(tLis)=names(kp)

mnt=lapply(tLis,colSums)
weigh=unlist(cdsTot)/1000000

adf=t(as.data.frame(mnt))
sw=sweep(adf,MARGIN=1,STATS=weigh,FUN="/")

meanLis=list()
for(state in 1:length(ct)){
meanLis[[state]]=colMeans(sw[grep(paste0("^",ct[state]),rownames(sw)),])
}

names(meanLis)=ct

colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")

maxy=max(unlist(lapply(meanLis,max)))


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"STOP700.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(-600:99,meanLis[[ci]],col=colz[ci],ylim=c(0,maxy),type="l",lwd=2,xlab="Distance from STOP",ylab="Normalised 5P Depth",bty="n")
} else {
print(ci)
points(-600:99,meanLis[[ci]],col=colz[ci],type="l",lwd=2)
}
}

legend("topright",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)
abline(v=0,lty=2)

dev.off()


########
########




sta=function(v){
return(as.numeric(v[1:700]))
}

kp2=kp
lens=unlist(lapply(kp2[[1]],length))

tLis=list()
for(i in 1:length(kp)){
tLis[[i]]=t(as.data.frame(lapply(kp2[[i]][which(lens>700)],sta)))
}

names(tLis)=names(kp)



names(tLis)=names(kp)

mnt=lapply(tLis,colSums)
weigh=unlist(cdsTot)/1000000

adf=t(as.data.frame(mnt))
sw=sweep(adf,MARGIN=1,STATS=weigh,FUN="/")

meanLis=list()
for(state in 1:length(ct)){
meanLis[[state]]=colMeans(sw[grep(paste0("^",ct[state]),rownames(sw)),])
}

names(meanLis)=ct

colz=c("black","firebrick3","chartreuse3","dodgerblue","darkorange","purple")

maxy=max(unlist(lapply(meanLis,max)))


pdf(paste0("~/Desktop/Not5_5P/Metagene_5P_",paste(ct,collapse="_"),"START700.pdf"))

for(ci in 1:length(ct)){

if(ci == 1){
print("TRUE")
print(ci)
plot(-100:599,meanLis[[ci]],col=colz[ci],ylim=c(0,maxy),type="l",lwd=2,xlab="Distance from START",ylab="Normalised 5P Depth",bty="n")
} else {
print(ci)
points(-100:599,meanLis[[ci]],col=colz[ci],type="l",lwd=2)
}
}

legend("topright",legend=ct,col=colz[1:length(ct)],lwd=5,bty="n",cex=1.2)
abline(v=0,lty=2)

dev.off()



}



DLD1
247,248 amino acid


741+utr

for(ii in 1:length(kp)){
print(names(kp)[ii])
print(kp[[ii]][grep("DLD1",names(kp[[ii]]))])
}


######




fl=readLines("~/Downloads/fivepseq_results02/fivepseq_not5�sol/fivepseq_counts/not5�sol1_04_0318_CTACAC/protein_coding/counts_FULL_LENGTH.txt")

st=strsplit(fl,split="\t")

gtf <- rtracklayer::import("~/Desktop/Not5_5P/Saccharomyces_cerevisiae.R64-1-1.94.gff3")

mr=gsub("_mRNA","",gsub("transcript:","",gtf[which(gtf$type=="mRNA")]$ID))
mrgen=gtf[match(paste0("gene:",mr),gtf$ID)]
mrgen=mrgen[-match(c("YJL047C-A","YLR161W","YEL033W","YLR156W","YLR159W","YLR157W-E"),as.character(mrgen$gene_id))]

tra=read.table("~/Downloads/fivepseq_results02/fivepseq_WTsol/fivepseq_counts/WTsol1_04_0318_CGTGCA/protein_coding/transcript_assembly.txt",header=TRUE)

ln=tra[,4]-tra[,3]
lan=unlist(lapply(st,length))






stlis=lapply(st,sta)
stTab=t(as.data.frame(stlis))







setwd("~/Downloads/fivepseq_results02/")
lf=list.files()
lf=grep("fivepseq_log",lf,invert=TRUE,value=TRUE)
#lf=grep("not5|WT",lf,value=TRUE)[c(3,4,1,2)]

ding=c()
meanL=list()
meanCod=c()
for(l in 1:length(lf)){

setwd(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts"))
lf2=list.files()


for(l2 in 1:length(lf2)){

ding=c(ding,paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts/",lf2[1],"/protein_coding/codon_pauses.txt"))
}
}

cfl=read.table(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts/",lf2[l2],"/protein_coding/codon_pauses.txt"))
cod=read.table(paste0("~/Downloads/fivepseq_results02/",lf[l],"/","fivepseq_counts/",lf2[l2],"/protein_coding/count_distribution.txt"))
suz=sum(cod[,1]*cod[,2])

if(l2==1){
meanTab=(cfl/(suz/1000000))
}

meanTab=meanTab+(cfl/(suz/1000000))

}

meanL[[l]]=meanTab/3

meanCod=cbind(meanCod,(meanTab/3)[,14])

}

names(meanL)=lf
colnames(meanCod)=lf