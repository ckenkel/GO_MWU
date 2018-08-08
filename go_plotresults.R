#by Carly Kenkel, ckenkel@usc.edu
# last modified: 8 August 2018

#These scripts pull from the completed GO_MWU analysis 

# Make sure your input files are in your working directory
# then edit these commands as needed to match your filenames

gg=read.table("plob_iso2gene.tab",header=F,sep="	")  #the iso2gene file for your species
dfull=read.csv("VSDandPVALSnobadsam_deseq1_6aug.csv") #read in VSD file for entire expression dataset
rownames(dfull)<-dfull$X #make gene names your rownames
efulldata=c(2:15) #the columns where your expression data are
d=read.csv("VSDs_GObinaryPatch.csv") #read in VSD file for categorical subset of interest - e.g. genes in WGCNA module, or those with sig expr diffs; NOTE ignore if using GO_MWU
rownames(d)<-d$X #make gene names your rownames
edata=c(2:15) #the columns where your expression data are in the subset

source("summarySE.R")


fileglob="GOpatch.csv" # put the name of your own GO input file here 
filego="plob_iso2go.tab" #put the name of your iso2go file here
division="BP" #and your desired GO category
cutoff=0.05 # p-value cutoff for the GO categories to plot as a tree; the worst ones will be printed in grey italic font.
level2=0.01 # # p-value cutoff for printing in black regular font.
level3=0.005 # p-value cutoff for printing in large black bold font.
adjusted=T # Use F to plot raw p-values; replace with T to plot adjusted p-values.
txtsize=1  # decrease this one to squeeze more GO descriptions on the same panel.

################################################################

in.mwu=paste("MWU",division,fileglob,sep="_")
in.dissim=paste("dissim",division,filego,sep="_")
in.term=paste(division,fileglob,sep="_")

library(ape)
cutoff=-log(cutoff,10)
pv=read.table(in.mwu,header=T)

#plot(sort(log(pv$pval,10)),main="log10(pvalues)")
#abline(h=(-cutoff),col="red")
row.names(pv)=pv$term
if (adjusted==TRUE) { pvals=pv$p.adj } else { pvals=pv$pval }
heat=data.frame(cbind("pval"=pvals)) 
row.names(heat)=pv$term
heat$pval=-log(heat$pval,10) 
goods=subset(heat,pval>=cutoff)
head(goods,100)
goods.names=row.names(goods)
goods.names

# reading and subsetting dissimilarity matrix
diss=read.table(in.dissim,sep="\t",header=T,check.names=F)
row.names(diss)=names(diss)
diss.goods=diss[goods.names,goods.names]

# renaming rows in the p-value list for plotting
row.names(goods)=paste(pv[pv$term %in% goods.names,]$nseqs,pv[pv$term %in% goods.names,]$name)
row.names(heat)=paste(pv$nseqs,pv$name)
row.names(diss.goods)=paste(pv[pv$term %in% goods.names,]$nseqs,pv[pv$term %in% goods.names,]$name)

row.names(diss.goods) # see if some descriptions are very long, if so, must truncate them for plotting



# clustering terms better than cutoff
GO.categories=as.dist(diss.goods)
cl.goods=hclust(GO.categories,method="average")
labs=cl.goods$labels[cl.goods$order] # saving the labels to order the barplot
goods=goods[labs,]


####Recording significances of top GO terms - rename file below if you want to write out
sig=data.frame(cbind(labs,goods))
sig$goods<-as.numeric(sig$goods)
sig$pval<-(10^(-(sig$goods)))
sig
write.csv(sig,file="Terms_BP_patch_GOMWU.csv",quote=F,row.names=F)
#************

# now, reprint GO term list to make sure it matches original fig
#highlight and run next block of code (lines 80-100) and execute all at once

step=100
left=1
top=step*(2+length(labs))
quartz()
plot(c(1:top)~c(1:top),type="n",axes=F,xlab="",ylab="")
ii=1
for (i in length(labs):1) {
	ypos=top-step*ii
	ii=ii+1
	if (goods[i]> -log(level3,10)) { 
		text(left,ypos,labs[i],font=2,cex=1*txtsize,adj=c(0,0)) 
	} else {
		if (goods[i]>-log(level2,10)) { 
			text(left,ypos,labs[i],font=1,cex=0.8* txtsize,adj=c(0,0))
		} else {
			if (goods[i]>1) { 
				text(left,ypos,labs[i],font=3,cex=0.8* txtsize,col="grey50",adj=c(0,0))
			} else { text(left,ypos,labs[i],font=3,cex=0.6* txtsize,col="grey75",adj=c(0,0)) }
		}
	}
}


###############plotting what significant GO genes are

gene=subset(pv,name=="chromatin organization") #change GO term of interest here
t=rownames(gene)
t

golist=read.table(in.term,header=T)


is=golist$seq[grep(t,golist$term,ignore.case=T)]
length(is) #total number of genes with this GO annotation
#should be the same length as the denominator in your GO plot


#####################first loop through genes matching GO term in expression subset
sel=c();gnms=c()
for ( i in is){
	if (i %in% d$X){
		sel=rbind(sel,d[d$X==i,])
		gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
	}
}
row.names(sel)=paste(gnms,sel$X,sep=".")

exp=sel[,edata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }

#####################then through genes matching GO term in whole dataset
sel2=c();gnms2=c()
for ( i in is){
	if (i %in% dfull$X){
		sel2=rbind(sel2,dfull[dfull$X==i,])
		gnms2=append(gnms2,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
	}
}
row.names(sel2)=paste(gnms2,sel2$X,sep=".")

exp2=sel2[,efulldata]
if (length(exp2[,1])==1) { exp2=rbind(exp2,exp2) }

########return number of genes matching GO terms in each
a=nrow(exp) #genes with GO annotation in specified data subset
b=nrow(exp2) #genes with GO annotation in entire dataset
c=(nrow(d)-a) #genes without GO annotation in specified data subset
f=(nrow(dfull)-b) #genes without GO annotation in entire dataset

###SANITY check!  run or re-do the fisher test just for fun
		sig=c(a,c) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(b,f) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("go","notgo"),sig=c("go","notgo")))
		ff=fisher.test(mm,alternative="greater")
ff #should be a significant value for categorical GO runs...
a #number of genes in the specified subset that are annotated with the GO term

#------- heatmap of specific GOgenes stacked by factor or continuous trait
library(pheatmap)

head(sel)
names(sel)
sorted=sel[,c(3,5,7,9,11,13,15,2,4,6,8,10,12,14)]
names(sorted)

#uncomment and modify following lines instead of those above for a categorical reorder
#if you have continuous trait data of interest

# traits=read.csv("WGCNA_traits.csv")
#rownames(traits)<-traits$Sample #make sample names the rownames
#traits=traits[order(traits$Sym_cell.cm2),] #change trait here to reorder by different factor
#sorted=exp[,rownames(traits)] #sort rows by values of a trait *** if you don't want to sort, just say sorted=exp

##########################################
means=apply(sorted,1,mean) # means of rows

expc=sorted-means #rescale expression data

library(RColorBrewer)
col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.7)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

quartz()
#pdf("HeatmapGOPatchChromatin.pdf",width=8,height=17)
pheatmap(expc,color=col,cluster_cols=F,clustering_distance_rows="correlation") #plot the heatmap
#dev.off()

####plotting individual gene-trait correlations instead of heatmap

texp=data.frame(t(sorted))
type=c(1:length(rownames(texp)))
type[grep("B",rownames(texp))]="Patch"
type[grep("N",rownames(texp))]="Normal"
texp$type<-type
sexp=stack(texp[,1:4]) #change columns to indicate total number of genes
sexp[,3]<-as.factor(texp$type)
head(sexp)
summary(sexp)
unique(sexp$ind) #unique returns a list of all the gene names in your GOterm for the module of choice


gene=subset(sexp,ind=="AT.rich.interactive.domain.containing.protein.4A.O.isogroup61766_g2") #insert name of interest from unique output
quartz()
plot(gene$values~gene$V3,main="isogroup61766_g2",ylab="log2(expr)",xlab=NULL) #can plot regression of trait value vs gene of interest
#abline(lm(traits$Sym_cell.cm2~gene$values))
#summary(lm(traits$Sym_cell.cm2~gene$values)) 

quartz()
barplot(o$Sym_cell.cm2, main="", cex.main=2,
        ylab="Symbiont Density",xlab="sample")

