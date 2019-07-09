
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.


# Edit these to match your data file names: 
input="GOpatchHost.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="plob_iso2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
	largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest=5,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results
quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
#	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	absValue=-log(0.05,10), 
	level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
	level2=0.01, # FDR cutoff to print in regular (not italic) font.
	level3=0.005, # FDR cutoff to print in large bold font.
	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5, # height of the hierarchical clustering tree
	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed "dodgerblue2","firebrick1","skyblue2","lightcoral"
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[order(results$pval),]

################################
####### Heatmaps of 'good genes'
################################

#Note: can run this script anytime after desired GO MWU has been run

library(pheatmap)
library(RColorBrewer)

gg=read.table("plob_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)  #the iso2gene file for your species 
dfull=read.csv("hostVSDandPVALS_no_g4_deseq1_4jun_plusPC1.csv") 
rownames(dfull)<-dfull$X #make gene names rownames
edata=c(2:15) #only columns with your expression data 
d=dfull[dfull$pval.p<=0.05 & !is.na(dfull$pval.p),]

golist=read.table("CC_GOpatchHost.csv",sep="	",header=T) #read in proper GO list from gomwu output - BP/MF/CC
d=read.csv("VSDs_GObinaryHostPatch.csv") #read in VSD file for sig DE genes



in.mwu=paste("MWU",goDivision,input,sep="_")
pv=read.table(in.mwu,header=T)

#######################

gene=subset(pv,name=="mitochondrial protein complex") #write GO term of interest here from your sig list#
t=gene[,5]
t

is=golist$seq[grep(t,golist$term,ignore.case=T)]
length(is) #should match denominator in MWU dendrogram plot for term

#####################loop through genes matching GO term in module or in "good expression" subset

sel=c();gnms=c()
for ( i in is){
	if (i %in% d$X){
		sel=rbind(sel,d[d$X==i,])
		gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
	}
}
row.names(sel)=paste(gnms,sel$X,sep=".")
nrow(sel)  
rownames(sel)


exp=sel[,edata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp)
#should match 'good genes' numerator in figure
rownames(exp)
#rownames(daat)

means=apply(exp,1,mean) # means of rows
expc=exp-means #rescale expression data so it's up and down relative to mean

#reorder columns to match healthy/wps
expc<-expc[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)] #host


col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.55)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

#quartz()
#pdf("HeatmapPatchUp.pdf",width=8,height=17)
pheatmap(expc,color=col,cluster_cols=F,clustering_distance_rows="correlation") #plot the heatmap
#dev.off()
