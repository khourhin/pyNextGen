## from Law et al 2016 in F1000Research
# Also used in Li et.al 2014 Nat.Bioteck
## RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR

## THIS FILE IS NOT A SCRIPT !!!!
## DO NOT LOAD IT TO R

## This is designed to be used using the C-c C-j command in ess to
## load one line at the time in order to be more interactive and
## verify the good functionning of the analysis

## TO IMPROVE:
## Check to add NoiseqBio, RUV for ex

library(edgeR)

files <- list.files("~/analysis/muchardt/os_erv/shrimp/counts/counts_for_edgeR", pattern="*_genes*",  full.name=T)

## Can import the data directly out of star/htseq (columns given are geneid and counts)
x <- readDGE(files, columns=c(1,2))

## Changing for nicer names
colnames(x) <- basename(colnames(x))

## 4 control and 4 treated with dioxin
group <- as.factor(c(rep("untreated", 3), rep("treated_30m", 3), rep("treated_2h", 3)))
x$samples$group <- group

## Can add also the lane to see the sequencing effect (see article)

## Normalizing counts
lcpm <- cpm(x, log=TRUE)

## Check how much genes are never expressed in any sample
table(rowSums(x$counts)==0)

## Keep only genes which are expressed in at least 1 group (i.e 4
## samples in our case), i.e with log-CPM > 0
keep.exprs <- rowSums(lcpm>1)>=3
x <- x[keep.exprs,keep.lib.sizes=FALSE]

## Genes left now:
dim(x)

## Trimmed mean of M-values normalization (TMM)
x <- calcNormFactors(x, method="TMM")
x$samples$norm.factors

# Creating design of the differential expression analysis
design <- model.matrix(~0+group) #(~0+group+othergroup+...)
colnames(design) <- gsub("group", "", colnames(design))
design

## Creating the contrasts for the pairwise comparisons
contr.matrix <- makeContrasts(
    untreatedVS30min = untreated-treated_30m,
    untreatedVS2hours = untreated-treated_2h, 
    levels = colnames(design)
)
contr.matrix

# Voom treatment (check article...)
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

## Summarize DEG
summary(decideTests(efit))

# Adding a Log Fold Change threshold of 1 
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

## Possible to do if more than 2 groups
##de.common <- which(dt[,1]!=0 & dt[,2]!=0)
##length(de.common)
##head(tfit$genes$SYMBOL[de.common], n=20)
##vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

## Output all results
write.fit(tfit, dt, file="results.txt")

## Examine individual DE genes (check sorted by pvalues)
tops <- topTreat(tfit, coef=1, n=Inf) ## coef=2 for next pairwise comparison ...
head(tops)

######################################################################
                                        # PLOTS
######################################################################

## To save them:
## at the beginning
## pdf()
## at the end
## dev.off()

samplenames <- colnames(x)
lcpm <- cpm(x, log=TRUE)

## Plotting expression (before and after filtering)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(2,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0, 0.60), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

## Plotting expression (before and after normalization)
boxplot(lcpm, col=col, main="", xaxt="n",)
axis(1, seq(1,8), labels=F)
text(seq(1,8), par("usr")[3]-1.5, labels=samplenames, srt=45, adj=c(1,1,1,1), xpd=T)
title(main="A. Filtered data", ylab="Log-cpm")

## Exploratory plot of differential expression
## Change of color (should have minium 3 groups for that)
#col.group <- group
#levels(col.groups) <- brewer.pal(nlevels(col.group), "Set1")
#col.group <- as.character(col.group)

plotMDS(lcpm, labels=group, col=as.numeric(group))
title(main="A. Sample groups")

## In case second grouping (for example by sequencing lanes)
##col.lane <- lane
##levels(col.lane) <- brewer.pal(nlevels(col.lane), "Set2")
##col.lane <- as.character(col.lane)
##plotMDS(lcpm, labels=lane, col=col.lane)
##title(main="B. Sequencing lanes")

# After bayesian correction
plotSA(efit)

## Making a heatmap of the 100 highest DEG
##library(gplots)
tops.topgenes <- tops$ENSEMBL[1:80]
i <- which(v$genes$ENSEMBL %in% tops.topgenes)

library(pheatmap)
d <- v$E[i,]
# That's the only option I found to not have rognated xlabels (labels_col parameter in pheatmap rognate them)
colnames(d) <- group
pheatmap(d, scale='row', labels_row=v$genes$SYMBOL[i], cex=0.8)

## Graphical representation of DEG
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))

######################################################################
                                        # ADVANCED FEATURES
######################################################################

# To add if necessary/usefull:
# For easier annotations
library(Homo.sapiens)
geneid <- rownames(x)

genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
head(genes)

# Fixing duplicates (TO ENSURE same gene order between annotation and counts)
dup <- genes$ENSEMBL[duplicated(genes$ENSEMBL)]
genes[genes$ENSEMBL %in% dup, ][1:10,]
# keep the first occurence of each gene ID
mat <- match(geneid, genes$ENSEMBL)
genes <- genes[mat,]
genes[genes$ENSEMBL %in% dup,][1:5,]

x$genes <- genes


# Cannot install this package ! Check the article for its use.
#library(Glimma)

## Camera method (check article)
## Seems to work with EntrezID, so have to figure out how to make the conversion
## load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p1.rdata"))
## idx <- ids2indices(Hs.c2,id=rownames(v))
## cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1], inter.gene.cor=0.01)
## head(cam.BasalvsLP)
## barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
