# USAGE
#Rscript ~/Programming/pyNextGen/R/junctionSeq.R TABLE-SAMPLE-CONDITION.tab QORT-FOLDER/ FLAT-GFF CORR-FILE

# By default, use 20 threads

## Check python good_practices pyensembl for creating CORR-FILE
## 

options(bitmapType="cairo")
library(JunctionSeq)

args = commandArgs(trailingOnly = TRUE)

# Table with columns: samplenames, groups
sampleData = args[1]
## Directory where to find the qort output folder. HAS TO FINISH WITH "/" !!!!!!!!
qortDir = args[2]
## Flattened GTF from qorts
gff.file = args[3]
## Correspondance file (ENS names / associated gene names)
corr.file = args[4]


decoder = read.table(sampleData, col.names=c('ID', 'GROUP'), stringsAsFactors=FALSE, comment="#")

print(decoder)

## Raw from sample group table
countFiles = paste0(qortDir, decoder$ID, '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz')
## OR with name modification
countFiles = paste0(qortDir, gsub('.bam', '',decoder$ID), '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz')


jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$ID,
                               condition=factor(decoder$GROUP),
                               flat.gff.file = gff.file,
                               nCores = 20,
                               analysis.type = "junctionsAndExons",
                               gene.names = corr.file
                               )

writeCompleteResults(jscs,outfile.prefix="./test",save.jscs = TRUE)

buildAllPlots(jscs=jscs,
              outfile.prefix = "./plots/",
              use.plotting.device = "png",
              FDR.threshold = 0.01)


