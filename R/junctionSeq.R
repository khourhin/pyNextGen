################################################################################
                                        # STUB
################################################################################

## Decoder is just a table with columns: samplenames, groups
sampleData = '/home/ekornobis/analysis/batsche/badeg/1.1/junctionSeq/sample_groups.tab'
## Directory where to find the qort output folder. HAS TO FINISH WITH "/" !!!!!!!!
qortDir = '/home/ekornobis/analysis/batsche/badeg/1.1/qorts/'



decoder = read.table(sampleData, col.names=c('ID', 'GROUP'), stringsAsFactors=FALSE)

countFiles = paste0(qortDir, decoder$ID, '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz')

## From qorts
gff.file = '/home/ekornobis/analysis/batsche/badeg/1.1/qorts/Homo_sapiens.GRCh37.75_FLAT.gtf'

library(JunctionSeq)
jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$ID,
                               condition=factor(decoder$GROUP),
                               flat.gff.file = gff.file,
                               nCores = 20,
                               analysis.type = "junctionsAndExons"
                               )

writeCompleteResults(jscs,outfile.prefix="./test",save.jscs = TRUE)

buildAllPlots(jscs=jscs,
              outfile.prefix = "./plots/",
              use.plotting.device = "png",
              FDR.threshold = 0.01)

