################################################################################
                                        # STUB
################################################################################

## Decoder is just a table with columns: samplenames, groups
sampleData = '/home/ekornobis/analysis/batsche/badeg/1.1/junctionSeq/sample_groups.tab'
## Directory where to find the qort output folder. HAS TO FINISH WITH "/" !!!!!!!!
qortDir = '/home/ekornobis/analysis/batsche/badeg/1.1/qorts/'



decoder = read.table(sampleData, col.names=c('ID', 'GROUP'))
countFiles = paste0(qortDir, decoder$ID, '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz')

## From qorts
gff.file = '/home/ekornobis/analysis/batsche/badeg/1.1/qorts/Homo_sapiens.GRCh37.75_FLAT.gtf'

library(JunctionSeq)
jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$sample.ID,
                               condition=factor(decoder$group.ID),
                               flat.gff.file = gff.file,
                               nCores = 1,
                               analysis.type = "junctionsAndExons"
                               );

