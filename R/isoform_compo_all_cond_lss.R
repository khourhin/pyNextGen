library(grid)
library(ggplot2)
library(gridExtra)
library(reshape2)

                                        #df = read.table('/home/ekornobis/analysis/allemand/gphn/2017_01_18/total_clusters_nom/table_isoName_exonCode_freq_sorted.tab',header=TRUE)
df = read.table('/home/ekornobis/analysis/allemand/gphn/2017_01_18/table_isoName_exonCode_freq_sorted.tab',header=TRUE)

s = df[,1:41]
s = melt(s, id='Name')

# To keep the order in the input table
s$Name = as.character(s$Name)
s$Name = factor(s$Name, levels=rev(unique(s$Name)))

colors=c('white','chartreuse4','red','red','red','red','red')
#chartreuse4

p1 = ggplot(s, aes(variable, Name, fill=factor(value))) +
    geom_tile(colour="black", size=0.25) +
    scale_fill_manual(values=colors) +
    guides(fill=FALSE) +
    xlab("Exons") + ylab("Isoforms") +
    scale_y_discrete(breaks = c("GPHN_1","GPHN_50", "GPHN_100", "GPHN_150", "GPHN_200", "GPHN_250", "GPHN_277")) +
    theme(axis.text.y = element_text(size=10), axis.title=element_text(size=20,face="bold"))

# NOT FULLY WORKING (DOES NOT SORT FREQUENCIES)

p2 = ggplot(df, aes(Name, freq)) +
    geom_bar(stat="identity") +
    xlab("") + ylab("Frequencies") + 
    coord_flip() +
    scale_y_log10() +
    scale_y_continuous(trans='log10', breaks = c(2,5,10,25,50,100,200,500,2000,5000,20000)) +
    theme(axis.ticks = element_blank(), axis.text.y=element_blank(), axis.title=element_text(size=20,face="bold"))

#p2 = qplot(iso_freq, geom="bar", xlab = "",) + coord_flip() + theme(axis.ticks = element_blank())

grid.arrange(p1,p2, widths=c(3/4, 1/4))
dev.copy2pdf()
