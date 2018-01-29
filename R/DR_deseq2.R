library(DESeq2)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(BiocParallel)


### SPECIFIC
# For now, only getting human hg19 annotations

## USING 10 threads by default
register(MulticoreParam(10))


check_counts_meta_tables = function(counts, meta){
    # Verify rownames in meta are the same as colnames in counts
    if (!isTRUE(all.equal(rownames(meta),colnames(counts)))){
    warning("Metadata doesn't seem to fit count matrix. Check both inputs")
    }
    else {
        message("OK: Count and meta tables seems to correspond")
    }
}


make_DR = function(counts, meta, design){

    check_counts_meta_tables(counts, meta)
    
    dds = DESeqDataSetFromMatrix(countData=counts, colData=meta, design=design)

    ## Filter out features with 1 or less read mapped in total (no need of filtering so much as in edgeR)
    dds <- dds[ rowSums(counts(dds)) > 1, ]

    ## Set factors levels of reference for later comparison
                                        #dds$group <- relevel(dds$group4, ref=group_ref)

    dds = DESeq(dds, parallel=TRUE)
    return(dds)
}

export_results = function(counts, res, annot){
    
    annot_counts = merge(annot,counts, by=0, all=TRUE)
    write.csv(annot_counts, paste(prefix, 'raw_counts_annot.csv', sep="_"), row.names=FALSE)

    ## Add the annotation produced earlier
    total_res =  merge(annot,as.data.frame(res), by=0, all.y=TRUE)

    ## UCSC specific !!! Might have to remove the "chr" at some point
    total_res$ucsc = paste('http://genome.ucsc.edu/cgi-bin/hgTracks?org=', 'human',
                           '&db=','hg19',
                           '&position=chr', total_res$chromosome_name, ':', total_res$start_position, '-', total_res$end_position,
                           sep='')

    write.csv(total_res, file=paste(prefix, 'degs_annotated_DESeq2.csv', sep="_"), row.names=FALSE)
    return(total_res)
}


pairwise_comparison = function(comps, meta_col_name){
    ## Perform all comparison specified in comps (a list of size 2
    ## vectors). Using the groups specified in metadata table column 'meta_col'
    
    compa_res = list()

    for (comp in comps){
        compa_name = paste(c(comp[1], comp[2]), collapse="_VS_")
        message(compa_name)
        res_i = results(dds, addMLE=TRUE, contrast=c(meta_col_name, comp[1], comp[2]), alpha=0.05, parallel=TRUE)
        summary(res_i)
        
        compa_res[[compa_name]] = res_i
    }
    return(compa_res)
}


export_pairwise = function(res, annot, prefix=''){
    ## For pairwise comparisons
    for (compa in names(res)){
        tmp_res = merge(annot, as.data.frame(res[[compa]]), by=0, all.y=TRUE)
        write.csv(tmp_res, file=paste(prefix, '_',  compa, '.csv', sep=''), row.names=FALSE)
    }
}

volcano_plot = function(res, title=''){
    ## From DESEQ2 result, draw a volcano plot
    with(res, plot(log2FoldChange, -log10(padj), pch=20, main=title))
    with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
    with(subset(res, padj<.05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
}

## SPECIFIC
annotate = function(counts){
    ## From a count matrix with ensembl IDs get the corresponding
    ## annotation from bioMart

    ## Select database and dataset 
    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

    ## Edit to fetch the correct database
    annot <- getBM(attributes = c("ensembl_gene_id", 
                                  "external_gene_name",
                                  "chromosome_name",
                                  "start_position", 
                                  "end_position", 
                                  "strand"), filter="ensembl_gene_id", values=rownames(counts),mart=grch37)

    rownames(annot) = annot$ensembl_gene_id
    annot = annot[,-1]
    return(annot)
}
