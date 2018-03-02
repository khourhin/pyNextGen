library(DESeq2)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(BiocParallel)


### LIMITATIONS
## For now, only getting human hg19 and mm38 annotations

## DEFAULTS
## Normalized counts using counts function in DESEQ (could use rlog and vst as well)

## IMPROVEMENTS
## - ucsc link ?
## For pairwise export: currently annotation is run for each comparison. Could be launched only once


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


export_counts = function(dds, prefix='', species=''){
    #rlog = rlog(dds, blind=FALSE)
    #vsd = varianceStabilizingTransformation(dds, blind=FALSE)

    counts = counts(dds)
    norm_counts = counts(dds, normalized=TRUE)

    if (species != ''){
        annot = annotate(counts, species=species)

        counts = merge(annot, counts, by=0, all.y=TRUE)
        norm_counts = merge(annot, norm_counts, by=0, all.y=TRUE)
    }
    
    write.csv(counts, paste0(prefix, 'counts_raw.csv'), row.names=FALSE)
    write.csv(norm_counts, paste0(prefix, 'counts_norm.csv'), row.names=FALSE)
}


export_results = function(res, species='', prefix=''){
    ## Export the DE results table with annotation if species provided

    if (species != ''){
        annot = annotate(as.data.frame(res), species=species)
        res =  merge(annot,as.data.frame(res), by=0, all.y=TRUE)
    }
    
    write.csv(res, file=paste(prefix, 'degs_DESeq2.csv', sep="_"))
    return(res)
}


pairwise_comparison = function(comps, meta_col_name){
    ## Perform all comparison specified in comps (a list of size 2
    ## vectors). Using the groups specified in metadata table column 'meta_col'
    
    compa_res = list()

    for (comp in comps){
        compa_name = paste(c(comp[1], comp[2]), collapse="_VS_")
        message(compa_name)
        res_i = results(dds, contrast=c(meta_col_name, comp[1], comp[2]), alpha=0.05, parallel=FALSE)
        summary(res_i)
        
        compa_res[[compa_name]] = res_i
    }
    return(compa_res)
}


export_pairwise = function(res, species='', prefix=''){
    ## For pairwise comparisons
    for (compa in names(res)){
        export_results(res[[compa]], species=species, prefix=paste(prefix, '_', compa, sep=''))
        #tmp_res = merge(annot, as.data.frame(res[[compa]]), by=0, all.y=TRUE)
        #write.csv(tmp_res, file=paste(prefix, '_',  compa, '.csv', sep=''), row.names=FALSE)
    }
}

volcano_plot = function(res, title=''){
    ## From DESEQ2 result, draw a volcano plot
    with(res, plot(log2FoldChange, -log10(padj), pch=20, main=title))
    with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
    with(subset(res, padj<.05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
}


annotate = function(df, species=''){
    ## From a dataframe with ensembl identifiers as rownames get the
    ## corresponding annotation from bioMart

    ## Select database and dataset
    if (species == 'human'){
        ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    }
    else if (species == 'mouse'){
        ensembl = useMart(biomart="ensembl",
                          dataset="mmusculus_gene_ensembl")
    }
    else stop('Unavailable or No species selected for annotation.')

    ## Edit to fetch the correct database
    annot <- getBM(attributes = c("ensembl_gene_id", 
                                  "external_gene_name",
                                  "chromosome_name",
                                  "start_position", 
                                  "end_position", 
                                  "strand"), filter="ensembl_gene_id", values=rownames(counts),mart=ensembl)

    rownames(annot) = annot$ensembl_gene_id
    annot = annot[,-1]
    return(annot)
}
