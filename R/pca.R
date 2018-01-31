library(edgeR)
library(ggplot2)
library(ggrepel)
library(repr)
library(biomaRt)


get.normalized.counts = function(counts, meta){
    ## Import data and normalization

    ## Create the EdgeR count object
    x = DGEList(counts, group=meta[,1])

    x <- calcNormFactors(x, method="TMM")

    ## As advised in EdgeR manual (MDS plot section)
    all_cpms = cpm(x)
    
    return(all_cpms)
}

get.loadings = function(counts, PC='PC1', n=10){
    # Plot loadings of each genes 
    
    pca = prcomp(t(counts))
    loadings = abs(pca$rotation)
    
    # Calculate the % represented by each loadings for each gene 
    # within a principal component (here by columns)
    loadings.perc = sweep(loadings, 2, colSums(loadings), "/")
    
    plot(sort(loadings.perc[,PC]))
    
    sup.loadings = sort(loadings.perc[,PC], decreasing=T)[0:n]
    res = cbind(annot[names(sup.loadings),], sup.loadings)
    
    #print(annot[names(sort(loadings.perc[,PC], decreasing=T)[0:10]),])
    
    return(res)
}

do.pca = function(counts, meta, file, gene_sel=NULL, lib_sel=NULL,
                  repel_annot=TRUE, legend=TRUE, annot_set=NULL, annot_cor=NULL, 
                  annot_cor_col=2, logTrans=F, col_group_choice=1,
                  shape_group_choice=1){
    
    # counts: the count matrix 
    # meta: a table with sample names as first columns and extra columns are groups 
    # file: The name of the outfile
    # gene_sel: A list of genes to select 
    # lib_sel: A list of libraries to select
    # repel_annot: boolean, if using or not repel_annot
    # legend: boolean, printing legend or not
    # annot_set: TO add further annotations to the plots
    # annot_cor: TO ADD
    # annot_cor_col: TO ADD
    # logTrans: boolean, if performing log2 transformation of the counts
    # col_group_choice: int, number of the column in "meta" to use as groups for colors on plot
    # shape_group_choice: int, number of the column in "meta" to use as groups for shapes on plot
    
    if (!is.null(gene_sel)){
        counts = counts[rownames(counts) %in% gene_sel,]
    }
    
    if (!is.null(lib_sel)){
        counts = counts[,colnames(counts) %in% lib_sel]
    }
    
    if (logTrans == TRUE){
        counts = log(counts + 1)
    }
    
    pca = prcomp(t(counts))
    print(summary(pca))
    groups = meta[rownames(pca$x),, drop=FALSE] 
    scores=data.frame(groups, pca$x[,1:3])
    
    # THIS SHOULD BE IMPROVED !
    # I think the proper way will be to have groups defined as columns in 'scores' rather than outside
    
    p = ggplot(scores, aes(PC1,PC2, guide=FALSE)) +
    geom_point(aes(color=factor(groups[,col_group_choice]),
                   shape=factor(groups[,shape_group_choice])),
               show.legend=legend) +
        labs(x='PC1', y='PC2') + 
        theme_classic() +
        theme(axis.text = element_text(size = 12))

    
    
    if (repel_annot){
        
      p = p + geom_text_repel(data = scores, aes(PC1,PC2, label = rownames(scores)))  
        
    } else if (!is.null(annot_set) & !is.null(annot_cor)){
        
        annot_scores = scores[rownames(scores) %in% annot_set,]
        rownames(annot_scores) = annot_cor[rownames(annot_scores), annot_cor_col]
        p = p + geom_text_repel(data = annot_scores, aes(PC1,PC2, label = rownames(annot_scores)))
    }    
    ggsave(file, p, device='pdf', width=10, height=10, units='cm')
    return(p)
}
