library(AnnotationHub)
library(topGO)

### TOFIX
## For now, this will only print the top 50 GOs according to Fisher
## test I would be better to have only a pvalue threshold but
## apparently Gentable cannot have a pvalue threshold

get.go.map = function(species, genes){
    ah = AnnotationHub()
    orgs <- subset(ah, ah$rdataclass == "OrgDb")
    org.obj = query(orgs, species)[[1]]

    GO.table = mapIds(org.obj, keys=genes, keytype="ENSEMBL", 
                  column="GO", multiVals=function(x){paste(x, collapse=", ")})

    write.table(as.data.frame(GO.table), file="GO_map.txt", col.names=FALSE, sep="\t", quote=FALSE)
    geneIDGO = readMappings(file="GO_map.txt")
    return(geneIDGO)
}

do.GO.analysis = function(ontology, assayed.genes, de.genes, geneIDGO, topNodes=topNodes, prefix="OUT_PREFIX_"){
    
    gene.vector=factor(as.integer(assayed.genes%in%de.genes))
    names(gene.vector)=assayed.genes

    table(gene.vector)
    
    GOdata <- new("topGOdata", ontology = ontology, allGenes = gene.vector ,
                  annot = annFUN.gene2GO, gene2GO = geneIDGO)
    
    # Test based on gene counts
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    # Tests based on gene scores (??? Which are those)
    #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    #resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   #classicKS = resultKS, elimKS = resultKS.elim,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes=topNodes)
    
    #allRes = allRes[allRes$classicFisher < pvalCutOff,]
    write.csv(allRes, paste(prefix, 'GOEA.csv', sep="_"))
    
    return(allRes)
}

do.GO.analysis.for.all.files = function(DE_genes_ids_files, ontology, assayed.genes, geneIDGO, topNodes=topNodes){
    
    
    for (gene_list_file in DE_genes_ids_files){
        print(gene_list_file)
        compa_name = basename(gene_list_file)
        de.genes = scan(gene_list_file, what='raw')
        
        res = do.GO.analysis(ontology, assayed.genes, de.genes, geneIDGO, topNodes=topNodes, prefix=compa_name)
    }
        
}

do.GO.analysis.for.all.list = function(DE_genes_list, ontology, assayed.genes, geneIDGO, topNodes=topNodes){
    
    for (i in 1:length(DE_genes_list)){
        compa_name = names(DE_genes_list)[i]
        message(paste0("GO analysis for: ", compa_name))
        de.genes = DE_genes_list[[i]]
        message(paste0("DE genes: ", length(de.genes), " | Assayed genes: ", length(assayed.genes)))
        
        res = do.GO.analysis(ontology, assayed.genes, de.genes, geneIDGO, topNodes=topNodes, prefix=compa_name)
    }
        
}
