library("data.table")
source("analysis/R/gafTools.R")
source("analysis/R/removeDuplicateAnnots.R")
source("analysis/R/getNonredDataset.R")

## Clean the input dataset downloaded from maizeGDB
maizev4CrossGenes <-fread("data/maize/gene_model_xref_v4.txt")
selCols <- c("v4_gene_model","v3_gene_model","W22","PH207","Mo17b")
filtCrossGenes <- maizev4CrossGenes[,selCols,with=F]
colnames(filtCrossGenes) <- c("B73v4","B73v3","W22","PH207","Mo17")
fwrite(filtCrossGenes,"data/maize/maizeCrossRef.tsv",sep = "\t",row.names = F)

maizev4CrossGenes <- fread("data/maize/maizeCrossRef.tsv")
colnames(maizev4CrossGenes)

maizeCrossRef_melt <- melt.data.table(maizev4CrossGenes,id.vars = c("B73v3"),measure.vars = c("B73v4","Mo17","PH207","W22"),variable.name = "inbred",value.name = "geneModel")
maizeCrossRef_melt
splitComma <- function(data){
  geneModels <- unlist(strsplit(data,","))
  return(geneModels)
}

maizeCrossRef_genes <- maizeCrossRef_melt[,list(geneModel=splitComma(geneModel)),by=list(inbred,B73v3)]
maizeCrossRef_genes[,list(count=.N),by=list(inbred)]
maizeCrossRef_melt[geneModel!="",list(count=.N),by=list(inbred)]


curatedMaizeAnnot <- read_gaf("data/curated/maize_v3.gold.gaf")


write_gaf(curatedMaizeAnnot[grep("^GRM|^AC",db_object_symbol)],"data/curated/B73v3.curated.gaf")
cat(colnames(curatedMaizeAnnot),sep = '", "')
maizeCrossRef_melt[geneModel!=""]
annot_merge <- merge.data.table(maizeCrossRef_genes,curatedMaizeAnnot,by.x="B73v3",by.y="db_object_symbol")
annot_merge[,db_object_symbol:=geneModel]
annot_merge[,db_object_id:=geneModel]

annot_merge[,list(count=.N,geneCount=length(unique(db_object_symbol))),by=list(inbred,aspect)]
out <- lapply(as.character(unique(annot_merge$inbred)),function(x){
    print(x)
    curGafFile <- paste("data/curated/",x,".curated.gaf",sep="")
    print(annot_merge[inbred==x])
    write_gaf(annot_merge[inbred==x,-c("B73v3","inbred","geneModel"),with=F],curGafFile)
})


