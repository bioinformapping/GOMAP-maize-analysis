source("analysis/R/cafaEval.R")
source("analysis/R/plotTools.R")
library("yaml")
library("data.table")
library("ggplot2")
library("naturalsort")
library("gtable")
library("ontologySimilarity")
#install.packages("Rgraphviz")
library("ontologyPlot")

go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yml")

commInbreds = unlist(lapply(config$predicted,function(x){
  if(x$source=="Community"){
    unique(x$inbred)
  }
}))
commInbreds

predData_list <- lapply(config$predicted,function(x){
  inputGaf <- file.path("data/clean",x$file)
  inbred = x$inbred
  if(inbred %in% commInbreds){
    input_gaf_data <- read_gaf(inputGaf)
    input_gaf_data[,inbred:=x$inbred]
    input_gaf_data[,source:=x$source]
    input_gaf_data
  }else{
    NULL
  }
})

predData_dt <- rbindlist(predData_list)
predData_dt
dcast(predData_dt[,list(count=length(unique(db_object_symbol))),list(inbred,aspect,source)],aspect~inbred+source)

curateData_list <- lapply(config$curated,function(x){
  inputGaf <- file.path("data/curated/",x$file)
  if(x$inbred %in% commInbreds){
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:="curated"]
  input_gaf_data
  }else{
    NULL
  }
})

curatedData_dt <- rbindlist(curateData_list)


curatedData_dt[,annotType:="curated"]
predData_dt[,annotType:="predicted"]

gafData_dt <- rbind(predData_dt,curatedData_dt)
gafData_dt[,list(count=.N),by=list(inbred,source,annotType)]

communityGenes <- unique(predData_dt[source=="Community",.(inbred,db_object_symbol,aspect)])
curatedGenes <- merge.data.table(communityGenes,unique(curatedData_dt[,.(inbred,db_object_symbol,aspect)]))
curatedGenes

tmpFilt <- gafData_dt[,list(count=length(unique(source))),by=list(inbred,aspect,db_object_symbol)][count==3]
tmpIdx=1
tmpFilt[tmpIdx]
tmpGene=tmpFilt[tmpIdx]$db_object_symbol
tmpAspect=tmpFilt[tmpIdx]$aspect

term_list <- lapply(split(gafData_dt[db_object_symbol == tmpGene & aspect==tmpAspect],by = c("source")),function(x){
  get_ancestors(go_obo,x$term_accession)
})

geneOntoPlot <- onto_plot(ontology = go_obo,term_sets = term_list,shape-"box",fixedsize = F,label = label_by_term_set,fillcolor=colour_by_frequency)
tmpOut <- file.path("paper/figures/dot",paste0(tmpGene,"-",tmpAspect,".dot"))
write_dot(geneOntoPlot,file = tmpOut)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("DiagrammeR")

V <- letters[1:10]
M <- 1:4
set.seed(225)
g1 <- randomGraph(V, M, 0.4)
g1 <- agopen(g1,name = "test")
#png("test.png") #,width = 5,height = 5,units = "in",res = 300)
z = graphLayout(g1)
z

Rgraphviz::
