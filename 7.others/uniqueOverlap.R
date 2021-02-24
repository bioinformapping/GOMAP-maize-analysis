source("analysis/R/plotTools.R")
source("analysis/R/gafTools.R")
source("analysis/R/oboTools.R")
source("analysis/R/goCompTools.R")
library("data.table")
library("naturalsort")
library("kableExtra")
library("gtable")
library("yaml")
library("DiagrammeR")

go_obo <- check_obo_data("data/go/go.obo")
config <- read_yaml("config.yml")
config

predData_list <- lapply(config$predicted[grep("B73v4|PH207",config$predicted)],function(x){
  inputGaf <- file.path("data/clean",x$file)
  inbred = x$inbred
  input_gaf_data <- read_gaf(inputGaf)
  input_gaf_data[,inbred:=x$inbred]
  input_gaf_data[,source:=x$source]
  input_gaf_data
})

predData_dt <- rbindlist(predData_list )
unique(predData_dt[,.(inbred,source)])
predData_dt[,source:=factor(source,levels = c("GOMAP","Community"),ordered = T)]


venn_term_acc = split(predData_dt[,.(inbred,aspect,source,term_accession)],by=c("inbred","aspect","source"),flatten=FALSE,keep.by = F,sorted = T)

distinctB73TermList = lapply(names(venn_term_acc$B73v4),function(x){
  allCommunityTerms = unique(unlist(get_ancestors(go_obo,venn_term_acc$B73v4[[x]]$Community$term_accession)))
  allGomapTerms = unique(unlist(get_ancestors(go_obo,venn_term_acc$B73v4[[x]]$GOMAP$term_accession)))
  uniqCommunityTerms = setdiff(allCommunityTerms,allGomapTerms)
  communityPredDt <- cbind(term_accession=uniqCommunityTerms,source="Community",aspect=x)
  uniqGomapTerms = setdiff(allGomapTerms,allCommunityTerms)
  gomapPredDt <- cbind(term_accession=uniqGomapTerms,source="GOMAP",aspect=x)
  data.table(rbind(communityPredDt,gomapPredDt))
})
distinctB73Terms= unique(rbindlist(distinctB73TermList))

gomapNodes = distinctB73Terms[aspect=="C" & source=="GOMAP"]$term_accession
communityNodes = distinctB73Terms[aspect=="C" & source=="Community"]$term_accession
allLeafNodes = unique(unlist(go_obo$ancestors[c(gomapNodes,communityNodes)]))
allNodes = unlist(go_obo$ancestors[allLeafNodes])

nodes_1 = data.table(cbind(id=1:length(allNodes),label=NA,type="both",style="filled",term=allNodes))
nodes_1[,id:=as.integer(id)]
nodes_1[term%in%gomapNodes,type:="gomap"]
nodes_1[term%in%communityNodes,type:="community"]
allNodes_id = 1:length(allNodes)
names(allNodes_id) = allNodes
RColorBrewer::brewer.pal(n = 3,name = "Accent")

nodeColors = RColorBrewer::brewer.pal(n = 3,name = "Accent")
names(nodeColors) <- c("gomap","both","community")
nodes_1[,fillcolor:=nodeColors[type]]
nodeColors



node_parents = lapply(allNodes,function(x){
  parent_nodes = go_obo$parents[[x]]
  if(length(parent_nodes)>0){
    out <- data.table(cbind(to=allNodes_id[[x]],from=unlist(allNodes_id[parent_nodes])))  
    return(out)
  }else{
    NULL
  }
  
})
edges_1 = rbindlist(node_parents)
edges_1

nodes_1

g <- create_graph(nodes_df = nodes_1,edges_df = edges_1,directed = T)

g %>%
  set_node_attr_to_display() %>%
  render_graph()
  #export_graph(file_name = "test.png")



