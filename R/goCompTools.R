library("ontologyIndex")
# library("ontologySimilarity")
library("data.table")
library("parallel")
library("ggplot2")
library("reshape2")
library("tools")
library("scales")

getJaccardIndex = function(set1Terms,set2Terms){
  return(length(intersect(set1Terms,set2Terms))/length(union(set1Terms,set2Terms)))
}



getAnnotComp <- function(data,go_obo){
  numSource = length(unique(data$source))
  out = list(intersect=as.integer(0),
             commUniq=as.integer(0),
             gomapUniq=as.integer(0),
             union=as.integer(0))
  if(numSource>1){
    comm_leaf_terms = unique(data[source=="Community"]$term_accession)
    gomap_leaf_terms = unique(data[source=="GOMAP"]$term_accession)
    comm_all_terms <- unique(unlist(go_obo$ancestors[comm_leaf_terms]))
    gomap_all_terms <- unique(unlist(go_obo$ancestors[gomap_leaf_terms]))
    out["intersect"] = as.integer(length(intersect(comm_all_terms,gomap_all_terms)))
    out["commUniq"] = as.integer(length(setdiff(comm_all_terms,gomap_all_terms)))
    out["gomapUniq"] = as.integer(length(setdiff(gomap_all_terms,comm_all_terms)))
    out["union"] = as.integer(length(union(gomap_all_terms,comm_all_terms)))
  }else{
    if(unique(data$source)=="GOMAP"){
      gomap_leaf_terms = unique(data$term_accession)
      gomap_all_terms <- unique(unlist(go_obo$ancestors[gomap_leaf_terms]))
      out["gomapUniq"] = as.integer(length(unique(gomap_all_terms)))
      out["union"] = as.integer(length(unique(gomap_all_terms)))
    }else if(unique(data$source)=="Community"){
      comm_leaf_terms = unique(data$term_accession)
      comm_all_terms <- unique(unlist(go_obo$ancestors[comm_leaf_terms]))
      out["commUniq"] = as.integer(length(unique(comm_all_terms)))
      out["union"] = as.integer(length(unique(comm_all_terms)))
    }
  }
  return(out)
}
