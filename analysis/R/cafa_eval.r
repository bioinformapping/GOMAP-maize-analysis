library("data.table")
library("ggplot2")
library("reshape2")
library("tools")
library("parallel")
library("sets")
library("pbmcapply")

source("R/obo_tools.r")
source("R/gen_utils.r")
source("R/gaf_tools.r")
source("R/get_nr_dataset.r")
# source("code/aigo_tools.r")

root_go  = c("GO:0008150","GO:0005575","GO:0003674")
other_go = c("GO:0005515")
excl_go  = c(root_go,other_go)

ont_conv = list(C="CC",F="MF",P="BP")

get_tool_cafa_eval = function(th,all_annots,go_obo,cores){
    # back_buff = paste(rep(" ",16),collapse = "")
    tool_name = unique(all_annots[assigned_by!="MaizeGDB"]$assigned_by)
    annot_source = unique(all_annots$assigned_by)
    all_annots[,num_types:=length(unique(type)),by=list(db_object_symbol,aspect)]
    
    cat("\nProcessing threshold:",th,"for",tool_name,"\n")
    filt_annots = all_annots[with>=th & ! term_accession %in% excl_go]
    setkey(filt_annots,db_object_symbol)
    
    all_genes = unique(filt_annots$db_object_symbol)
    common_genes = intersect(all_genes,filt_annots[type=="pred"]$db_object_symbol)
    gold_genes = setdiff(all_genes,filt_annots[type=="pred"]$db_object_symbol)
    
    out_list <- pbmclapply(all_genes,function(in_gene){
    #out_list <- mclapply(all_genes,function(in_gene){
        gene_annots=filt_annots[in_gene]
        aspects=unique(gene_annots[type=="gold"]$aspect)
        aspect_cafa_list <- lapply(aspects,function(in_aspect){
            gene_aspect_annots = gene_annots[aspect==in_aspect]
            cafa_evals = as.data.table(list(cafa_pr=0,cafa_rc=0,num_gold=0,num_pred=0))
            if(NROW(gene_aspect_annots[type=="gold"])>0){
                gold_leaf_terms = minimal_set(go_obo,gene_aspect_annots[type=="gold"]$term_accession)
                cafa_evals$num_gold = length(gold_leaf_terms)
                if(NROW(gene_aspect_annots[type=="pred"])>0){
                    pred_leaf_terms = minimal_set(go_obo,gene_aspect_annots[type=="pred"]$term_accession)
                    cafa_evals = get_cafa_eval(pred_leaf_terms,gold_leaf_terms)
                }
            }
            aspect_cafa_dt <- cbind(aspect=in_aspect,as.data.table(cafa_evals))
        })
        aspect_cafa_dt = do.call(rbind,aspect_cafa_list)
        gene_aspect_cafa_dt <- cbind(aspect_cafa_dt,gene=in_gene)
    #})
    },mc.cores = cores)
    
    cafa_eval_dt <-  do.call(rbind,out_list)
    cafa_eval_dt
    out_dt = cbind(cafa_eval_dt,tool_name,th)
    return(out_dt)
}

get_cafa_eval <- function(pred_leaf_terms,gold_leaf_terms){
    prop_pred_terms = setdiff(get_ancestors(go_obo,pred_leaf_terms),excl_go)
    prop_gold_terms = setdiff(get_ancestors(go_obo,gold_leaf_terms),excl_go)
    cafa_pr = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_pred_terms)
    cafa_rc = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_gold_terms)
    cafa_pr = ifelse(is.nan(cafa_pr),0,cafa_pr)
    cafa_rc = ifelse(is.nan(cafa_rc),0,cafa_rc)
    out_list = list(cafa_pr=cafa_pr,cafa_rc=cafa_rc,num_gold=length(gold_leaf_terms),num_pred=length(pred_leaf_terms))
    return(out_list)
}