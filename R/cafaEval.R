library("data.table")
library("ggplot2")
library("reshape2")
library("tools")
library("parallel")
library("sets")

source("analysis/R/oboTools.R")
source("analysis/R/genUtils.R")
source("analysis/R/gafTools.R")
source("analysis/R/getNonredDataset.R")
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

getCafaEval  = function(data){

}

getCafaPr <- function(pred_leaf_terms,inbred_txt,aspect_txt,db_object_symbol_txt,go_obo,curatedData_dt){
  gold_leaf_terms = curatedData_dt[inbred==inbred_txt & aspect == aspect_txt & db_object_symbol==db_object_symbol_txt]$term_accession
  out_list = list(pr=0.0,num_gold=0.0,num_pred=0.0,evalType="CAFA")
  if(NROW(gold_leaf_terms)>0 & NROW(pred_leaf_terms)>0){
    prop_pred_terms = setdiff(get_ancestors(go_obo,pred_leaf_terms),excl_go)
    prop_gold_terms = setdiff(get_ancestors(go_obo,gold_leaf_terms),excl_go)
    cafa_pr = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_pred_terms)
    out_list[["pr"]] = cafa_pr
    out_list[["num_gold"]] = as.numeric(NROW(prop_gold_terms))
    out_list[["num_pred"]] = as.numeric(NROW(prop_pred_terms))
  }else if(NROW(tmpPredicted)==0){
    cat("=============================================\nError\n")
    cat(inbred_txt,aspect_txt,db_object_symbol_txt,"\n")
    cat("=============================================\n")
  }

  return(out_list)
}

getCafaRc <- function(gold_leaf_terms,inbred_txt,aspect_txt,db_object_symbol_txt,go_obo,pred_data){
  pred_leaf_terms = pred_data[inbred==inbred_txt & aspect == aspect_txt & db_object_symbol==db_object_symbol_txt]$term_accession
  out_list = list(rc=0.0,num_gold=0.0,num_pred=0.0,evalType="CAFA")
  if(NROW(gold_leaf_terms)>0 & NROW(pred_leaf_terms)>0){
    prop_pred_terms = setdiff(get_ancestors(go_obo,pred_leaf_terms),excl_go)
    prop_gold_terms = setdiff(get_ancestors(go_obo,gold_leaf_terms),excl_go)
    cafa_pr = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_gold_terms)
    out_list[["rc"]] = cafa_pr
    out_list[["num_gold"]] = as.numeric(NROW(prop_gold_terms))
    out_list[["num_pred"]] = as.numeric(NROW(prop_pred_terms))
  }else if(NROW(gold_leaf_terms)==0){
    cat("=============================================\nError\n")
    cat(inbred_txt,aspect_txt,db_object_symbol_txt,"\n")
    cat("=============================================\n")
  }

  return(out_list)
}

getGamerEval <- function(data){
  #print(data)
  tmpCurated <- data[annotType=="curated"]
  tmpPredicted <- data[annotType=="predicted"]

  out_list = list(pr=0.0,rc=0.0,num_gold=0.0,num_pred=0.0,evalType="GAMER")
  if(NROW(tmpCurated)>0 & NROW(tmpPredicted)>0){

    gold_leaf_terms <- minimal_set(go_obo,tmpCurated$term_accession)
    pred_leaf_terms <- minimal_set(go_obo,tmpPredicted$term_accession)

    gamerPr <- calcGamerPr(pred_leaf_terms,gold_leaf_terms)
    gamerRc <- calcGamerRc(pred_leaf_terms,gold_leaf_terms)

    out_list[["pr"]] = gamerPr
    out_list[["rc"]] = gamerRc
    out_list[["num_gold"]] = as.numeric(NROW(tmpCurated))
    out_list[["num_pred"]] = as.numeric(NROW(tmpPredicted))

  }else if(NROW(tmpCurated)>0 & NROW(tmpPredicted)==0){
    out_list[["num_gold"]] = as.numeric(NROW(tmpCurated))
  }

  return(out_list)
}

calcGamerPr <- function(pred_terms,gold_terms){
  pr_list <- sapply(pred_terms,function(x) sapply(gold_terms,function(y){
    #prop_pred_terms = setdiff(get_ancestors(go_obo,x),excl_go)
    prop_pred_terms = get_ancestors(go_obo,x)
    #prop_gold_terms = setdiff(get_ancestors(go_obo,y),excl_go)
    prop_gold_terms = get_ancestors(go_obo,y)
    if(length(prop_pred_terms)>0 & length(prop_gold_terms)>0){
      tmp_pr = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_pred_terms)
    }else{
      tmp_pr=0
    }
    tmp_pr
  }))
  gamerPr <- mean(unlist(pr_list))
  return(gamerPr)
}

calcGamerRc <- function(pred_terms,gold_terms){
  rc_list <- sapply(gold_terms,function(x) sapply(pred_terms,function(y){
    # prop_gold_terms = setdiff(get_ancestors(go_obo,x),excl_go)
    # prop_pred_terms = setdiff(get_ancestors(go_obo,y),excl_go)
    prop_gold_terms = get_ancestors(go_obo,x)
    prop_pred_terms = get_ancestors(go_obo,y)
    if(length(prop_pred_terms)>0 & length(prop_gold_terms)>0){
      tmp_rc = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_gold_terms)
    }else{
      tmp_rc=0
    }
    tmp_rc
  }))
  gamerRc <- mean(unlist(rc_list))
  return(gamerRc)
}

getAigoEval <- function(data){
  #print(data)
  tmpCurated <- data[annotType=="curated"]
  tmpPredicted <- data[annotType=="predicted"]

  out_list = list(pr=0.0,rc=0.0,num_gold=0.0,num_pred=0.0,evalType="GAMER")
  if(NROW(tmpCurated)>0 & NROW(tmpPredicted)>0){

    gold_leaf_terms <- minimal_set(go_obo,tmpCurated$term_accession)
    pred_leaf_terms <- minimal_set(go_obo,tmpPredicted$term_accession)

    gamerPr <- calcAigoPr(pred_leaf_terms,gold_leaf_terms)
    gamerRc <- calcAigoRc(pred_leaf_terms,gold_leaf_terms)

    out_list[["pr"]] = gamerPr
    out_list[["rc"]] = gamerRc
    out_list[["num_gold"]] = as.numeric(NROW(tmpCurated))
    out_list[["num_pred"]] = as.numeric(NROW(tmpPredicted))

  }else if(NROW(tmpCurated)>0 & NROW(tmpPredicted)==0){
    out_list[["num_gold"]] = as.numeric(NROW(tmpCurated))
  }

  return(out_list)
}

calcAigoPr <- function(pred_terms,gold_terms){
  pr_list <- sapply(pred_terms,function(x) {
    tmp_prs <- sapply(gold_terms,function(y){
    prop_pred_terms = get_ancestors(go_obo,x)
    prop_gold_terms = get_ancestors(go_obo,y)
    if(length(prop_pred_terms)>0 & length(prop_gold_terms)>0){
      tmp_pr = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_pred_terms)
    }else{
      tmp_pr=0
    }
    tmp_pr
  })
    return(max(tmp_prs))
  })
  gamerPr <- mean(unlist(pr_list))
  return(gamerPr)
}

calcAigoRc <- function(pred_terms,gold_terms){
  rc_list <- sapply(gold_terms,function(x) {
    tmp_rcs <- sapply(pred_terms,function(y){
    prop_gold_terms = get_ancestors(go_obo,x)
    prop_pred_terms = get_ancestors(go_obo,y)
    if(length(prop_pred_terms)>0 & length(prop_gold_terms)>0){
      tmp_rc = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_gold_terms)
    }else{
      tmp_rc=0
    }
    tmp_rc
  })
    return(max(tmp_rcs))
  })
  gamerRc <- mean(unlist(rc_list))
  return(gamerRc)
}
?get_ancestors


