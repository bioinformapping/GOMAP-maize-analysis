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

root_go=c("GO:0008150","GO:0005575","GO:0003674")
other_go=c("GO:0005515")
excl_go =  c(root_go,other_go)

ont_conv = list(C="CC",F="MF",P="BP")

if(F){
  tmp_data <- fread("zcat maize.B73.AGPv3.r1/a.mm_gaf/fanngo-0.0.gaf.gz",nrows = 10)
  tmp_data <- read_gaf("maize.B73.AGPv3.r1/a.mm_gaf/fanngo-0.0.gaf.gz",rows = 10)
  tmp_data
}

get_tool_aigo_eval <- function(th,all_annots,go_obo,cores){
    tool_name = unique(all_annots[assigned_by!="MaizeGDB"]$assigned_by)
    annot_source = unique(all_annots$assigned_by)
    all_annots[,num_types:=length(unique(type)),by=list(db_object_symbol,aspect)]
    
    cat("\nProcessing threshold:",th,"for",tool_name,"\n")
    filt_annots = all_annots[with>=th & !term_accession %in% excl_go]
    setkey(filt_annots,db_object_symbol)
    
    all_genes = unique(filt_annots$db_object_symbol)
    common_genes = intersect(all_genes,filt_annots[type=="pred"]$db_object_symbol)
    gold_genes = setdiff(all_genes,filt_annots[type=="pred"]$db_object_symbol)
    
    out_list <- pbmclapply(all_genes,function(in_gene){
    # out_list <- lapply(all_genes,function(in_gene){
        gene_annots=filt_annots[in_gene]
        aspects=unique(gene_annots[type=="gold"]$aspect)
        aspect_aigo_list <- lapply(aspects,function(in_aspect){
            gene_aspect_annots = gene_annots[aspect==in_aspect]
            if(NROW(gene_aspect_annots[type=="gold"])>0){
                gold_leaf_terms = minimal_set(go_obo,gene_aspect_annots[type=="gold"]$term_accession)
                # aigo_evals$num_gold = length(gold_leaf_terms)
                
                if(NROW(gene_aspect_annots[type=="pred"])>0){
                    pred_leaf_terms = minimal_set(go_obo,gene_aspect_annots[type=="pred"]$term_accession)
                    aigo_evals = get_aigo_eval(pred_leaf_terms,gold_leaf_terms)
                    val_dt <- cbind(gene=in_gene,aspect=in_aspect,aigo_evals)
                }else{
                    gold_term=rep(gold_leaf_terms,each=2)
                    metric=c("pr","rc")
                    val_dt <- cbind(gene=in_gene,aspect=in_aspect,gold_term,pred_term="",value=0,metric)
                }
                data.table(val_dt)
            }
        })
        aspect_aigo_dt = rbindlist(aspect_aigo_list,use.names = T)
    # })
    },mc.cores = cores)
    aigo_eval_dt <-  rbindlist(out_list,use.names = T)
    out_dt = data.table(cbind(aigo_eval_dt,tool_name,th))
    return(out_dt)
}

get_aigo_eval <- function(pred_leaf_terms,gold_leaf_terms){
    
    pr_list <- lapply(pred_leaf_terms,function(pred_term){
        prop_pred_terms = get_ancestors(go_obo,pred_term)
        pr_val_list <- lapply(gold_leaf_terms,function(gold_term){
            prop_gold_terms = get_ancestors(go_obo,gold_term)
            pr_val = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_pred_terms)
            cbind(gold_term,value=pr_val)
        })
        
        data.table(cbind(pred_term,do.call(rbind,pr_val_list),metric="pr"))
    })
    pr_dt <- rbindlist(pr_list,use.names = T)
    
    rc_list <- lapply(gold_leaf_terms,function(gold_term){
        prop_gold_terms = get_ancestors(go_obo,gold_term)
        pr_val_list <- lapply(pred_leaf_terms,function(pred_term){
            prop_pred_terms = get_ancestors(go_obo,pred_term)
            rc_val = length(intersect(prop_pred_terms,prop_gold_terms))/length(prop_gold_terms)
            cbind(pred_term,value=rc_val)
        })
        
        data.table(cbind(gold_term,do.call(rbind,pr_val_list),metric="rc"))
    })
    rc_dt <- rbindlist(rc_list,use.names = T)
    out_dt <- rbind(pr_dt,rc_dt)
    return(out_dt)
}


get_tool_score = function(th,tool_data,gold_data,go_obo){
    back_buff = paste(rep(" ",16),collapse = "")
    tool_name = unique(tool_data$assigned_by)
    
    cat("\nProcessing threshold:",th,"for",tool_name,"\n")
    tool_data_filt = tool_data[with>=as.numeric(th)]
    common_genes = unique(intersect(tool_data_filt$db_object_symbol,gold_data$db_object_symbol))
    tool_data_eval = tool_data_filt[db_object_symbol %in% common_genes]
    gold_only_genes = unique(setdiff(gold_data$db_object_symbol,tool_data_filt$db_object_symbol))
    
    if(NROW(tool_data_eval)>0){
        tool_data_list = gaf2list(tool_data_filt)
        rownames(tool_data_filt) = NULL
        
        cat("Number of Rows ",NROW(tool_data_filt),"\n\n")
        cat(back_buff,"\n")
        
        unit_perc = 1
        unit_size = ceiling(nrow(gold_data) / (100/unit_perc))
        
        
        tool_score_perf = gold_data[sort(common_genes),get_aigo_score(db_object_symbol,aspect,.SD,unit_size,.I,unit_perc,tool_data_list,go_obo,tool_name),by=list(db_object_symbol,aspect)]
        #tool_score_perf = tool_data_filt[sort(com_genes),get_annot_sim_score(db_object_symbol,term_accession,.SD,unit_size,.I,unit_perc),by=list(db_object_symbol,term_accession)]
        tool_score_filt = tool_score_perf[tool_score_perf[, .I[hPr > -1 & hRc > -1 & max_hPr > -1 & max_hRc > -1],]]
    }
    if(length(gold_only_genes)>0){
        gold_only_score = cbind(gold_data[gold_only_genes,.(db_object_symbol,aspect)],score=0,hPr=0,hRc=0,max_hPr=0,max_hRc=0)
    }
    
    all_score <- rbind(tool_score_filt,gold_only_score)
    all_score = cbind(all_score,th)
    return(all_score)
}

get_aigo_score = function(gene,aspect,data,unit_size,idx,unit_perc,tool_data_list,go_obo,tool_name){
    
    print_dt_progress(unit_size,idx,unit_perc,tool_name)
    out_list = list(score=-1,hPr=-1,hRc=-1,max_hPr=-1,max_hRc=-1)
    #gold_ans = unique(unlist(go_obo$ancestors[gold_term]))
    
    tool_terms = unlist(unique(tool_data_list[[gene]][[aspect]]))
    gold_terms = unlist(unique(data$term_accession))
    
    #if(nrow(gold_dt)>0 | !is.null(gold_terms)){
    
    if(length(gold_terms) > 0 & length(tool_terms)>0){
        
        gold_terms = minimal_set(go_obo,unlist(unique(data$term_accession)))
        tool_terms = minimal_set(go_obo,unlist(tool_data_list[[gene]][[aspect]]))
        
        hPrs = lapply(tool_terms,function(term){
            get_hPr(term,gold_terms,go_obo)
        })
        
        out_list$hPr = mean(unlist(lapply(hPrs,mean)))
        out_list$max_hPr = mean(unlist(lapply(hPrs,max)))
        
        hRcs = lapply(gold_terms,function(term){
            get_hRc(term,tool_terms,go_obo)
        })
        
        out_list$hRc = mean(unlist(lapply(hRcs,mean)))
        out_list$max_hRc = mean(unlist(lapply(hRcs,max)))
        
        if(F & aspect!="C"){
            print(aspect)
            print(gold_terms)
            print(tool_terms)
            print(hPrs)
            print(hRcs)
        }
    }
    
    out_list = lapply(out_list,function(x){
        ifelse(is.na(x),-1,x)
    })
    return(out_list)
}

gaf2list = function(gaf){
    
    tmp2list = function(gene,aspect,terms){
        out = list()
        out[[gene]] = list()
        out[[gene]][[aspect]] = terms
        return(out)
    }
    
    tmp_gaf = gaf[,list(terms=tmp2list(db_object_symbol,aspect,unique(term_accession))),by=list(db_object_symbol,aspect)]
    
    list2list = function(gene,lists){
        out = list()
        out[[gene]] = unlist(lists,recursive = F)
        out
    }
    
    final_gaf = tmp_gaf[,list(final_list=list2list(db_object_symbol,terms)),by=db_object_symbol]
    final_list = final_gaf$final_list
    names(final_list) = final_gaf$db_object_symbol
    return(final_list)
}

get_hPr_hRc = function(gold_ans,tool_term){
    #gold_ans = unique(unlist(go_obo$ancestors[gold_term]))
    tool_ans =  unique(unlist(go_obo$ancestors[tool_term]))
    hPr = length(intersect(tool_ans,gold_ans))/length(tool_ans)
    hRc = length(intersect(tool_ans,gold_ans))/length(gold_ans)
    return(list(hPr=hPr,hRc=hRc))
}

get_hPr = function(tool_term,gold_terms,go_obo){
    tool_ns = go_obo$namespace[[tool_term]]
    tool_ans = unique(unlist(go_obo$ancestors[tool_term]))
    #tool_ans = tool_ans[go_obo$namespace[tool_ans] == tool_ns]
    tool_ans = intersect(tool_ans,go_obo$ns2go[[tool_ns]])
    tmp_list = lapply(gold_terms,function(term){
        gold_ans = unique(unlist(go_obo$ancestors[term]))
        #gold_ans = gold_ans[go_obo$namespace[gold_ans] == tool_ns]
        gold_ans = intersect(gold_ans,go_obo$ns2go[[tool_ns]])
        hPr = length(intersect(tool_ans,gold_ans))/length(tool_ans)
        return(hPr)
    })
    return(unlist(tmp_list))
}
get_hRc = function(gold_term,tool_terms,go_obo){
    tool_ns = go_obo$namespace[[gold_term]]
    gold_ans = unique(unlist(go_obo$ancestors[gold_term]))
    #gold_ans = gold_ans[go_obo$namespace[gold_ans] == tool_ns]
    gold_ans = intersect(gold_ans,go_obo$ns2go[[tool_ns]])
    tmp_list = lapply(tool_terms,function(term){
        tool_ans = unique(unlist(go_obo$ancestors[term]))
        #tool_ans = tool_ans[go_obo$namespace[tool_ans] == tool_ns]
        tool_ans = intersect(tool_ans,go_obo$ns2go[[tool_ns]])
        hRc = length(intersect(tool_ans,gold_ans))/length(gold_ans)
        return(hRc)
    })
    return(unlist(tmp_list))
}