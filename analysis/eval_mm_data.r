source("code/cafa_eval.r")
source("code/aigo_eval.r")

#install.packages("pbmcapply")

sel_cols=c("db_object_symbol","term_accession","evidence_code","with","aspect","assigned_by")

go_file="go_data/plants.obo"
gold_file="gold/maize_v3.gold.gaf"

#check and read/load the go.obo file
go_obo = check_obo_data(go_file)

setdiff(gold_annots$term_accession,go_obo$id)
length(unique(gold_annots$term_accession))

#read the gold standard file
gold_annots = read_gaf(gold_file,sel_cols)
gold_annots = gaf_check_simple(go_obo,gold_annots)
gold_annots = gold_annots[grep("GRMZM|AC",db_object_symbol)]
gold_annots$with = as.numeric(gold_annots$with)
gold_annots[,with:=2]
setkey(gold_annots,db_object_symbol)
cores = 8

if(F){
  pred_annots <- read_gaf(infiles[3],sel_cols)
  pred_filt <- pred_annots[db_object_symbol %in% gold_annots$db_object_symbol]
  all_annots <- rbind(gold_annots,pred_filt)
  all_annots[,type:=ifelse(assigned_by=="MaizeGDB","gold","pred")]
  get_tool_cafa_eval(0,all_annots[db_object_symbol=="GRMZM2G312201"],go_obo,cores)
}

infiles=dir("maize.B73.AGPv3.r1/a.mm_gaf",pattern = "*.gaf.gz",full.names = T)

cafa_log <- lapply(infiles, function(infile){
    tmp_annots <- read_gaf(infile,sel_cols,rows = 100)
    tool_name = unique(tmp_annots$assigned_by)
    outfile = paste("tables/cafa_eval.",tool_name,".tsv",sep="")
    if(!file.exists(outfile)){
        cat("Processing ",infile,"and writing ",outfile,"\n")
        pred_annots <- read_gaf(infile,sel_cols)
        pred_annots = gaf_check_simple(go_obo,pred_annots)
        pred_filt <- pred_annots[db_object_symbol %in% gold_annots$db_object_symbol]
        rm(pred_annots)
        gc()
        
        all_annots <- rbind(gold_annots,pred_filt)
        all_annots[,type:=ifelse(assigned_by=="MaizeGDB","gold","pred")]
        cafa_eval_list <- lapply(seq(0,1.0,0.05),function(th){
            get_tool_cafa_eval(th,all_annots,go_obo,cores)
        })
        cafa_eval_dt <-  do.call(rbind,cafa_eval_list)
        write.table(cafa_eval_dt,outfile,quote = F,sep = "\t",row.names = F)    
    }
    else{
        cat("The oufile",outfile,"exists delete to remake","\n")
    }
})

aigo_log <- lapply(infiles, function(infile){
    tmp_annots <- read_gaf(infile,sel_cols,rows = 100)
    tool_name = unique(tmp_annots$assigned_by)
    outfile = paste("tables/aigo_eval.",tool_name,".tsv",sep="")
    if(!file.exists(outfile)){
        cat("Processing ",infile,"and writing ",outfile,"\n")
        pred_annots <- read_gaf(infile,sel_cols)
        pred_annots = gaf_check_simple(go_obo,pred_annots)
        pred_filt <- pred_annots[db_object_symbol %in% gold_annots$db_object_symbol]
        rm(pred_annots)
        gc()
        
        all_annots <- rbind(gold_annots,pred_filt)
        all_annots[,type:=ifelse(assigned_by=="MaizeGDB","gold","pred")]
        aigo_eval_list <- lapply(seq(0,1.0,0.05),function(th){
            get_tool_aigo_eval(th,all_annots,go_obo,cores)
        })
        aigo_eval_dt <-  do.call(rbind,aigo_eval_list)
        write.table(aigo_eval_dt,outfile,quote = F,sep = "\t",row.names = F)    
    }
    else{
        cat("The oufile",outfile,"exists delete to remake","\n")
    }
})
