library("data.table")
library("ggplot2")
library("ontologyIndex")
library(datapasta)
library(scales)
source("code/gaf_tools.r")
source("code/obo_tools.r")
source("code/aigo_eval.r")

infiles <- c(dir("maize.B73.AGPv3.r1/d.non_red_gaf",".*gaf.gz",full.names = T),dir("maize.B73.AGPv3.r1/e.agg_data","nr.gaf.gz",full.names = T))
infiles

datasets <- fread("public_datasets.tsv")
datasets$Type = factor(datasets$Type, levels=datasets$Type)
datasets$dataset = factor(datasets$dataset, levels=datasets$dataset)
datasets

infiles <- unlist(sapply(datasets$file, function(x){
    infiles[grep(x,infiles)]
}))
infiles

refset <- fread("maize.B73.AGPv3.r1/maize_refset.txt",header = F)$V1
num_genes <- NROW(refset)

go_file="obo/go.obo"
gold_file="gold/maize_v3.gold.nr.gaf"

#check and read/load the go.obo file
go_obo = check_obo_data(go_file)

gold_file="gold/maize_v3.gold.nr.gaf"
#read the gold standard file
gold_annots = read_gaf(gold_file)
gold_annots = gaf_check_simple(go_obo,gold_annots)
gold_annots = gold_annots[grep("GRMZM|AC",db_object_symbol)]
gold_annots$with = as.numeric(gold_annots$with)
gold_annots[,with:=2]
gold_annots[,type:="gold"]
setkey(gold_annots,db_object_symbol)
cores = 4

data_list <- lapply(infiles,read_gaf)
all_data <- rbindlist(data_list,use.names = T)
all_data[,with:=as.numeric(with)]
all_data[,with:=1]
cov_dt <- all_data[,list(value=length(unique(db_object_symbol))/num_genes*100,metric="coverage"),by=list(assigned_by,aspect)]
cov_dt

all_aigo_eval <- lapply(datasets$dataset,function(x){
    th=0
    tmp_annots = all_data[assigned_by==x]
    tmp_annots[,type:="pred"]
    all_annots = rbind(tmp_annots,gold_annots)
    get_tool_aigo_eval(th,all_annots,go_obo,cores)
})

aigo_eval_dt <- rbindlist(all_aigo_eval)
aigo_eval_dt[,value:=as.numeric(value)]
write.table(aigo_eval_dt,"tables/aigo_eval.public.nr.tsv",quote = F,sep = "\t",row.names = F)
aigo_eval_dt <- fread("tables/aigo_eval.public.nr.tsv")

tmp_pr <- aigo_eval_dt[metric=="pr",list(aigo_pr=mean(value)),by=list(gene,aspect,tool_name)]
tmp_rc <- aigo_eval_dt[metric=="rc",list(aigo_rc=mean(value)),by=list(gene,aspect,tool_name)]
aigo_eval_dt <- merge(tmp_pr,tmp_rc,by=c("gene","aspect","tool_name"))
aigo_eval_dt[,aigo_fscore:=ifelse((aigo_pr+aigo_rc>0),2*(aigo_pr*aigo_rc)/(aigo_pr+aigo_rc),0)]
aigo_eval_dt
fmax_dt <- aigo_eval_dt[,.(value=mean(aigo_fscore),metric="fmax"),by=list(tool_name,aspect)]
colnames(fmax_dt) <- colnames(cov_dt)

fscore_dt <- aigo_eval_dt[,list(value=mean(aigo_fscore),metric="fscore"),by=list(tool_name,aspect)]
num_annot_dt <- all_data[,list(value=.N,metric="annotations"),by=list(assigned_by,aspect)]

get_spec <- function(in_list){
    anc_list <- mclapply(in_list,get_ancestors,ontology=go_obo,mc.cores = cores)
    anc_len <- lapply(anc_list,length)
    spec <- mean(unlist(anc_len))
    return(spec)
}
spec_dt <- all_data[,list(value=get_spec(term_accession),metric="specificity"),by=list(assigned_by,aspect)]
spec_dt
all_metrics <- rbind(cov_dt,num_annot_dt,spec_dt,fmax_dt)
plot_data <- merge(all_metrics,datasets,by.x = "assigned_by",by.y = "dataset")
plot_data$assigned_by = factor(plot_data$assigned_by,levels = datasets$dataset)
plot_data$Type = factor(plot_data$Type,levels=unique(datasets$Type))
plot_data[assigned_by=="Aggregate"]$assigned_by="maize-GAMER"


public_p

public_p <- ggplot(plot_data,aes(x=assigned_by,y=value,fill=Type))
public_p <- public_p + geom_bar(stat="identity",color="black")
public_p <- public_p + facet_grid(metric~aspect,scales = "free") + labs(x="Dataset")
public_p <- public_p + scale_fill_manual(values=unique(datasets$color)) + scale_y_continuous(labels = comma)
public_p <- public_p + theme_bw() + theme(axis.title.y = element_blank(),strip.text = element_blank())
public_p
ggsave("plots/public-assess.png",width = 12,height = 5,units = "in",dpi = 300)
?scale_fill_manual()
