library("ggplot2")
library("data.table")

cores = 8

infiles=dir("tables",pattern = "cafa.*tsv",full.names = T)


all_data <- lapply(infiles,fread,sep="\t")
data <- rbindlist(all_data)[!(aspect=="C" & tool_name=="FANN-GO")]
data[aspect=="C"]

data_summary <-  data[,list(pr=mean(cafa_pr),rc=mean(cafa_rc)),by=list(aspect,th,tool_name)]
data_summary
pr <- data[num_pred>0,list(pr=mean(cafa_pr)),by=list(aspect,th,tool_name)]
rc <- data[num_gold>0,list(rc=mean(cafa_rc)),by=list(aspect,th,tool_name)]
plot_data <- merge(pr,rc,by =c("aspect","th","tool_name"),all = T)
plot_data[,fscore:=2*(pr*rc)/(pr+rc)]
plot_data

fmax_data <- plot_data
fmax_data[is.na(plot_data)]=0
fmax_cafa <- fmax_data[,list(fscore=max(fscore),th=th[which.max(fscore)],pr=pr[which.max(fscore)],rc=rc[which.max(fscore)]),by=list(aspect,tool_name)]
fmax_cafa

p <- ggplot(plot_data,aes(x=rc,y=pr,color=tool_name))
p <- p + geom_line() + geom_point(data=fmax_cafa,color="#000000",size=3) + geom_point()
p <- p + facet_grid(aspect~.) + labs(x="Recall",y="Precision")
p <- p + theme(legend.position = "bottom")
p
ggsave("plots/mm-cafa-pr-rc.png",width = 5,height = 6)

p <- ggplot(plot_data,aes(x=th,y=fscore,color=tool_name))
p <- p + geom_line() + geom_point(data=fmax_cafa,color="#000000",size=3) + geom_point()
p <- p + facet_grid(aspect~.) + labs(x="Score",y="F-score")
p <- p + theme(legend.position = "bottom")
p
ggsave("plots/mm-cafa-fscore.png",width = 6,height = 5)



data[ num_gold==2 & num_pred == 3][order(th,decreasing = T)]

all_annots[db_object_symbol=="GRMZM5G842058" & with >= 0.75 & aspect == "P"]

pred_terms=minimal_set(go_obo,all_annots[db_object_symbol=="GRMZM5G842058" & with >= 0.75 & aspect == "P" & assigned_by=="PANNZER"]$term_accession)
gold_terms=minimal_set(go_obo,all_annots[db_object_symbol=="GRMZM5G842058" & with >= 0.75 & aspect == "P" & assigned_by=="MaizeGDB"]$term_accession)
pred_terms
gold_terms

pr_out <- sapply(pred_terms,function(x){
  sapply(gold_terms,function(y){
    a <- setdiff(get_ancestors(go_obo,x),excl_go)
    b <- setdiff(get_ancestors(go_obo,y),excl_go)  
    length(intersect(a,b))/length(a)
  })
})
write.table(pr_out,"tables/aigo-pr-eg.tsv",sep = "\t")

rc_out <- sapply(gold_terms,function(x){
  sapply(pred_terms,function(y){
    a <- setdiff(get_ancestors(go_obo,x),excl_go)
    b <- setdiff(get_ancestors(go_obo,y),excl_go)  
    length(intersect(a,b))/length(a)
  })
})
write.table(rc_out,"tables/aigo-rc-eg.tsv",sep = "\t")