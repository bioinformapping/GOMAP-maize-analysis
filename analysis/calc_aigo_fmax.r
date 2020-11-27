library("ggplot2")
library("data.table")

infiles=dir("tables",pattern = "aigo_eval.*tsv",full.names = T)
all_data <- lapply(infiles,fread,sep="\t")
data <- rbindlist(all_data)[!(aspect=="C" & tool_name=="FANN-GO")]
data[gene=="AC149810.2_FG004" & tool_name == "PANNZER"]

tmp_pr <- data[metric=="pr",list(aigo_pr=mean(value)),by=list(gene,aspect,th,tool_name)]
tmp_rc <- data[metric=="rc",list(aigo_rc=mean(value)),by=list(gene,aspect,th,tool_name)]
aigo_data <- merge(tmp_pr,tmp_rc,by=c("gene","aspect","th","tool_name"))

aigo_data[,aigo_fscore:=2*(aigo_pr*aigo_rc)/(aigo_pr+aigo_rc)]
aigo_data[gene=="AC149810.2_FG004"]
aigo_data[is.nan(aigo_fscore),aigo_fscore:=0]

data_summary <-  aigo_data[,list(pr=mean(aigo_pr),rc=mean(aigo_rc),fscore=mean(aigo_fscore)),by=list(aspect,th,tool_name)]
plot_data <- data_summary
fmax_data <- plot_data
fmax_cafa <- fmax_data[,list(fscore=max(fscore),th=th[which.max(fscore)],pr=pr[which.max(fscore)],rc=rc[which.max(fscore)]),by=list(aspect,tool_name)]
dcast(fmax_cafa,aspect~tool_name,value.var = "th")
dcast(fmax_cafa,aspect~tool_name,value.var = "fscore")

p <- ggplot(plot_data,aes(x=rc,y=pr,color=tool_name))
p <- p + geom_line() + geom_point(data=fmax_cafa,color="#000000",size=3) + geom_point()
p <- p + facet_grid(.~aspect) + labs(x="Recall",y="Precision")
p <- p + scale_color_discrete(guide_legend(title = "CAFA Tool"))
# p <- p + theme(legend.position = "bottom")
p
ggsave("plots/mm-cafa-pr-rc.png",width = 12,height = 3)

p <- ggplot(plot_data,aes(x=th,y=fscore,color=tool_name))
p <- p + geom_line() + geom_point(data=fmax_cafa,color="#000000",size=3) + geom_point()
p <- p + facet_grid(.~aspect) + labs(x="Score",y="F-score")
p <- p + scale_color_discrete(guide_legend(title = "CAFA Tool"))
# p <- p + theme(legend.position = "bottom")
p
ggsave("plots/mm-cafa-fscore.png",width = 12,height = 3)

data[ num_gold==2 & num_pred == 3][order(th,decreasing = T)]

all_annots[db_object_symbol=="GRMZM5G842058" & with >= 0.75 & aspect == "P"]

pred_terms=minimal_set(go_obo,all_annots[db_object_symbol=="GRMZM5G842058" & with >= 0.75 & aspect == "P" & assigned_by=="PANNZER"]$term_accession)
gold_terms=minimal_set(go_obo,all_annots[db_object_symbol=="GRMZM5G842058" & with >= 0.75 & aspect == "P" & assigned_by=="MaizeGDB"]$term_accession)

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