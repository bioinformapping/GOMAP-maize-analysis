source("R/plotHelper.R")
source("R/plotTools.R")
library("data.table")
library("ggplot2")

#library("scales")
library("kableExtra")

sraData = fread("data/sra/sra_query.tsv")
sraData[,year:=format(as.Date(sraData$pub_date),format = "%Y")]
sraSummary <- sraData[!is.na(year) & year < 2020,list(studies=length(unique(study)),runs=.N),by=year][order(year)]
sraSummary
sraSummaryMelt <- melt.data.table(sraSummary,id.vars = "year")
sraSummaryMelt[,variable:=factor(variable,levels = c("studies","runs"),labels = c("Studies","Runs"))]
sraSummaryMelt


sraPlot <-  ggplot(sraSummaryMelt,mapping = aes(x=year,y=value,fill=year),color="black") + 
            geom_bar(stat = "identity",color="black") + 
            facet_grid(variable~.,scales = "free") + 
            scale_fill_brewer(type = "seq",palette = 15) +  theme_bw() + theme_plot +
            theme(axis.text.x = element_text(hjust = 0.5,angle = 0)) + 
            labs(fill="Year",x="Year",y="Count") + guides(fill=F)
sraPlot

ggsave("figures/sraDatasets.png",width = 5.5,height = 3)

colSums(sraSummary[,.(studies,runs)])
colnames(sraSummary) <- c("Year","Number of Studies","Num of Datasets")
cat(kable(sraSummary,format = "markdown",caption = "SRA Plant Datasets"),file = "tables/README.md",sep="\n", append=TRUE)

