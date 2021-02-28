library("data.table")
library("ggplot2")
library("naturalsort")

walltime <- fread("data/walltime/walltime.csv")
wall_melt <- melt.data.table(walltime,"Step",c("Mo17","W22","PH207","B73v4"))

plotData <- wall_melt[,c("h","m"):=tstrsplit(wall_melt$value,"[hm]")]
plotData[,h:=as.numeric(h)]
plotData[,m:=as.numeric(m)]
plotData[,totMin:=h*60+m]

stepLvls <- c('seqsim','domain','mixmeth-blast','fanngo','mixmeth-preproc','mixmeth','aggregate')
plotData[,Step:=factor(Step,levels = stepLvls)]
plotData[,variable:=naturalfactor(variable)]
initSteps = c("seqsim","domain","mixmeth-blast","fanngo")
plotData[Step %in% initSteps,ymin:=0]

plotData[Step=="mixmeth-preproc"]$ymin = plotData[Step=="mixmeth-blast"]$totMin
plotData[Step=="mixmeth"]$ymin = plotData[Step=="mixmeth-blast"]$totMin + plotData[Step=="mixmeth-preproc"]$totMin
plotData[Step=="aggregate"]$ymin = plotData[Step=="mixmeth-blast"]$totMin + plotData[Step=="mixmeth-preproc"]$totMin + plotData[Step=="mixmeth"]$totMin
plotData[,ymax:=ymin+totMin]

cbpallet <- c("#b3de69","#fdb462","#80b1d3","#fb8072","#bebada","#f2e6d9","#8dd3c7")


p <-  ggplot(plotData,aes(x=variable,ymin=ymin/60,ymax=ymax/60,color=Step)) +
      geom_linerange(stat = "identity",position = position_dodge2(width = 1),size=4) +
      geom_hline(mapping = aes(yintercept=24),linetype = 2) + 
      geom_label(mapping = aes(x=3.1,y=24,label="24h"),color="black") + 
      scale_color_manual(values = cbpallet) +
      scale_y_continuous(limits = c(0,24)) +
      theme_bw(base_size = 18) +
      coord_flip() +  
      labs(x="Maize Genomes",y="Time (Hours)",fill="GOMAP Step")
p
ggsave("figures/walltime.png",width = 8,height = 4,units = "in",dpi = 300)

