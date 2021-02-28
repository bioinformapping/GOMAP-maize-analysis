library("data.table")
library("ggplot2")
library("reshape2")
library("tools")
library("scales")
library("gridExtra")
library("VennDiagram")
library("DiagrammeR")

aspect2one = function(aspect){
    tmp_list = list(MF="F",BP="P",CC="C")
    return(unlist(tmp_list[aspect]))
}

plot_disc = function(plot_data,tmp_measure,outfile,save_plot=T){
    print(tmp_measure)
    ppt_f = paste("plots/ppt/",outfile,".png",sep="")
    poster_f = paste("plots/",outfile,".png",sep="")
    lims = c(0,max(plot_data[measure==tmp_measure]$value*1.1))

    cbpallete = c("#a6d854","#fdcdac")
    p = ggplot(plot_data[measure==tmp_measure],aes(x=tool,y=value,fill=Type))
    p = p + geom_bar(stat="identity",color="#000000",size=0.5) #+ geom_errorbar()
    p = p + xlab("Dataset") + scale_y_continuous(expand=c(0,0),labels = comma,limits = lims)
    p = p + scale_fill_manual(values=cbpallete) + theme_bw() + ylab(measures[[tmp_measure]])
    p = p + facet_grid(.~aspect,scales = "free_y",labeller = labeller(aspect=aspect_lbl))
    p = p + theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))
    print(p)
    if(save_plot){
        ggsave(ppt_f,width = 10.5,height=6,dpi=300)
        ggsave(poster_f,width = 10.5,height=3,dpi=300)
    }

}

aspect_lbl = function(aspect){
    tmp_lbl = list(C="Cellular Component",P="Biological Process",F="Molecular Function")
    return(tmp_lbl[aspect])
}

theme_plot <- theme(strip.text = element_text(colour = "#000000"),
                    strip.background = element_rect(fill = "#FFFFFF"),
                    axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5),
                    plot.title = element_text(hjust = 0.5))

marginCommon = margin(0, 0.2, 0, 0.2, "cm")
panelSpaceCom = unit(0.1,"cm")
theme_strip_nxaxis <- theme( strip.text = element_text(colour = "#000000"),
                          strip.background = element_rect(fill = "#FFFFFF"),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank(),
                          plot.title = element_blank(),
                          panel.spacing = panelSpaceCom,
                          plot.margin = marginCommon)

theme_strip_nyaxis <- theme( strip.text = element_text(colour = "#000000"),
                             strip.background = element_rect(fill = "#FFFFFF"),
                             axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank(),
                             plot.title = element_blank(),
                             panel.spacing = panelSpaceCom,
                             plot.margin = marginCommon)

theme_nstrip_nxaxis <- theme( strip.text = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank(),
                           plot.title = element_blank(),
                           panel.spacing = panelSpaceCom,
                           plot.margin = marginCommon)

theme_nstrip_xaxis <- theme( strip.text = element_blank(),
                           plot.title = element_blank(),
                           panel.spacing = panelSpaceCom,
                           plot.margin = marginCommon)

annotType_lbl <- function(annotType){
  tmp_list <- list(curated="Curated",predicted="Predicted")
  return(tmp_list[annotType])
}

geneType_lbl <- function(geneType){
  tmp_list <- list(predCount="Prediction & Curation ",curCount="Only Curation")
  return(tmp_list[annotType])
}

geneType_list <- list(predCount="Prediction & Curation ",curCount="Curation Only")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plotMultiVenn <- function(plotFile="test.png",b73Venns=NULL,ph207Venns=NULL,legendCols,mainTitle=NULL){
  png(plotFile,width = 8,height = 5.5,units = "in",res = 300)
  rctGpar = gpar(lty="solid",lwd=1,col="black")
  rectPlace = rectGrob(width = 1, height = 1,gp = rctGpar,default.units = "npc")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 4,widths = c(1,10,10,10),heights = c(1,1,10,10,1))))
  pushViewport(plotViewport(layout.pos.row = 1,
                            layout.pos.col = 2:4))
  grid.text(label=mainTitle,x = 0.5,y = 0.5,gp=gpar( col="black",fontface="bold"))
  popViewport()
  for(i in 1:3){
    pushViewport(plotViewport(layout.pos.row = 2,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.text(label=unlist(aspect_lbl(sort(unique(predData_dt$aspect))[i])),
              x = 0.5,y = 0.5,
              gp=gpar( col="black"))

    popViewport()
    pushViewport(plotViewport(layout.pos.row = 3,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.draw(b73Venns[[i]])
    print(i)
    print(b73Venns[[i]])
    popViewport()
    pushViewport(plotViewport(layout.pos.row = 4,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.draw(ph207Venns[[i]])
    popViewport()
  }
  pushViewport(plotViewport(layout.pos.row = 3,
                            layout.pos.col = 1))
  grid.text(label="B73v4",x = 0.5,y = 0.5,rot = 90,gp=gpar( col="black"))
  popViewport()
  pushViewport(plotViewport(layout.pos.row = 4,
                            layout.pos.col = 1))
  grid.text(label="PH207",x = 0.5,y = 0.5,rot = 90,gp=gpar( col="black"))
  popViewport()
  pushViewport(plotViewport(layout.pos.row = 5,
                            layout.pos.col = 1:4))
  grid.legend(labels = c("GOMAP","Both","Community"),gp=gpar(fill=legendCols,lwd=0),pch=22,nrow=1)
  popViewport()
  dev.off()
}

getVenn = function(data,fillScale=c("#66c2a5","#fc8d62","#8da0cb"),catPos=c(-20,20,180)){
  venn.plot = venn.diagram(data,
               filename = NULL,
               fill=fillScale,
               scaled=F,
               lwd=0,
               margin=0.05,
               cat.pos=catPos,
               ext.pos=180,
               ext.length=0.8,
               cat.cex=1,
               cat.default.pos="outer",
               cat.dist=0.1,
               rotation.degree=0,
               euler.d =F)
  venn.gTree = gTree(children=venn.plot)
  return(venn.gTree)
}

?venn.diagram

plotInbredVenn <- function(plotFile="test.png",geneVenns=NULL,termVenns=NULL,annotVenns=NULL,legendCols,mainTitle=NULL){
  png(plotFile,width = 6,height = 6,units = "in",res = 300)
  rctGpar = gpar(lty="solid",lwd=1,col="black")
  rectPlace = rectGrob(width = 1, height = 1,gp = rctGpar,default.units = "npc")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, 4,widths = c(1,10,10,10),heights = c(1,1,10,10,10))))
  pushViewport(plotViewport(layout.pos.row = 1,
                            layout.pos.col = 2:4))
  grid.text(label=mainTitle,x = 0.5,y = 0.5,gp=gpar( col="black",fontface="bold"))
  popViewport()
  for(i in 1:3){
    pushViewport(plotViewport(layout.pos.row = 2,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.text(label=unlist(aspect_lbl(sort(unique(predData_dt$aspect))[i])),
              x = 0.5,y = 0.5,
              gp=gpar( col="black"))

    popViewport()
    pushViewport(plotViewport(layout.pos.row = 3,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.draw(geneVenns[[i]])
    print(i)
    popViewport()
    pushViewport(plotViewport(layout.pos.row = 4,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.draw(termVenns[[i]])
    popViewport()
    pushViewport(plotViewport(layout.pos.row = 5,
                              layout.pos.col = i+1))
    grid.draw(rectPlace)
    grid.draw(annotVenns[[i]])
    popViewport()
  }
  pushViewport(plotViewport(layout.pos.row = 3,
                            layout.pos.col = 1))
  grid.text(label="Genes",x = 0.5,y = 0.5,rot = 90,gp=gpar( col="black"))
  popViewport()
  pushViewport(plotViewport(layout.pos.row = 4,
                            layout.pos.col = 1))
  grid.text(label="GO Terms",x = 0.5,y = 0.5,rot = 90,gp=gpar( col="black"))
  popViewport()
  pushViewport(plotViewport(layout.pos.row = 5,
                            layout.pos.col = 1))
  grid.text(label="Annotations",x = 0.5,y = 0.5,rot = 90,gp=gpar( col="black"))
  popViewport()
  # pushViewport(plotViewport(layout.pos.row = 6,
  #                           layout.pos.col = 1:4))
  # grid.legend(labels = c("GOMAP","Both","Community"),gp=gpar(fill=legendCols,lwd=0),pch=22,nrow=1)
  # popViewport()
  dev.off()
}
