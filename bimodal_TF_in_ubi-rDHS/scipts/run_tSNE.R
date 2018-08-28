
# -- Kaili
# This script is for runing t-SNE and get all the figures.
# EXP: Rscript run_tSNE.R /data/public_html_users/fankaili/tSNE_10921/ /data/zusers/fankaili/ccre/tf/tsne/matrix/ /data/zusers/fankaili/ccre/hg19_ubi_ccRE_ccreid_agnostic.txt /data/zusers/fankaili/ccre/tf/tsne/hg19_ubi-rDHS_DNAme_matrix.txt GM12878


library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(grid)

args<-commandArgs(T)
outDir = args[1]
matrixDir = args[2]
agnosticFile = args[3]
dnameFile = args[4]
celltype0 = args[5]

if(startsWith(celltype0, "2_")){
  celltype = substring(celltype0,3)
}else{
  celltype = celltype0
}

##### local path for test
# outDir="/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/tSNE/tSNE_10921/"
# matrixDir="/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/tSNE/matrix/"
# agnosticFile="/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/tSNE/hg19_ubi_ccRE_ccreid_agnostic.txt"
# dnameFile="/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/tSNE/hg19_ubi-rDHS_DNAme_matrix.txt"
# celltype="GM12878"
# cellList = as.vector(read.table("/Users/kaili/Dropbox (UMass Medical School)/Project/ccre/bimodal_tf_rdhs/data_figure/tSNE/hg19_cellline_four_dimension_list.txt")[,1])

###################

setwd(outDir)
data = read.table(paste(matrixDir,"hg19_ubi-rDHS_",celltype0,"_five_dimension_matrix.txt", sep=""), header = TRUE, row.names = 1)
agnostic = read.table(agnosticFile)
dname = read.table(dnameFile, header = TRUE, row.names = 1)

cellList = read.table("/data/zusers/fankaili/ccre/tf/hg19_cellline_four_dimension_list.txt")

################### main
set.seed(1)
tSNE <- Rtsne(data[,1:5], dim = 2, perplexity=30)
out = data.frame(tSNE$Y)

if(celltype=="GM12878"){
  write.table(out, paste("/data/zusers/fankaili/ccre/tf/tsne/tsne_table/hg19_ubi-rDHS_", celltype0,"_tNSE.txt", sep=""), 
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}


for (i in 0:nrow(cellList)){
  if(i==0){
    type = celltype
  }else if(as.character(cellList[i,1])==celltype){
    type = "NA"
  }else{
    type = as.character(cellList[i,1])
  }
  #
  if(type!="NA"){
    if(startsWith(celltype0, "2_")){
      data2=read.table(paste(matrixDir,"hg19_ubi-rDHS_2_",type,"_five_dimension_matrix.txt", sep=""), header = TRUE, row.names = 1)
    }else{
      data2=read.table(paste(matrixDir,"hg19_ubi-rDHS_",type,"_five_dimension_matrix.txt", sep=""), header = TRUE, row.names = 1)
    }
    pdfName = paste("./",celltype,"/hg19_ubi-rDHS_", celltype, "_tSNE_in_",type,".pdf", sep="")
    pdf(pdfName, height = 12, width = 18)
    ### agnostic
    a = ggplot(out, aes(x = X1, y = X2, color = agnostic$V3)) + 
      geom_point(alpha = 0.3) + theme_gray() + 
      ggtitle("agnostic") + 
      scale_color_manual(limits=c("Promoter-like_Signatures","Enhancer-like_Signatures","CTCF-only"),
                         values=c("#ff0000","#ffcd00","#00b0f0"), name = 'agnostic',
                         breaks=c("Promoter-like_Signatures","Enhancer-like_Signatures","CTCF-only"),
                         labels=c("PLS (n=9009)","ELS (n=1826)","CTCF (n=86)"))
    ### cell-type specific
    aa = unlist(summary(data2$type))
    bb = c(paste("PLS"," (n=",aa["Promoter-like_Signatures"], ")", sep=""),
           paste("ELS"," (n=",aa["Enhancer-like_Signatures"], ")", sep=""),
           paste("CTCF"," (n=",aa["CTCF-only"], ")", sep=""),
           paste("DNase"," (n=",aa["DNase-only"], ")", sep=""),
           paste("Inactive"," (n=",aa["Inactive"], ")", sep=""))
    b = ggplot(out, aes(x = X1, y = X2, color = data2$type)) + 
      geom_point(alpha = 0.3) + theme_gray() + ggtitle(type) + 
      scale_color_manual(limits=c("Promoter-like_Signatures","Enhancer-like_Signatures","CTCF-only","DNase-only","Inactive"),
                         values=c("#ff0000","#ffcd00","#00b0f0","#06da93","#e1e1e1"),
                         name = paste(substring(type,1,20),"\nspecific five-group", sep=""), 
                         breaks=c("CTCF-only","DNase-only","Enhancer-like_Signatures","Inactive","Promoter-like_Signatures"),
                         labels=bb)
    ### distance
    c = ggplot(out, aes(x = X1, y = X2, color = data2$distance)) + 
      geom_point(alpha = 0.3) + theme_gray() + ggtitle("distance2TSS") + 
      scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "distance2TSS\n(log10)")
    ### DNase
    d = ggplot(out, aes(x = X1, y = X2, color = data2$dnase)) + 
      geom_point(alpha = 0.3) + theme_gray() + ggtitle("DNase")  + 
      scale_color_gradient2(low = "#00b0f0", mid= "#f7f7f7", high = "#ff0000", midpoint = 0, name = "DNase\n(z-score)")
    # H3K4me3
    e = ggplot(out, aes(x = X1, y = X2, color = data2$h3k4me3)) + 
      geom_point(alpha = 0.3) + theme_gray() + ggtitle("H3K4me3")  + 
      scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), 
                            values = rescale(c(-10,-5,0,5)), name = "H3K4me3\n(z-score)")
    # H3K27ac
    f = ggplot(out, aes(x = X1, y = X2, color = data2$h3k27ac)) + 
      geom_point(alpha = 0.3) + theme_gray() + ggtitle("H3K27ac")  + 
      scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), 
                            values = rescale(c(-10,-5,0,5)), name = "H3K27ac\n(z-score)")
    # CTCF
    g = ggplot(out, aes(x = X1, y = X2, color = data2$ctcf)) + 
      geom_point(alpha = 0.3) + theme_gray() + ggtitle("CTCF")  + 
      scale_color_gradientn(colours = c("#00b0f0","#00b0f0","#f7f7f7","#ff0000"), 
                            values = rescale(c(-10,-5,0,5)), name = "CTCF\n(z-score)")
    # DNAme
    if(is.element(type,colnames(dname))){
      h = ggplot(out, aes(x = X1, y = X2, color = dname[,type])) + 
        geom_point(alpha = 0.3) + theme_gray() + ggtitle("DNAme") +
        scale_color_gradientn(colours = c("#00b0f0","#f7f7f7","#ff0000"), 
                              values = rescale(c(0,0.5,1)), name = "DNAme",
                              limits = c(0,1))
      title = paste(celltype, "(ubi-rDHS, n=",nrow(data),"), color by ", type, sep="")
      grid.arrange(a, b, c, d, e, f, g, h, ncol=4, newpage = TRUE,
                   top=textGrob(title,gp=gpar(fontsize=20,font=3)))
    }
    else{
      title = paste(celltype, "(ubi-rDHS, n=",nrow(data),"), color by ", type, sep="")
      grid.arrange(a, b, c, d, e, f, g, nrow=2, ncol=4, newpage = TRUE,
                   top=textGrob(title,gp=gpar(fontsize=20,font=3)))
    }
    dev.off()
  }
}


