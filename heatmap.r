

library(pheatmap)
library(preprocessCore)
library(RColorBrewer)
library(qdapRegex)
library(optparse)
library(data.table)
library(splitstackshape)

args = commandArgs(trailingOnly=TRUE)

filename <- args[1]
print (filename)
tmp <- strsplit(filename,"_")[[1]][c(1,4)]
mark <- tmp[1]
tmp1 <- paste0(tmp[1],"_",tmp[2])
sample <- strsplit(tmp1,".",fixed=T)[[1]][1]
which.type <- strsplit(sample,"_")[[1]][2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_most_variable_peaks <- function (data,n){
  # estimate the variance for each row in the matrix
  var_peaks <- apply(data, 1, var)
  head(var_peaks)
  
  # Get the peaks for the top 10000 most variable genes
  select_var <- (order(var_peaks, decreasing=TRUE))[1:n]
  head(select_var)
  
  # Subset matrix
  var_data <- data[select_var,]
  dim(var_data)
  return (var_data)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kmeans_clust <- function (data,k){
  #set.seed(1000)
  k.data <- kmeans(data,centers = k)
  data$id <- row.names(data)
  cluster.assignments <- data.frame(k.data$cluster)
  colnames(cluster.assignments) <- c("cluster.number")
  cluster.assignments$id <- row.names(cluster.assignments)
  merged <- merge(data,cluster.assignments,by="id")
  # print (merged[merged$cluster.number==1,])
  merged.sorted <- merged[order(as.integer(as.character(merged$cluster.number))),]
  row.names(merged.sorted) <- merged.sorted$id
  merged.sorted <- merged.sorted[,c(2:(ncol(merged.sorted)-1))]
  return (merged.sorted)
  #print (head(merged.sorted))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot_heatmap <- function (data, mark, title, n=NA){
  breaksList = seq(-1,1,by=0.10)
  cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(length(breaksList)-1)
  
  annotation_col <- data.frame(patient=colnames(data), stringsAsFactors=F)
  annotation_col <- data.frame(cSplit(annotation_col,"patient","-"),stringsAsFactors = F)
  names(annotation_col) <- c("Patient","Status")

  if (length(unique(annotation_col$Status))==2){
    rownames(annotation_col) = colnames(data)
    #annotation_col$patient <- gsub("\\-.*","",annotation_col$patient)
    
    library(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector_p <- col_vector[-c(14,16)]
    n.color <- sample(col_vector_p,(length(names(data))/2))
    ann_colors = list(Patient = c(n.color),Status=c(col_vector[c(14,16)]))
    
    names(ann_colors$Patient) = unique(annotation_col$Patient)
    names(ann_colors$Status)= unique(annotation_col$Status)
    
    ann_colors
    
    pheatmap(data, color = cell_colors, scale='row', border_color = NA, 
             cluster_rows = F, cluster_cols = T,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean", 
             clustering_method = "complete",
             fontsize = 5,
             fontsize_col = 6,
             #cellwidth = 6.5,
             breaks = breaksList,
             kmeans_k = NA,
             show_colnames = T,
             show_rownames = F,
             main=paste0(mark,":",title),annotation_col = annotation_col,
             annotation_colors = ann_colors)
  } else if(length(unique(annotation_col$Status))==1) {
    pheatmap(data, color = cell_colors, scale='row', border_color = NA, 
             cluster_rows = F, cluster_cols = T,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean", 
             clustering_method = "complete",
             fontsize = 7,
             fontsize_col = 7,
             #cellwidth = 7,
             breaks = breaksList,
             kmeans_k = NA,
             show_colnames = T,
             show_rownames = F,
             main=paste0(mark,":",title))
  }
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

generate_plot <- function(data,
                          get.top.peaks=TRUE,
                          perform.quantile.normalization=FALSE,
                          perform.kmeans.clustering=TRUE,
                          mark, type.values , filename) {
  
  if (perform.quantile.normalization==TRUE) { 
    c.names <- names(data)
    r.names <- row.names(data)
    data <- as.data.frame(normalize.quantiles(as.matrix(data)),stringsAsFactors = F)
    names(data) <- c.names
    row.names(data) <- r.names
  }
  if (get.top.peaks==TRUE) { data.var <- get_most_variable_peaks(data,10000) }
  if (perform.kmeans.clustering==TRUE) { data.var <-  kmeans_clust(data.var,5) }
  #png(filename,height = 1000, width=1000, res = 100)
  plot_heatmap(data.var,mark,type.values,n=NA)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#diff.raw <- read.table("~/Downloads/peaks_matrix.tsv",sep="\t",header=T)
diff.raw <- read.table(filename,sep="\t",header=T)
diff.raw$X <- as.character(diff.raw$X)
diff.raw <- diff.raw[grep("chrX",diff.raw$X,invert = T),]
diff.raw <- diff.raw[grep("chrY",diff.raw$X,invert = T),]
if (grep("D-H3K9AC_peaks_matrix",filename)==1){diff.raw <- (diff.raw[,grep(".R.H3K9AC.",names(diff.raw),invert = T)])}

row.names(diff.raw) <- diff.raw$X
diff <- diff.raw[,-1]

diff <- (diff[,order(colnames(diff))])
names(diff) <- unlist(rm_between(names(diff), "bam.", ".ali", extract=TRUE))
names(diff)

drop.cols <- grep("SAH", colnames(diff))
if (length(drop.cols)>0){
  print ("SAH")
  diff <- diff[, -(drop.cols)]
}
drop.cols <- NULL
drop.cols <- grep("GHW", colnames(diff))
if (length(drop.cols)>0){
  print ("GHW")
  diff <- diff[, -(drop.cols)]
}
drop.cols <- NULL
drop.cols <- grep("PAVIFZ", colnames(diff))
if (length(drop.cols)>0){
  print ("PAVIFZ")
  diff <- diff[, -(drop.cols)]
}

for (i in 1:length(names(diff))){
  tmp.name <- (strsplit(names(diff)[i],".",fixed=T)[[1]][1:2])
  new.name <- paste0(tmp.name[1],"-",tmp.name[2])
  names(diff)[i] <- new.name
}
names(diff)

outfile <- paste0(sample,".pdf")

#plot.new()
pdf(outfile,onefile=F)
generate_plot(diff,
              perform.quantile.normalization=FALSE,
              get.top.peaks=TRUE,
              perform.kmeans.clustering=TRUE,
              mark=mark, 
              type.values=which.type , 
              filename=NA)

dev.off()


outfile <- paste0(sample,"_quantile.norm",".pdf")
pdf(outfile,onefile=F)
generate_plot(diff,
              perform.quantile.normalization=TRUE,
              get.top.peaks=TRUE,
              perform.kmeans.clustering=TRUE,
              mark=mark, 
              type.values=paste0(which.type,".quantile.norm") , 
              filename=NA)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


