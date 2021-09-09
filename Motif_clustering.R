
# load motifs
load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
# correct names
TF_clusters_PWMs$metadata$motif_description2 <- as.character(TF_clusters_PWMs$metadata$motif_description)
TF_clusters_PWMs$metadata$motif_description2 <- sapply(strsplit(TF_clusters_PWMs$metadata$motif_description2,"_"), `[`, 1)
TF_clusters_PWMs$metadata$motif_description2 <- sapply(strsplit(TF_clusters_PWMs$metadata$motif_description2,"/"), `[`, 1)
TF_clusters_PWMs$metadata$motif_description2 <- sapply(strsplit(TF_clusters_PWMs$metadata$motif_description2,"\\("), `[`, 1)
TF_clusters_PWMs$metadata$motif_description2 <- sapply(strsplit(TF_clusters_PWMs$metadata$motif_description2,"\\["), `[`, 1)
TF_clusters_PWMs$metadata$motif_description2[complete.cases(TF_clusters_PWMs$metadata$Dmel)] <- TF_clusters_PWMs$metadata$Dmel[complete.cases(TF_clusters_PWMs$metadata$Dmel)]

TF_clusters_PWMs$metadata$motif_description2 <- gsub("ttk-PF", "ttk", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("tramtrack", "ttk", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("Kruppel", "Kr", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("activating-protein 1 \\(FOS-JUN heterodimer\\)", "kay_Jra", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("homer__GCTGATAASV_Unknown5", "GATAA", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("cisbp__M2335", "GATAA", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("TIFDMEM0000091", "Ohler6", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("ENSG00000235187", "Ets", TF_clusters_PWMs$metadata$motif_description2)
TF_clusters_PWMs$metadata$motif_description2 <- gsub("Ohler5\\(E-box\\)", "Ohler5", TF_clusters_PWMs$metadata$motif_description2)

#####
# Step 3: Hierarchically cluster motifs by similarity
#####
# based on Jeff Vierstra https://raw.githubusercontent.com/jvierstra/motif-clustering/master/hierarchical.py

sim_file <- "tomtom.all.treated.txt"
sim <- read.delim(sim_file)
simsq <- data.table::dcast(sim, Query_ID~Target_ID, value.var = "E.value")
rownames(simsq) <- simsq$Query_ID
simsq <- simsq[,-1]
simsq[is.na(simsq)] <- 100
simsq[1:20,1:20]

mat = -log10(simsq)
mat[mat == Inf]=10 # the ones that had 0 euclidean distance in the beginning

# We then performed hierarchical clustering using Pearson correlation as the distance metric and complete linkage
# comparison between columns
tmp_cor <- cor(mat, method="pearson")
Z = hclust(as.dist(1-tmp_cor), method = 'complete')

### test clustering at a range of tree heights (0.5-1)
step=0.1
start=0.5
end=1.0

thresholds=seq(start, end, step)

pdf("Hierarchical_clusters_diff_thresholds.pdf", width = 20, height = 5)
for(thresh in thresholds){
  
  cl = dendextend:::cutree(Z, h=thresh, order_clusters_as_data = FALSE)
  df = data.frame(Motifs=names(mat),
                  Cluster=cl[match(names(mat), names(cl))])
  write.table(df, paste0("clusters.", thresh,".txt"), sep="\t", row.names = F, quote=F)

  plot(Z, labels=FALSE, main=paste0("tree height: ",thresh, " - ", length(unique(df$Cluster)), " clusters"))
  abline(h=thresh, col="red")
  
  print(paste0("tree height: ",thresh, " - ", length(unique(df$Cluster)), " clusters"))
  
}
dev.off()


### choose final clusters cutting the dendrogram at height 0.8
thresh=0.8
cl = dendextend:::cutree(Z, h=thresh, order_clusters_as_data = FALSE)
df = data.frame(Motifs=names(mat),
                Cluster=cl[match(names(mat), names(cl))],
                Order_dendogram=match(names(mat), Z$labels[Z$order]))
df <- merge(df, TF_clusters_PWMs$metadata[,c(1,13,10)], by=1)
df <- df[order(df$Order_dendogram),]
length(unique(df$Cluster))
sort(table(df$Cluster))
save(Z, mat, file = paste0("All_motifs_data_and_hclust_objects.Rdata"))
saveRDS(df, paste0("All_motifs_final_clusters_thresh", thresh, ".rds"))


### plot hierarchical clustering heatmap of motifs clustered by simililarity and clusters identified cutting the dendrogram at height 0.8
# Notice that image interprets the z matrix as a table of f(x[i], y[j]) values, so that the x axis corresponds to row number and the y axis to column number,
# with column 1 at the bottom, i.e. a 90 degree counter-clockwise rotation of the conventional printed layout of a matrix.
# that's why I need to reverse the order of the columns
# top-right should represent cluster1. - and so on
png(paste0("All_motifs_hierarchically_clustered_heatmap_pairwise_similarity_scores.png"), type="cairo", width = 2000, height = 2000, res = 300)
out <- tmp_cor[Z$order,rev(Z$order)]
image(out, col=c("white", "black"), # colorRampPalette(c("grey100", "grey0"))(100)
      las=1, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
dev.off()

# highlight specific clusters on the hierarchical clustering heatmap
for(c in 1:length(unique(df$Cluster))){
  png(paste0("Highlight_cluster_",c,".png"), type="cairo", width = 2000, height = 2000, res = 50)
  highli <- names(cl)[which(cl==c)]
  image(out, col=c("white", "black"),
        las=1, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  axis(1, at=seq(0,1,length.out = nrow(out))[rownames(out) %in% highli], labels=rep("",length(highli)), col = "red")
  axis(2, at=seq(0,1,length.out = nrow(out))[colnames(out) %in% highli], labels=rep("",length(highli)), col = "red")
  dev.off()
  print(c)
}



#####
# Step 4: Annotation of motif clusters: Clusters were manually curated and annotated with the respective motif types
#####

### annotation of clusters with names

# if cluster has 1-5 elements -> use most common (or random) motif name
Cluster_IDs <- unique(df$Cluster)
names(Cluster_IDs) <- Cluster_IDs
Cluster_IDs <- sapply(Cluster_IDs, function(c){
  tmp <- df[df$Cluster == c,]
  if(nrow(tmp)<6){
    # first get the most common drosophila name
    if(length(which(complete.cases(tmp$Dmel)))>0) out <- names(sort(table(tmp$Dmel), decreasing = TRUE)[1])
    if(length(which(complete.cases(tmp$Dmel)))==0) out <- names(sort(table(tmp$motif_description2), decreasing = TRUE)[1])
    return(out)
  }else(return(NA))
})

df$Cluster_name <- Cluster_IDs[match(df$Cluster, names(Cluster_IDs))]
# for the remaining clusters do it through manual curation

# save table
out <- df[order(df$motif_group),]
write.table(out[!duplicated(out$Cluster),c(2,6,7)], paste0("All_final_clusters_thresh", thresh, "_annotated.txt"), sep="\t", row.names = F, quote=F)



#####
# Step 5: Plot motif logos aligned for each motif cluster
#####
# For each of the clusters, we then selected a seed motif model (top with absolute enrichment in dev or hk) to which we aligned all other motifs within cluster (both position and orientation; order motifs by absolute enrichment in dev or hk).

# load motif cluster names
TF_motif_clusters_manual_annotation <- read.csv("TF_motif_clusters_threshold0.8_manual_annotation.csv")

### done in R script Plot_motif_logos.R
# add name to PDF filename

# plot clusters
library(motifStack)

# prepare motifs in motifStack format
PWM_candidates <- TF_clusters_PWMs$metadata[TF_clusters_PWMs$metadata$X..motif_collection_name %in% c("bergman",
                                                                                                      "cisbp",
                                                                                                      "flyfactorsurvey",
                                                                                                      "homer",
                                                                                                      "jaspar",
                                                                                                      "stark",
                                                                                                      "idmmpmm"),]
TF_clusters_PWMs$All_pwms_perc_selected <- TF_clusters_PWMs$All_pwms_perc[name(TF_clusters_PWMs$All_pwms_perc) %in% PWM_candidates$motif_name]
All_motifs <- lapply(1:length(TF_clusters_PWMs$All_pwms_perc_selected), function(x){
  new("pfm", mat=TF_clusters_PWMs$All_pwms_perc_selected[[x]]@profileMatrix,
      name=paste0(name(TF_clusters_PWMs$All_pwms_perc_selected)[x], " (",
                  TF_clusters_PWMs$metadata$motif_description2[match(name(TF_clusters_PWMs$All_pwms_perc_selected)[x], TF_clusters_PWMs$metadata$motif_name)], ") "))
})
names(All_motifs) <- name(TF_clusters_PWMs$All_pwms_perc_selected)

saveRDS(All_motifs, "All_motifs_PWMs_motifStack_format.rds")
# All_motifs <- readRDS("All_motifs_PWMs_motifStack_format.rds")

for(c in 1:max(df$Cluster)){
  # c=4
  tmp <- df[df$Cluster == c,]
  PWM_tmp <- All_motifs[match(tmp$Motifs, names(All_motifs))]
  out <- paste0("Clusters_logos/Cluster", c, "_",
                TF_motif_clusters_manual_annotation$Cluster_name[TF_motif_clusters_manual_annotation$Cluster==c], "_",
                nrow(tmp), "motifs.pdf")
  if(nrow(tmp)>2){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/3
    h=length(PWM_tmp)
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp, layout = "tree", xaxis=F,yaxis=F, xlcex=0, ylcex=0, ncex=0.9)
    
  }else if(nrow(tmp)>=200){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/3
    h=length(PWM_tmp)
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp, xaxis=F,yaxis=F, xlcex=0, ylcex=0, ncex=0.9)
    
  }else if(nrow(tmp)==2){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/3
    h=length(PWM_tmp)*2
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp, xaxis=F,yaxis=F, xlcex=0, ylcex=0)
    
  }else if(nrow(tmp)==1){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/2
    h=3
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp[[1]], xaxis=F,yaxis=F, xlcex=0, ylcex=0)
    
  }
  dev.off()
  print(paste0("Cluster ", c))
}


#####
# Step 6: Curate metadata information with cluster information and save PWM models into single R object (markdown Create_consensus_TF_motif_database.Rmd)
#####
