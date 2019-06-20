
library(stringr)
library(dplyr)
library(readr)
library(plyr)
library(igraph)
library(GenomicRanges)
library(tidyr)
library(data.table)

big_function_RNA_seq_ass <- function(main_data, cell_type) {
  main_data$gene_id <- gsub("\\..*", "", main_data$gene_id)


  ID <- grepl("^ID\\d", colnames(annotation))
  anno_id <- colnames(annotation[, ID])
  mean_value <- data_frame()
  for (id in anno_id) {
    anno <- merge(main_data, annotation, by.y = id, by.x = "gene_id")
    mean_value <- bind_rows(mean_value, anno)
  }

  mean_value <- select(mean_value, -c(gene_id, paste("ID", seq(1, 6), sep = ""), region)) %>%
    group_by(ID, gene_name, Chr, Start, End) %>%
    summarise_all("mean")
  mean_value
  write_csv(mean_value, path = paste(cell_type, "mean_value_rnaseq_with_chromosome_ID.csv"))

  bedepi <- with(mean_value, GRanges(Chr, ranges = IRanges(Start, End), strand = Rle(strand(rep("*", nrow(mean_value))))))
  mcols(bedepi) <- mean_value[, -c(1:5)]
  names(bedepi) <- mean_value$ID
  PCHiC_map <- read_tsv("/home/ngoc/Desktop/Ngoc_CRCT/1st week/data/monocyte/PCHiCJavierreCell2016_merged_samples_12Apr2015_full.txt", skip = 5, col_types = cols(oeChr = col_character(), baitChr = col_character()))
  PCHiC_wt <- PCHiC_map[, c(12:28)]
  PCHiC_wt <- PCHiC_map[any(PCHiC_map[, c(12:28)] > 5), ]
  PCHiC_wt$baits <- paste(PCHiC_wt$baitChr, PCHiC_wt$baitStart, PCHiC_wt$baitEnd, sep = "_")
  PCHiC_wt$OEs <- paste(PCHiC_wt$oeChr, PCHiC_wt$oeStart, PCHiC_wt$oeEnd, sep = "_")
  PCHiC_baits <- unique(PCHiC_wt$baits)
  PCHiC_OE_all <- unique(PCHiC_wt$OEs)

  PCHiC_PromOE <- PCHiC_OE_all[which(PCHiC_OE_all %in% PCHiC_baits)]
  PCHiC_nonPromOE <- PCHiC_OE_all[-which(PCHiC_OE_all %in% PCHiC_baits)]
  PCHiC_wholenet <- graph_from_data_frame(PCHiC_wt[, c(29:30)], directed = FALSE)
  PCHiC_PP <- delete.vertices(PCHiC_wholenet, V(PCHiC_wholenet)[which(V(PCHiC_wholenet)$name %in% PCHiC_nonPromOE)])
  PCHiC_POE <- delete.vertices(PCHiC_wholenet, V(PCHiC_wholenet)[which(V(PCHiC_wholenet)$name %in% PCHiC_PromOE)])

  colnames(PCHiC_wt)[c(1:3, 6:8)] <- rep(c("chr", "start", "end"), 2)
  PCHiC_bed <- unique(rbind(PCHiC_wt[, c(1:3)], PCHiC_wt[, c(6:8)]))
  PCHiC_GRange <- with(PCHiC_bed, GRanges(chr, IRanges(start, end)))
  PCHiC_GRange$ID <- paste(PCHiC_bed$chr, PCHiC_bed$start, PCHiC_bed$end, sep = "_")

  overlaps <- findOverlaps(PCHiC_GRange, bedepi)
  match_hit <- data.frame(PCHiC_GRange$ID[queryHits(overlaps)], as.data.frame(mcols(bedepi)[subjectHits(overlaps), ]), stringsAsFactors = T)
  colnames(match_hit)[1] <- "fragment"

  data.dt <- data.table(match_hit)
  setkey(data.dt, fragment) # sorts ascending by fragment
  agchic <- data.frame(data.dt[, lapply(.SD, mean), by = fragment]) # mean of ChIP-seq features by fragment
  agchic <- data.frame(agchic[, -1], row.names = agchic[, 1])

  calc_assort <- function(G, data) {
    names <- colnames(data)
    ass <- list()
    G_epi <- list()
    for (i in c(1:ncol(data))) {
      G <- set.vertex.attribute(G, names[i], value = data[V(G)$name, names[i]])
      attsel <- which(names(vertex.attributes(G)) == names[i])
      G_epi[[names[i]]] <- delete.vertices(G, V(G)[is.na(vertex.attributes(G)[[attsel]])])
      ass[[names[i]]] <- assortativity(G_epi[[names[[i]]]], types1 = vertex.attributes(G_epi[[names[i]]])[[attsel]], directed = F)
    }

    return(ass)
  }

  ab_PCHiC_all <- colMeans(agchic[which(rownames(agchic) %in% V(PCHiC_wholenet)$name), ])
  ab_PCHiC_PP <- colMeans(agchic[which(rownames(agchic) %in% V(PCHiC_PP)$name), ])
  ab_PCHiC_POE <- colMeans(agchic[which(rownames(agchic) %in% V(PCHiC_POE)$name), ])

  ass_PCHiC_all <- calc_assort(PCHiC_wholenet, agchic)
  ass_PCHiC_PP <- calc_assort(PCHiC_PP, agchic)
  ass_PCHiC_POE <- calc_assort(PCHiC_POE, agchic)

  namesvec <- as.data.frame(colnames(agchic))
  namesvec$names <- gsub("(\\.).*", "", namesvec[, 1])
  names <- unique(namesvec$names)
  annotation_feature <- read_tsv("/home/ngoc/Desktop/Ngoc_CRCT/1st week/data/monocyte/selected dataset_Blueprint_RNAseq.tsv")[, 1:3]
  annotation_feature <- annotation_feature[which(annotation_feature$ID %in% names), ]
  annotation_feature$Description <- c("CD14+,CD16-", "none", "T=0day", "monocyte_none")
  annotation_feature
  namesvec <- merge(namesvec, annotation_feature, by.x = "names", by.y = "ID")
  colnames(namesvec)[2] <- "ID"
  namesvec
  data_anno <- cbind(ass_PCHiC_all, ab_PCHiC_all, ass_PCHiC_PP, ab_PCHiC_PP, ass_PCHiC_POE, ab_PCHiC_POE)

  data_anno <- merge(data_anno, namesvec, by.x = "row.names", by.y = "ID")
  nodecats <- namesvec[, c(2, 4)]
  colnames(nodecats)[1] <- "ID"
  namesvec <- colnames(agchic)
  nodecats
  cols10 <- rainbow(4)
  names(cols10) <- annotation_feature$ID
  cats <- as.vector(unique(nodecats[, 2]))
  cols10 <- rainbow(11)
  names(cols10) <- cats
  cols10["Other"] <- "grey"
  cols4plot <- rep("grey", 19)
  names(cols4plot) <- namesvec
  for (c in cats) {
    rel <- which(nodecats[, 2] == c)
    cols4plot[rownames(nodecats)[rel]] <- cols10[c]
  }

  rownames(data_anno) <- data_anno$Row.names

  rownames(nodecats) <- nodecats[, 1]
  cats <- as.vector(unique(nodecats[, 2]))
  cols10 <- rainbow(4)
  names(cols10) <- cats


  cols10["Other"] <- "grey"
  cols4plot <- rep("grey", 19)
  names(cols4plot) <- namesvec
  for (c in cats) {
    rel <- which(nodecats[, 2] == c)
    cols4plot[rownames(nodecats)[rel]] <- cols10[c]
  }

  pdf(paste(cell_type, "PCHiC all assortativity vs abundance.pdf"))
  labs <- c(1:ncol(agchic))
  labs2 <- rep("", ncol(agchic))
  sel <- which(ass_PCHiC_all > 0.05 | ass_PCHiC_all < 0 | ab_PCHiC_all > 0.03)
  labs2[sel] <- namesvec[sel]
  plot(ab_PCHiC_all, unlist(ass_PCHiC_all), pch = 20, col = cols4plot, main = "PCHi-C network", xlab = "Abundance", ylab = "ChAs")
  text((ab_PCHiC_all), jitter(unlist(ass_PCHiC_all), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
  abline(h = 0)
  legend("bottomright", legend = as.vector(unique(nodecats[, 2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[, 2]))], cex = 0.7)
  dev.off()

  pdf(paste(cell_type, "PCHiC Promoter-Promoter assortativity vs abundance.pdf"))
  labs <- c(1:ncol(agchic))
  labs2 <- rep("", ncol(agchic))
  sel <- which(ass_PCHiC_PP > 0.001 | ass_PCHiC_PP < 0)
  labs2[sel] <- names[sel]
  plot(ab_PCHiC_PP, unlist(ass_PCHiC_PP), pch = 20, col = cols4plot, main = "PCHi-C P-P network", xlab = "Abundance", ylab = "ChAs")
  text((ab_PCHiC_PP), jitter(unlist(ass_PCHiC_PP), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
  abline(h = 0)
  legend("bottomright", legend = as.vector(unique(nodecats[, 2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[, 2]))], cex = 0.7)
  dev.off()

  pdf(paste(cell_type, "PCHiC Promoter-OEs assortativity vs abundance.pdf"))
  labs <- c(1:ncol(agchic))
  labs2 <- rep("", ncol(agchic))
  sel <- which(ass_PCHiC_POE > 0.04 | ass_PCHiC_POE < 0 | ab_PCHiC_POE > 0.03)
  labs2[sel] <- namesvec[sel]
  plot(ab_PCHiC_POE, unlist(ass_PCHiC_POE), pch = 20, col = cols4plot, main = "PCHi-C P-OE network", xlab = "Abundance", ylab = "ChAs")
  text((ab_PCHiC_POE), jitter(unlist(ass_PCHiC_POE), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
  abline(h = 0)
  legend("bottomright", legend = as.vector(unique(nodecats[, 2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[, 2]))], cex = 0.7)
  dev.off()


  pdf(paste(cell_type, "PCHiC Promoter-Promoter vs Promoter-OEs assortativity.pdf"))
  labs <- c(1:ncol(agchic))
  labs2 <- rep("", ncol(agchic))
  sel <- which(ass_PCHiC_POE > 0.04 | ass_PCHiC_POE < 0 | ass_PCHiC_PP > 0.03)
  labs2[sel] <- namesvec[sel]
  plot(unlist(ass_PCHiC_PP), unlist(ass_PCHiC_POE), pch = 20, col = cols4plot, main = "PCHi-C P-P vs PO ChAs", xlab = "P-P", ylab = "P-OE")
  text(unlist(ass_PCHiC_PP), jitter(unlist(ass_PCHiC_POE), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
  abline(a = 0, b = 1, h = 0, v = 0)
  legend("topleft", legend = as.vector(unique(nodecats[, 2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[, 2]))], cex = 0.7)
  dev.off()
}



##
setwd('/home/ngoc/Desktop/Ngoc_CRCT/2nd week/data_selected/csvfile')

temp_file<-list.files('/home/ngoc/Desktop/Ngoc_CRCT/2nd week/data_selected/csvfile')
temp_path<-paste('/home/ngoc/Desktop/Ngoc_CRCT/2nd week/data_selected/csvfile',temp_file,sep = '/')

  
  
for (temp in temp_path){
  assign(x=paste(str_match(temp, "RNA_seq_(\\w+).csv$")[2],'file',sep = '_'),value=read_csv(temp),envir = .GlobalEnv)
  }

for (cell_type in cell_types){
  big_function_RNA_seq_ass(cell_type,cell_type)
}





