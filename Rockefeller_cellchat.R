# Pcyt REV interactome modeling
#try CellChat
library(CellChat)
library(patchwork)

{
  #Check differences between Pcyt and Ctr in cell chat
  Idents(Robj) <- Robj$Treatment
  
  #Pcyt Cellchat
  data.input <- GetAssayData(subset(Robj,  idents = 'Pcyt'), assay = "SCT", slot = "data") # normalized data matrix
  labels <- as.factor(subset(Robj,  idents = 'Pcyt')$CellType)
  meta <- data.frame(CellType = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat_ARD <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
  cellchat_ARD <- addMeta(cellchat_ARD, meta = meta, meta.name = "CellType")
  cellchat_ARD <- setIdent(cellchat_ARD, ident.use = "CellType") # set "labels" as default cell identity
  levels(cellchat_ARD@idents) # show factor levels of the cell labels
  
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  cellchat_ARD@DB <- CellChatDB
  
  cellchat_ARD <- subsetData(cellchat_ARD, features = NULL) # This step is necessary even if using the whole database
  cellchat_ARD <- identifyOverExpressedGenes(cellchat_ARD)
  cellchat_ARD <- identifyOverExpressedInteractions(cellchat_ARD)
  
  cellchat_ARD <- computeCommunProb(cellchat_ARD)
  cellchat_ARD <- filterCommunication(cellchat_ARD, min.cells = 5)
  
  cellchat_ARD <- computeCommunProbPathway(cellchat_ARD)
  cellchat_ARD <- aggregateNet(cellchat_ARD)
  groupSize <- as.numeric(table(cellchat_ARD@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat_ARD@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_ARD@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  mat <- cellchat_ARD@net$weight
  par(mfrow = c(3,5), xpd=TRUE)
  
  cellchat_ARD <- netAnalysis_computeCentrality(cellchat_ARD, slot.name = "netP") 
  netAnalysis_signalingRole_network(cellchat_ARD, width = 8, height = 2.5, font.size = 10)
  
  cellchat_ARD <- computeNetSimilarity(cellchat_ARD, type = "functional")
  #cellchat_ARD <- computeNetSimilarity(cellchat_ARD, type = "structural")
  #cellchat_ARD <- netEmbedding(cellchat_ARD, type = "functional")
  #cellchat_ARD <- netEmbedding(cellchat_ARD, type = "structural")
  
  #Ctr Cellchat
  data.input <- GetAssayData(subset(Robj,  idents = 'Ctr'), assay = "SCT", slot = "data") # normalized data matrix
  labels <- as.factor(subset(Robj,  idents = 'Ctr')$CellType)
  meta <- data.frame(CellType = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat_Control <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
  cellchat_Control <- addMeta(cellchat_Control, meta = meta, meta.name = "CellType")
  cellchat_Control <- setIdent(cellchat_Control, ident.use = "CellType") # set "labels" as default cell identity
  levels(cellchat_Control@idents) # show factor levels of the cell labels
  
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  cellchat_Control@DB <- CellChatDB
  
  cellchat_Control <- subsetData(cellchat_Control, features = NULL) # This step is necessary even if using the whole database
  cellchat_Control <- identifyOverExpressedGenes(cellchat_Control)
  cellchat_Control <- identifyOverExpressedInteractions(cellchat_Control)
  
  cellchat_Control <- computeCommunProb(cellchat_Control)
  cellchat_Control <- filterCommunication(cellchat_Control, min.cells = 5)
  
  cellchat_Control <- computeCommunProbPathway(cellchat_Control)
  cellchat_Control <- aggregateNet(cellchat_Control)
  groupSize <- as.numeric(table(cellchat_Control@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat_Control@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_Control@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  mat <- cellchat_Control@net$weight
  par(mfrow = c(3,5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  
  cellchat_Control <- netAnalysis_computeCentrality(cellchat_Control, slot.name = "netP") 
  #cellchat_Control <- computeNetSimilarity(cellchat_Control, type = "functional")
  #cellchat_Control <- netEmbedding(cellchat_Control, type = "functional")
  #cellchat_Control <- computeNetSimilarity(cellchat_Control, type = "structural")
  #cellchat_Control <- netEmbedding(cellchat_Control, type = "structural")
}
#merge objects and evaluate differences
#differential
#saveRDS(cellchat_ARD, file = "/Users/vyom/analysis/ARD_Cancer_analysis_output/Pcyt/ARD_ARD_CANCER_cellchat.rds")
#saveRDS(cellchat_Control, file = "/Users/vyom/analysis/ARD_Cancer_analysis_output/Ctr/Control_ARD_CANCER_cellchat.rds")

#cellchat_ARD <- readRDS(paste0("/Users/vyom/analysis/ARD_Cancer_analysis_output/Pcyt/ARD_ARD_CANCER_cellchat.rds"))
#cellchat_Control <- readRDS(paste0("/Users/vyom/analysis/ARD_Cancer_analysis_output/Ctr/Control_ARD_CANCER_cellchat.rds"))

object.list <- list(Ctr = cellchat_Control, Pcyt = cellchat_ARD)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave(file = paste0('Differential_interactions_barplot.pdf'), width=10, height=5, units="in")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, sources.use = c(13))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", sources.use = c(13))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T)
ggsave(file = paste0('Differential_interactions_network.pdf'), width=10, height=5, units="in")


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
ggsave(file = paste0('Differential_interactions_heatmap.pdf'), width=10, height=5, units="in")


#comparing Major soures and targets in 2d space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # Ctr the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
ggsave(file = paste0('Differential_interaction_strength.pdf'), width=10, height=5, units="in")

Robj$CellType
cellchat@idents
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Th1 CD4 T")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endo.")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cancer", signaling.label = c('MHC-I', 'MHC-II', 'FN1', 'SPP1', 'COLLAGEN', 'SEMA4', 'PECAM2', 'PECAM1', 'CDH', 'JAM', 'ADGRE', 'APP', 'MIF', 'THBS','LAMININ','IL2'))          
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 T", signaling.label = c('MHC-I', 'MHC-II'))
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cancer")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "DP T")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mono.")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8 T")
patchwork::wrap_plots(plots = list(gg1,gg2))
ggsave(file = paste0('ARD_CANCER_Signaling_changes.pdf'),plot = gg1, width=7.5, height=5.5, units="in")


Robj$CellType
#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")
c('#1b9e77' ,'#d95f02')
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) 
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave(file = paste0('Information_flow.pdf'), width=8, height=6, units="in")

gg1 <- rankNet(cellchat, mode = "comparison",pairLR = pairLR.use.up ,stacked = T, do.stat = TRUE) + scale_fill_manual(values = c('#d95f02', '#1b9e77'))
gg2 <- rankNet(cellchat, mode = "comparison",pairLR = pairLR.use.up , stacked = F, do.stat = TRUE) + scale_fill_manual(values = c('#d95f02', '#1b9e77'))
gg1 + gg2
ggsave(file = paste0('Information_flow.pdf'), width=8, height=6, units="in")

gg1 <- rankNet(cellchat,sources.use= 2, mode = "comparison", stacked = T, do.stat = TRUE) 
gg2 <- rankNet(cellchat, sources.use= 2, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave(file = paste0('Information_flow_Th1.pdf'), width=6, height=5, units="in")

gg1 <- rankNet(cellchat,targets.use= 6, mode = "comparison", stacked = T, do.stat = TRUE) + scale_fill_manual(values = c('#d95f02', '#1b9e77'))
gg2 <- rankNet(cellchat, targets.use= 6, mode = "comparison", stacked = F, do.stat = TRUE) + scale_fill_manual(values = c('#d95f02', '#1b9e77'))
gg1 + gg2
ggsave(file = paste0('Information_flow.pdf'), width=8, height=6, units="in")


library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 25)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 25)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
ggsave(file = paste0('Signaling_patterns.pdf'), width=8, height=6, units="in")


netVisual_bubble(cellchat, sources.use = 13,  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 2,  comparison = c(1, 2), angle.x = 45)

netVisual_bubble(cellchat, targets.use = 1,  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, targets.use = 8,  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, targets.use = 9,  comparison = c(1, 2), angle.x = 45)


gg1 <- netVisual_bubble(cellchat, sources.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Pcyt", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 1,  comparison = c(1, 2),pairLR.use	= as.data.frame(c('MHC-I','MHC-II', 'Prostaglandin')))
#> Comparing communications on a merged object
gg2
tuft_interactoins <- gg2$data
gg1 + gg2
ggsave(file = paste0('interaction_DE_tuft_ARD.pdf'), width=15, height=10, units="in")

#use DE to verify changes

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Pcyt"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0, thresh.fc = 0, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Pcyt
net.up <- subsetCommunication(cellchat, net = net, datasets = "Pcyt",ligand.logFC = 0.01, receptor.logFC = 0.01)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Ctr, i.e.,downregulated in Pcyt
net.down <- subsetCommunication(cellchat, net = net, datasets = "Ctr",ligand.logFC = -0.01, receptor.logFC = -.01)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))         
#> Comparing communications on a merged object
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[1]], sources.use = 1, slot.name = 'netP', net = net.up, lab.cex = .2, small.gap = 0, big.gap = 0,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[2]], sources.use = 1, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway

#Show chord diagram
pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # Ctr the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram
pathways.show <- c("MHC-I") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram
pathways.show <- c("MHC-II") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

#visualize up signaling
# Chord diagram
par( xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 1,  slot.name = 'netP', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))       
netVisual_chord_gene(object.list[[1]], sources.use = 4,  slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway


