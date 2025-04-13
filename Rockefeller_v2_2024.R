#rockefeller 2024 analysis

library(Rmagic)
library(plyr)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(scales)
library(metR)
library(ggpubr)
library(rstatix)
library(svglite)
library(viridis)
library(reshape)
library(harmony)
library(nichenetr)
library(RColorBrewer)
library(Libra)
library(Nebulosa)
library(SeuratDisk)
library(reticulate)

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())
options(future.globals.maxSize = 4000 * 1024^5000)
Sys.setenv('R_MAX_VSIZE'=32000000000000000)
options(Seurat.object.assay.version = "v4")


Robj <- ReadParseBio("~/data/gokhan")
metadata <- read.csv('~/data/gokhan/cell_metadata.csv')
Robj <- CreateSeuratObject(Robj, project = "Robj")



sample_mapping = data.frame(bc1 = paste0(c(rep("A",12),rep("B",12),rep("C",12),rep("D",12)),rep(1:12,4)), samp = rep(c("Ctr1","Ctr2","Ctr3","Ctr4","Pcyt1","Pcyt2","Pcyt3","Pcty4"),each=6))
Robj$orig.ident = sample_mapping$samp[match(metadata$bc1_well,sample_mapping$bc1)]

VlnPlot(Robj, group.by = 'orig.ident', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
Robj <- subset(Robj, subset = nCount_RNA > 500 & nCount_RNA < 50000 & nFeature_RNA > 100 & nFeature_RNA < 6500)


#Robj <- RNAransform(Robj, ncells = 5000, variable.features.n = 5000)
Robj@assays$RNA <- NULL
DefaultAssay(Robj) <- 'RNA'
Robj <- NormalizeData(Robj)
Robj <- FindVariableFeatures(Robj, selection.method = "vst", nfeatures = 200)
Robj <- ScaleData(Robj)
#use_python("/Users/vyom/miniconda3/bin/python")
#Robj <- magic(Robj)

VlnPlot(Robj, group.by = 'orig.ident', layer = 'data', features = c("Ptprc", 'H2-D1'),pt.size = 0 ,ncol = 2)

VlnPlot(Robj, group.by = 'orig.ident', layer = 'count', features = c("nFeature_RNA", 'nCount_RNA'),pt.size = 0 ,ncol = 2)

VlnPlot(Robj, group.by = 'orig.ident',layer = 'data', features = c("nFeature_RNA", "nCount_RNA"),pt.size = 0 ,ncol = 2)


Robj <- RunPCA(Robj)


DimPlot(Robj, reduction = "pca")
ElbowPlot(Robj, ndims = 50, reduction = "pca")
Robj <- Robj %>% 
  RunHarmony("orig.ident", assay.use="RNA", plot_convergence = TRUE)

Robj <- RunUMAP(Robj, dims = 1:50, reduction='harmony', n.neighbors = 5, n.epochs = 500)
Robj <- FindNeighbors(Robj, reduction='harmony', dims = 1:50, k.param = 30, compute.SNN = TRUE)
Robj <- FindClusters(Robj, resolution = .75)
DimPlot(Robj, reduction = "umap", label = TRUE)

DimPlot(Robj, reduction = "umap", label = TRUE, group.by = 'orig.ident')

DotPlot_Sig <- unique(c('Epcam','Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Jchain','Ighm','Ighg1','S100a8','S100a9','Clec4d','Cd14','Cd68','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Col3a1','Sparc'))
DotPlot_Sig <- unique(c('Epcam','Ptprc','Cd3g','Cd3e','Cd8a','Cd8b1','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Jchain','Ighm','Ighg1','S100a8','S100a9','Ly6a','Cd14','Fcgr4','Cd68','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','Xcr1','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Itgax','Itgam','Itgae','Irf4','Klf4','Irf7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc', 'Cd11b', 'Csf1r', 'Cxcl9','Cxcl10','Ccl5','Cxcl13','Il10','Ccl1','Ccl17','Ifng'))     
DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il5','Il13','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Jchain','Ighm','Ighg1','Itgam','Itgax','S100a8','S100a9','Clec4d','Cd14','Cd68','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Col3a1','Sparc'))
FeaturePlot(Robj, features = 'Cd4')

DotPlot(Robj, features = DotPlot_Sig, assay = 'RNA') + labs(y= "Cell Treatment", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, vjust= .01), axis.title.x = element_blank(), axis.title.y = element_blank())
DotPlot(Robj, features = DotPlot_Sig, assay = 'RNA') + labs(y= "Cell Treatment", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 5)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())


my_vector <- c("Th1 CD4 T", "Mac.", "CD8 T", "CD8 T", "CD8 T", "B", "Mac.", "Mac.", "DC", "Mac.", "Mono.", "NK", "DP T", "Neut.", "DC", "DC", "Cancer", "CD8 T", "Mono.", "Stromal", "Th17 CD4 T", "Treg", "Mac.", "DC", "CD8 T", "Cancer", "B", "Neut.", "Mac.", "Cancer", "Treg", "Cancer", "Neut.")


Robj[["CellType"]] <- Idents(Robj)
names(my_vector) <- levels(Robj)
Robj <- RenameIdents(Robj, my_vector)
Robj[["CellType"]] <- Idents(Robj)
DimPlot(Robj)

unique(Robj$orig.ident)
my_vector <- c("Ctr",  "Ctr",  "Ctr",  "Ctr",  "Pcyt", "Pcyt", "Pcyt", "Pcyt")
unique(my_vector)
Idents(Robj) <- Robj$orig.ident
Robj[["Treatment"]] <- Idents(Robj)
names(my_vector) <- levels(Robj)
Robj <- RenameIdents(Robj, my_vector)
Robj[["Treatment"]] <- Idents(Robj)
DimPlot(Robj)

my_levels <- c("CD8 T", "Th1 CD4 T", "Th17 CD4 T", "Treg", "DP T", "NK", "B", "Neut.", "Mac.", "Mono.", "DC", "Stromal", "Cancer")
Robj$CellType <- factor(Robj$CellType, levels = my_levels)

DimPlot(Robj, reduction = "umap", label = TRUE, group.by = 'CellType')

prop.table(x = table(Robj$CellType, Robj$Treatment), margin = 2)
Prop_table<- prop.table(x = table(Robj$CellType, Robj$Treatment), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = FALSE)
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1, levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = c('Ctr','Pcyt'))
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  scale_fill_manual(values = c('#7CA1CC' ,'#FF4902')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Treatment") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="stacked_bar_prop.pdf",  width=3.5, height=2)

my_levels <-  c("CD8 T", "Th1 CD4 T", "Th17 CD4 T", "Treg", "DP T", "NK", "B", "Neut.", "Mac.", "Mono.", "DC", "Stromal", "Cancer")
my_levels1 <- c('Ctr','Pcyt')

Prop_table<- prop.table(x = table(Robj$CellType, Robj$orig.ident), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = FALSE)
Prop_Table
Prop_Table$Sample <- gsub('[[:digit:]]+', '', Prop_Table$Var2)
Prop_Table$replicate <- extract_numeric(Prop_Table$Var2)
Prop_Table$Sample[53:104] <- 'Pcyt'

Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1, levels = my_levels)
Prop_Table1$Sample <- factor(Prop_Table1$Sample, levels = my_levels1)

plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Sample)) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width=0.9), width = 0.5) +
  #ggpubr::geom_pwc( aes(x = Var1, y = Freq, group = Sample), method = "t_test", hjust = .5,  tip.length = 0.002, vjust = .6, y.position = .5, label.size = 2.5, label = "p.signif", bracket.nudge.y = .01, step.increase = 0.05, p.adjust.method = "BH") +
  #stat_compare_means(label = "p.signif", method = "wilcox.test") +
  #geom_bar(position="dodge", stat="identity", na.rm = TRUE) +  
  scale_fill_manual(values = c('#7CA1CC' ,'#FF4902')) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Treatment") + 
  ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="PCYT_proportion_with_SE_only_NK_sig.pdf",  width=3.25, height=2.25)

celltypes <- unique(Robj$CellType)

Idents(Robj) <- Robj$CellType
for(i in celltypes){
  subset_cell <- subset(Robj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "Pcyt", ident.2 = "Ctr",  test.use = "MAST",recorrect_umi=TRUE, logfc.threshold = .1, min.pct = .1, assay = 'RNA')
  DE_subset<- DE_subset[!grepl("Rik", rownames(DE_subset)), ]
  DE_subset<- DE_subset[!grepl("mt-", rownames(DE_subset)), ]
  DE_subset<- DE_subset[!grepl("Rpl", rownames(DE_subset)), ]
  DE_subset<- DE_subset[!grepl("Rps", rownames(DE_subset)), ]
  DE_subset<- DE_subset[!grepl("Atp", rownames(DE_subset)), ]
  DE_subset<- DE_subset[!grepl("AY0", rownames(DE_subset)), ]
  DE_subset<- DE_subset[!grepl("Gm", rownames(DE_subset)), ]
  write.csv(DE_subset, paste0(i,'_scrna_DE.csv')) 
  DE_subset<- DE_subset[!grepl("Xist", rownames(DE_subset)), ]
  if (max(-log10(DE_subset$p_val_adj)) < 321){
    ylim_num <- max(-log10(DE_subset$p_val_adj))} else {
      ylim_num <- 320}
  EnhancedVolcano(DE_subset, lab = rownames(DE_subset), x = 'avg_log2FC', y = 'p_val_adj', title = 'VP16 vs WT',
                  pCutoff = .05, FCcutoff = 0.1, xlim = c(min(DE_subset$avg_log2FC)-.05,max(DE_subset$avg_log2FC)+.05),ylim = c(0,ylim_num+5),
                  subtitle = i, gridlines.minor = FALSE, gridlines.major = FALSE)
  ggsave(file = paste0(i,'_DE_peak_volcano_SCRNA.pdf'), width=6, height=6, units="in")
}

celltypes <- c("CD8 T", "Th1 CD4 T", "Th17 CD4 T", "Treg", "DP T", "NK", "B", "Neut.", "Mono.",'Mac.', "DC", "Stromal", "Cancer")

for (i in celltypes) {
  res_format <- read.csv(paste0('/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Atomic 2024-2025/ARD immunotherapy SCRNA/DE/',i, '_immunotherapy_ARD_scrna_DE.csv'))
  res_format <- as.data.frame(res_format[complete.cases(res_format),])
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
  enriched <- enrichr(unique(res_format[(res_format$p_val_adj < .05 & res_format$avg_log2FC < -0.1 ),]$X), dbs)
  Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
  Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
  #plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
  write.csv(Enriched_filter, paste0(i,'_enrichR_Down_gsea.csv'))
  Sys.sleep(3)
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "MSigDB_Hallmark_2020")
  enriched <- enrichr(unique(res_format[(res_format$p_val_adj < .05 & res_format$avg_log2FC > 0.1 ),]$X), dbs)
  Enriched1 <- rbind(enriched[[1]],enriched[[2]],enriched[[3]],enriched[[4]] )
  Enriched_filter <- Enriched1[Enriched1$Adjusted.P.value < .05,]
  #plotEnrich(Enriched_filter, showTerms = 100, numChar = 100, y = "Count", orderBy = "Adjusted.P.value")
  write.csv(Enriched_filter, paste0(i,'_enrichR_Up_gsea.csv'))
  Sys.sleep(3)
}

DefaultAssay(Robj) <- 'SCT'
use_python("/Users/vyom/miniconda3/bin/python")
All_Genes <- c(rownames(Robj@assays$RNA))
gene.list <- unique(c('Ccl5','Ccl5','Ccl4','Ifi208',	'Hspa1a',	'Hspa1b',	'Jun',	'Klf2',	'Gzmk',	'H2-D1',	'Ifit3',	'Ifi209',	'Btg2',	'Ifit1bl1',	'Bst2',	'Ifit1',	'Fos',	'Isg15',	'Cxcr3',	'Ccl9',	'Irf7',	'Apobec1',	'Slfn1',	'Rps27',	'Bcl2',	'Ms4a4c',	'Ier5',	'Ly6c2',	'Itgb7',	'Ifi206', 'Ifng','Ifna1','Ifnb1','Cxcl10','Cd86','Ccl12','Icam1','Pdcd1','Cd28','Cxcr3','Ifngr1', 'Cd8a','Cd8b1','B2m', 'Cxcl16', 'H2-D1', 'H2-K1', 'H2-Ab1', 'Cd74', 'Ccl25', 'Il15', 'Il18', 'Il33', 'Il22ra1', 'Il18bp', 'Tnfrsf14', 'Igtp', 'Ifit1bl1', 'Ifit1bl2', 'Gbp7', 'Irgm2', 'Cd74','H2-D1','H2-Ab1','H2-Aa','H2-Q1','Ccr9','Ciita','H2-Eb1','H2-T3', 'Gzma','Lag3','Jund', 'Foxp3', 'Il10', 'Ebi3', 'Ctla4','Pdcd1','Il2ra', 'Il2rb', 'Entpd1', 'Tgfb1', 'Cd8a','Ifng','Klrg1','Gzmf','Gzmd','Gzme','Gzmk','Serpinb1a','Fasl', 'Gzmb','Gzma', 'Ccr7','Tnf','Il17ra','Bcl2', 'Lag3', 'Havcr2', 'Pdcd1', 'Cxcr3','Irf4','Klrb1', 'Rorc','Ccr6','Il10', 'Il9', 'Prf1'))      

gene.list <- intersect(All_Genes, gene.list)
Robj <- Rmagic::magic(Robj)

gene.list <- unique(c('Ifi208',	'Hspa1a',	'Hspa1b',	'Jun',	'Klf2',	'Gzmk',	'H2-D1',	'Ifit3',	'Ifi209',	'Btg2',	'Ifit1bl1',	'Bst2',	'Ifit1',	'Fos',	'Isg15',	'Cxcr3',	'Ccl9',	'Irf7',	'Apobec1',	'Slfn1',	'Rps27',	'Bcl2',	'Ms4a4c',	'Ier5',	'Ly6c2',	'Itgb7',	'Ifi206', 'Ifng','Ifna1','Ifnb1','Cxcl10','Cd86','Ccl12','Icam1','Pdcd1','Cd28','Cxcr3','Ifngr1', 'Cd8a','Cd8b1','B2m', 'Cxcl16', 'H2-D1', 'H2-K1', 'H2-Ab1', 'Cd74', 'Ccl25', 'Il15', 'Il18', 'Il33', 'Il22ra1', 'Il18bp', 'Tnfrsf14', 'Igtp', 'Ifit1bl1', 'Ifit1bl2', 'Gbp7', 'Irgm2', 'Cd74','H2-D1','H2-Ab1','H2-Aa','H2-Q1','Ccr9','Ciita','H2-Eb1','H2-T3', 'Gzma','Lag3','Jund', 'Foxp3', 'Il10', 'Ebi3', 'Ctla4','Pdcd1','Il2ra', 'Il2rb', 'Entpd1', 'Tgfb1', 'Cd8a','Ifng','Klrg1','Gzmf','Gzmd','Gzme','Gzmk','Serpinb1a','Fasl', 'Gzmb','Gzma', 'Ccr7','Tnf','Il17ra','Bcl2', 'Lag3', 'Havcr2', 'Pdcd1', 'Cxcr3','Irf4','Klrb1', 'Rorc','Ccr6','Il10', 'Il9', 'Prf1', 'Ccl5','Ccl5','Ccl4','Ifi208',	'Hspa1a',	'Hspa1b',	'Jun',	'Klf2',	'Gzmk',	'H2-D1',	'Ifit3',	'Ifi209',	'Btg2',	'Ifit1bl1',	'Bst2',	'Ifit1',	'Fos',	'Isg15',	'Cxcr3',	'Ccl9',	'Irf7',	'Apobec1',	'Slfn1',	'Rps27',	'Bcl2',	'Ms4a4c',	'Ier5',	'Ly6c2',	'Itgb7',	'Ifi206', 'Ifng','Ifna1','Ifnb1','Cxcl10','Cd86','Ccl12','Icam1','Pdcd1','Cd28','Cxcr3','Ifngr1', 'Cd8a','Cd8b1','B2m', 'Cxcl16', 'H2-D1', 'H2-K1', 'H2-Ab1', 'Cd74', 'Ccl25', 'Il15', 'Il18', 'Il33', 'Il22ra1', 'Il18bp', 'Tnfrsf14', 'Igtp', 'Ifit1bl1', 'Ifit1bl2', 'Gbp7', 'Irgm2', 'Cd74','H2-D1','H2-Ab1','H2-Aa','H2-Q1','Ccr9','Ciita','H2-Eb1','H2-T3', 'Gzma','Lag3','Jund', 'Foxp3', 'Il10', 'Ebi3', 'Ctla4','Pdcd1','Il2ra', 'Il2rb', 'Entpd1', 'Tgfb1', 'Cd8a','Ifng','Klrg1','Gzmf','Gzmd','Gzme','Gzmk','Serpinb1a','Fasl', 'Gzmb','Gzma', 'Ccr7','Tnf','Il17ra','Bcl2', 'Lag3', 'Havcr2', 'Pdcd1', 'Cxcr3','Irf4','Klrb1', 'Rorc','Ccr6','Il10', 'Il9', 'Prf1'))      
i = 'Ctse'
for(i in gene.list) {
  DefaultAssay(Robj) <- 'RNA'
  selected_cells <- names(Robj$CellType)
  vln_data <- FetchData(Robj,
                        vars = c(i,"Treatment", "CellType"),
                        cells = selected_cells,
                        layer = "data")
  vln_data$Treatment
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(Robj, split.by = "Treatment",group.by = 'CellType', features = i, pt.size = 0, assay = "RNA",  cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = CellType, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('Split_vln_', i, '.pdf'), plot=All, width=2.5, height=3, units="in")
}

DimPlot(Robj, reduction = "umap", group.by= 'CellType', pt.size = .0001)
ggsave(file = paste0('scrna_UMAP.pdf'), width=4.5, height=4, units="in")


source("~/analysis/split_violin.R")
gene.list <-  c('Ciita',	'H2-Aa',	'H2-Ab1',	'H2-Ea',	'H2-Eb1',	'H2-DMb1',	'H2-DMa',	'Cd74',	'Tap1',	'Psmb8',	'Sectm1b',	'Vnn1',	'Btnl2',	'Lrrc19',	'Reg4')


gene.list <- c('Ciita', 'Cd74','H2-T3', 'H2-D1', 'H2-Q1', 'H2-Q10', 'H2-Q2', 'H2-Ab1','H2-Aa')
gene.list <- unique(c('Cd8a','Cd8b1','Ifng','Prf1', 'Bcl2', 'Gzmk','Fasl','Tnf','Il17ra','Irf4'))
gene.list <- unique(c('Cd8a','Ifng','Klrg1','Gzmb','Gzma','Gzmf','Gzmd','Gzme','Gzmk','Serpinb1a','Fasl', 'Gzmb','Gzma', 'Ccr7','Tnf','Il17ra','Bcl2', 'Lag3', 'Havcr2', 'Pdcd1', 'Cxcr3','Irf4','Klrb1', 'Rorc','Ccr6','Il10', 'Il9', 'Prf1'))
gene.list <- unique(c('Cxcr3','Ifng', 'Gzmk', 'Tnf','Irf7','Cd69', 'Itga4'))
gene.list <- c( 'Epcam','Itga4', 'Itga1','Itgad','Itgae', 'Itgal', 'Itgam', 'Itgax', 'Cdh1', 'Cd69')
gene.list <- unique(c('Cd8a','Cd8b1','Ifng','Prf1', 'Bcl2', 'Gzmk','Fasl','Tnf','Il17ra','Irf4'))

gene.list <- unique(c('Cd8a','Gzmb','Gzma','Ifng','Irf1', 'Irf7', 'Ccl4','Ccl5', 'Itgae'))

gene.list <- unique(c('Mki67','Ccnd3','Pcna','Il2','Runx3', 'Bcl2', 'Bcl6'))
gene.list <- unique(c('Tnfrsf1b','Tgfbr1','Ifngr1', 'Ccl5','Cxcl9','Bcl2', 'H2-D1','H2-K1','Ccr5'))
gene.list <- unique(c('Gzma', 'Gzmb', 'Ifng', 'Prf1', 'Klrc1', 'Neat1', 'Ccr5',  'Ccl7', 'Klf2', 'Itgal','Fos', 'Junb'))
gene.list <- unique(c('Gzmb', 'Ifng', 'Prf1', 'Klrc1', 'Neat1', 'Ccr5',  'Ccl7', 'Klf2', 'Itgal','Fos', 'Junb'))
gene.list <- c('Ifng','Irf1', 'Gzmb','Prf1', 'Cxcr3', 'Ccr5', 'Cd44', 'Cd69', 'Itgal', 'Junb', 'Hmgb2', 'Rora')

gene.list <- c('Ifng', 'Gzmb', "Tnf", 'Cxcr3', "Tbx21", "Stat1", "Stat4", "Ccr5",  "Cd40lg",'Cd44', "Fasl")

DefaultAssay(subset_CD8) <- "MAGIC_SCT"
selected_cells <- names(subset_CD8$CellType[subset_CD8$CellType %in% c("CD8 T")])
vln_data <- FetchData(subset_CD8,
                      vars = c(gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
#ggsave(file = 'Robj_suppression_TReg_vln.pdf', plot=All, width=3, height=2, units="in")
ggsave(file = 'CD8_T_cell_effector_function_collab.pdf', plot=All, width=3, height=2, units="in")

DefaultAssay(subset_CD4) <- "MAGIC_SCT"
selected_cells <- names(subset_CD4$CellType[subset_CD4$CellType %in% c("Th1 CD4 T")])
vln_data <- FetchData(subset_CD4,
                      vars = c(gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
#ggsave(file = 'Robj_suppression_TReg_vln.pdf', plot=All, width=3, height=2, units="in")
ggsave(file = 'CD4_T_cell_effector_function_collab.pdf', plot=All, width=3, height=2, units="in")


DefaultAssay(subset_cancer) <- "MAGIC_SCT"
selected_cells <- names(subset_cancer$CellType[subset_cancer$CellType %in% c("Cancer")])
vln_data <- FetchData(subset_cancer,
                      vars = c(gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
#ggsave(file = 'Robj_suppression_TReg_vln.pdf', plot=All, width=3, height=2, units="in")
ggsave(file = 'Cancer_cell_function_collab.pdf', plot=All, width=3, height=2, units="in")



#saveRDS(Robj, file = "/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/Rockefeller_v2_2025.rds")
#Robj <- readRDS('/Users/vyom/CSHL Dropbox Team Dropbox/Vyom Shah/Seurat_Objects/Rockefeller_v2_2025.rds', refhook = NULL)

#RERUN NORMALIZATION!!!

Robj@assays$RNA <- NULL
Robj@assays$integrated <- NULL
Robj@assays$SCT <- NULL

DefaultAssay(Robj) <- 'RNA'
Robj <- NormalizeData(Robj)
Robj <- FindVariableFeatures(Robj, selection.method = "vst", nfeatures = 2000)
Robj <- ScaleData(Robj)

DefaultAssay(Robj) <- 'SCT'
Idents(Robj) <- Robj$CellType
subset_cancer <- subset(Robj,  idents = 'Cancer')
use_python("/Users/vyom/miniconda3/bin/python")
subset_cancer <- magic(subset_cancer)

Idents(Robj) <- Robj$CellType
subset_CD8 <- subset(Robj,  idents = 'CD8 T')
use_python("/Users/vyom/miniconda3/bin/python")
subset_CD8 <- magic(subset_CD8)

subset_CD4 <- subset(Robj,  idents = 'Th1 CD4 T')
use_python("/Users/vyom/miniconda3/bin/python")
subset_CD4 <- magic(subset_CD4)

#instead of running it on the whole object subset into cd8 and th1 and just make the plots - you can troubleshoot the magic situation later
