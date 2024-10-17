#!/usr/bin/Rscript
#################
####Directory####
#################
DIR="/Users/aimedbio-a10174/Downloads/STREAM/NGS/"
Data="/Users/aimedbio-a10174/Downloads/STREAM/NGS/01_RNA/00_Data/"
Fig="/Users/aimedbio-a10174/Downloads/STREAM/NGS/01_RNA/01_Results/00_Figure/"
Table="/Users/aimedbio-a10174/Downloads/STREAM/NGS/01_RNA/01_Results/01_Table/"
###############
####library####
###############
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(edgeR)
library(DESeq2)
###################
#### Data load ####
###################
clinical.data <- read.csv(paste0(DIR, "Clinical.data.NGS.csv"), header=T)
clinical.data$Pathology.No2 <- gsub("SHH-activated", "SHH-activated/TP53-mutant", clinical.data$Pathology.No2)
clinical.data <- clinical.data[,-4]
colnames(clinical.data) <- c("ID", "Pathology", "Subpathology", "WHO", "EOR")
clinical.data$Gender <- c("M", "F", "F", "M", "M", "M", "F", "F", "M", "M", "M", "M", "F", "F", "F", "M", "F", "M", "M", "M", "F", "F")
gene.annotation.data <- read.table(paste0(DIR,"01_RNA/", "gencode.gene.annotation"), header=T)
rsem.data <- list()
for(i in 1:dim(clinical.data)[1]){
    tmp <- read.table(paste0(Data,clinical.data$ID[i],".genes.results"), sep="\t", header=T, check.names=F)
    rsem.data[[clinical.data$ID[i]]] <- tmp
}

# count data matrix
count.list <- list()
for(i in 1:dim(clinical.data)[1]){
    tmp <- read.table(paste0(Data,clinical.data$ID[i],".genes.results"), sep="\t", header=T, check.names=F)
    tmp.count <- tmp$expected_count
    count.list[[i]] <- tmp.count
}
tmp.count.mat <- do.call(cbind, count.list)
tmp.count.mat <- as.data.frame(tmp.count.mat)
count.mat <- as.data.frame(cbind(gene.annotation.data, tmp.count.mat))
colnames(count.mat) <- c("gene_id","genes",clinical.data$ID)
protein_coding.genes <- read.table(paste0(DIR, "01_RNA/protein.coding.genes"), sep="\t", header=T)

# length data matrix
length.list <- list()
for(i in 1:dim(clinical.data)[1]){
    tmp <- read.table(paste0(Data,clinical.data$ID[i],".genes.results"), sep="\t", header=T, check.names=F)
    tmp.length <- tmp$effective_length
    length.list[[i]] <- tmp.length
}
tmp.length.mat <- do.call(cbind, length.list)
tmp.length.mat <- as.data.frame(tmp.length.mat)
length.mat <- as.data.frame(cbind(gene.annotation.data, tmp.length.mat))
colnames(length.mat) <- c("gene_id","genes",clinical.data$ID)

groups <- factor(clinical.data$Pathology)
groups.other <- factor(gsub("Subependymal_giant_cell_astocytoma|Pediatric_type_diffuse_low_grade_glioma|Posterior_fossa_ependymoma","Other",groups))

medullo.groups <- factor(gsub("Supratentorial_ependymoma|Atypical_teratoid_rhabdoid_tumor|Subependymal_giant_cell_astocytoma|Pilocytic_astrocytoma|Pediatric_type_diffuse_low_grade_glioma|Posterior_fossa_ependymoma","Other",groups))
rhabdoid.groups <- factor(gsub("Supratentorial_ependymoma|Medulloblastoma|Subependymal_giant_cell_astocytoma|Pilocytic_astrocytoma|Pediatric_type_diffuse_low_grade_glioma|Posterior_fossa_ependymoma","Other",groups))
supraten.groups <- factor(gsub("Atypical_teratoid_rhabdoid_tumor|Medulloblastoma|Subependymal_giant_cell_astocytoma|Pilocytic_astrocytoma|Pediatric_type_diffuse_low_grade_glioma|Posterior_fossa_ependymoma","Other",groups))
#########################
#### Data processing ####
#########################
# DEGs.analysis <- function(count, group, gene_length){
#     res <- list()
#     # Create a DGEList object
#     only.counts <- as.matrix(count[,-(1:2)])
#     only.gene.length <- as.matrix(gene_length[,-(1:2)])
#     dge <- DGEList(counts = only.counts, groups = group, genes=count[,1:2])
#     dge$samples$group <- group
#     # Filter out lowly expressed genes
#     keep <- filterByExpr(dge)
#     dge <- dge[keep, , keep.lib.sizes = FALSE]
#     # Normalize the data
#     dge <- calcNormFactors(dge)
#     # Create a design matrix
#     design <- model.matrix(~0+group)
#     colnames(design) <- levels(group)
#     rownames(design) <- colnames(dge)
#     # Estimate common, trended, and tagwise dispersions
#     dge <- estimateDisp(dge, design)
#     # Fit the model and perform likelihood ratio tests
#     fit <- glmFit(dge, design)
#     lrt <- glmLRT(fit)
#     lrt$table$FDR <- p.adjust(lrt$table$PValue, 'fdr')
#     # Calculate library sizes in millions
#     lib_sizes <- dge$samples$lib.size / 1e6
#     # Calculate FPKM
#     #fpkm <- t(t(dge$counts) / lib_sizes) / (only.gene.length / 1000)  # Effective lengths should be in kilobases
#     res[["DEGs_analysis"]] <- lrt$table
#     #res[["FPKM"]] <- fpkm
#     return(res)
# }

# using DESeq2
DEGs.analysis <- function(count, group, gene_length, pc_genes){
    final.res <- list()
    # Create a object data
    only.counts <- as.matrix(count[,-(1:2)])
    only.gene.length <- as.matrix(gene_length[,-(1:2)])
    rsem_counts <- cbind(only.counts, only.gene.length)
    # Step 3: Prepare Sample Information
    sample_group <- c(as.character(group))
    col_data <- data.frame(
        row.names=colnames(only.counts),
        condition=sample_group
    )
    # Step 4: Create DESeq2 Dataset
    dds <- DESeqDataSetFromMatrix(
        countData = round(only.counts),
        colData = col_data,
        design = ~ condition
    )
    # Step 5: Run DESeq2 Analysis
    dds <- DESeq(dds)

    # Step 6: Extract Results
    res <- results(dds)
    res$gene_id <- count[, 1]
    res$genes <- count[,2]
    res <- res[order(res$padj), ]
    res <- res[res$baseMean!=0, ]
    res <- res[res$gene_id %in% protein_coding.genes$gene_id & res$genes %in% protein_coding.genes$Genes, ]
    res <- res[is.na(res$padj)!=TRUE,]
    
    # Step 7: Filter Significant Genes
    padj_threshold <- 0.05
    log2FC_threshold <- 1
    res_sig <- res[res$padj < padj_threshold & abs(res$log2FoldChange) > log2FC_threshold, ]

    # Step 8:  Calculate Nomalized expr values
    normalized_counts <- counts(dds, normalized=TRUE)
    normalized_counts <- data.frame(cbind(count[,1:2], normalized_counts))
    pc.normalized_counts <- normalized_counts[normalized_counts$gene_id %in% protein_coding.genes$gene_id, ]

    # Step 8: Calculate FPKM Values
    # Function to calculate FPKM using sample-specific effective lengths
    # calculate_fpkm <- function(gene.count, gene.length, pc.genes) {
    #     pc_counts <- gene.count[gene.count$gene_id %in% pc.genes$gene_id & gene.count$genes %in% pc.genes$Genes,]
    #     pc_lengths <- gene.length[gene.length$gene_id %in% pc.genes$gene_id & gene.length$genes %in% pc.genes$Genes,]
    #     only.pc.count <- pc_counts[,-(1:2)]
    #     only.pc.length <- pc_lengths[,-(1:2)]
    
    #     fpkm.list <- list()
    #     for(i in 1:ncol(only.pc.count)){
    #         rate <- only.pc.count[,i] / only.pc.length[,i]
    #         print(paste("SUM :", sum(rate, na.rm=TRUE)))
    #         final.rate <- rate / sum(rate, na.rm=TRUE) * 1e6
    #         fpkm.list[[i]] <- final.rate
    #     }
    #     fpkm <- do.call('cbind', fpkm.list)
    #     rownames(fpkm) <- pc_counts[,1]
    #     fpkm <- as.matrix(fpkm)
    #     return(fpkm)
    # }
    # # Calculate FPKM for each sample
    # # fpkm_values <- matrix(0, nrow = nrow(only.counts), ncol = ncol(only.counts))
    # # for (i in 1:ncol(only.counts)) {
    # #   fpkm_values[, i] <- calculate_fpkm(count, gene_length, protein_coding.genes, i)
    # # }
    # # colnames(fpkm_values) <- colnames(only.counts)
    # # rownames(fpkm_values) <- rownames(only.counts)
    # fpkm_values <- calculate_fpkm(count, gene_length, protein_coding.genes)


    # Combine FPKM values with results
    # res_fpkm <- cbind(as.data.frame(res), FPKM = fpkm_values)

    final.res[["DESeq2_Results"]] <- res
    final.res[["DESeq2_Significant_Results"]] <- res_sig
    final.res[["Normalized_Counts"]] <- pc.normalized_counts
    final.res[["DESeq2"]] <- dds

    # final.res[["FPKM"]] <- res_fpkm
    return(final.res)
}

medullo.DEGs <- DEGs.analysis(count.mat, medullo.groups, length.mat, protein_coding.genes)
rhabdoid.DEGs <- DEGs.analysis(count.mat, rhabdoid.groups, length.mat, protein_coding.genes)
supraten.DEGs <- DEGs.analysis(count.mat, supraten.groups, length.mat, protein_coding.genes)

medullo.DEGs.mat <- data.frame(medullo.DEGs$DESeq2_Results)
medullo.signif.DEGs.mat <- data.frame(medullo.DEGs$DESeq2_Significant_Results)
rhabdoid.DEGs.mat <- data.frame(rhabdoid.DEGs$DESeq2_Results)
rhabdoid.signif.DEGs.mat <- data.frame(rhabdoid.DEGs$DESeq2_Significant_Results)
supraten.DEGs.mat <- data.frame(supraten.DEGs$DESeq2_Results)
supraten.signif.DEGs.mat <- data.frame(supraten.DEGs$DESeq2_Significant_Results)

# DEGs condition = abs(log2FC) > 2 & Padj < 0.01
medullo.DEGs.mat$Category <- ifelse((abs(medullo.DEGs.mat$log2FoldChange)>1 & medullo.DEGs.mat$padj<0.05), "DEGs","NS")
rhabdoid.DEGs.mat$Category <- ifelse((abs(rhabdoid.DEGs.mat$log2FoldChange)>1 & rhabdoid.DEGs.mat$padj<0.05), "DEGs","NS")
supraten.DEGs.mat$Category <- ifelse((abs(supraten.DEGs.mat$log2FoldChange)>1 & supraten.DEGs.mat$padj<0.05), "DEGs","NS")

medullo.DEGs.mat$diffExpr <- ifelse((medullo.DEGs.mat$log2FoldChange>2 & medullo.DEGs.mat$padj<0.01), "Signif Up", ifelse((medullo.DEGs.mat$log2FoldChange > 1 & medullo.DEGs.mat$padj<0.05), "Up",
                                    ifelse((medullo.DEGs.mat$log2FoldChange < -2 & medullo.DEGs.mat$padj<0.01), "Signif Down", ifelse((medullo.DEGs.mat$log2FoldChange < -1 & medullo.DEGs.mat$padj<0.05), "Down", "NS"))))
rhabdoid.DEGs.mat$diffExpr <- ifelse((rhabdoid.DEGs.mat$log2FoldChange>2 & rhabdoid.DEGs.mat$padj<0.01), "Signif Up", ifelse((rhabdoid.DEGs.mat$log2FoldChange > 1 & rhabdoid.DEGs.mat$padj<0.05), "Up",
                                    ifelse((rhabdoid.DEGs.mat$log2FoldChange < -2 & rhabdoid.DEGs.mat$padj<0.01), "Signif Down", ifelse((rhabdoid.DEGs.mat$log2FoldChange < -1 & rhabdoid.DEGs.mat$padj<0.05), "Down", "NS"))))
supraten.DEGs.mat$diffExpr <- ifelse((supraten.DEGs.mat$log2FoldChange>2 & supraten.DEGs.mat$padj<0.01), "Signif Up", ifelse((supraten.DEGs.mat$log2FoldChange > 1 & supraten.DEGs.mat$padj<0.05), "Up",
                                    ifelse((supraten.DEGs.mat$log2FoldChange < -2 & supraten.DEGs.mat$padj<0.01), "Signif Down", ifelse((supraten.DEGs.mat$log2FoldChange < -1 & supraten.DEGs.mat$padj<0.05), "Down", "NS"))))


top.genes <- function(input){
    input <- input[order(input$log2FoldChange, decreasing=T), ]
    top <- c(head(input$genes, n=10), tail(input$genes, n=10))
    return(top)
}

medullo.DEGs.top10 <- top.genes(medullo.DEGs.mat)
rhabdoid.DEGs.top10 <- top.genes(rhabdoid.DEGs.mat)
supraten.DEGs.top10 <- top.genes(supraten.DEGs.mat)

medullo.DEGs.mat$genelabel <- ifelse(medullo.DEGs.mat$genes %in% medullo.DEGs.top10, medullo.DEGs.mat$genes, "")
rhabdoid.DEGs.mat$genelabel <- ifelse(rhabdoid.DEGs.mat$genes %in% rhabdoid.DEGs.top10, rhabdoid.DEGs.mat$genes, "")
supraten.DEGs.mat$genelabel <- ifelse(supraten.DEGs.mat$genes %in% supraten.DEGs.top10, supraten.DEGs.mat$genes, "")

library(ggrepel)

volcano.plot <- function(input, title){
    down <- "Down"
    sig.down <- "Signif Down"
    up <- "Up"
    sig.up <- "Signif Up"
    ns <- "NS"
    my_color <- setNames(c("#1d3557", "#a8dadc", "grey40", "#fed9b7", "#f07167"), c(sig.down, down, ns, up, sig.up))
    plot <- ggplot(data=input, aes(x=log2FoldChange, y=-log10(padj), col=diffExpr, label=genelabel)) + geom_point(size=0.8) +
            geom_text_repel(size=4) + scale_color_manual(values=my_color, name='diffExpr') +
            ggtitle(title) + theme(axis.text=element_text(size=12)) + theme_minimal()
    return(plot)
}

pdf(paste0(Fig,"Fig1.DEGs.volcano.plot.pdf"))
volcano.plot(medullo.DEGs.mat, "Medulloblastoma")
volcano.plot(rhabdoid.DEGs.mat, "rhabdoid tumor")
volcano.plot(supraten.DEGs.mat, "Supratentorial ependymoma")
dev.off()

write.table(medullo.DEGs$Normalized_Counts, paste0(Table,"Table1.DESeq2.Normalized.Counts.protein_coding.genes.txt"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(medullo.DEGs.mat[medullo.DEGs.mat$diffExpr %in% c("Signif Down", "Signif Up"), ]$genes, paste0(Table,"Table2.Medulloblastoma.Significant.DEGs.list.txt"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(rhabdoid.DEGs.mat[rhabdoid.DEGs.mat$diffExpr %in% c("Signif Down", "Signif Up"), ]$genes, paste0(Table,"Table2.Rhabdoid.Significant.DEGs.list.txt"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(supraten.DEGs.mat[supraten.DEGs.mat$diffExpr %in% c("Signif Down", "Signif Up"), ]$genes, paste0(Table,"Table2.Supratentolial.Significant.DEGs.list.txt"), sep="\t", quote=F, row.names=F, col.names=F)

medullo.DEGs.rank.mat <- medullo.DEGs.mat[, colnames(medullo.DEGs.mat) %in% c("log2FoldChange", "genes")][which(medullo.DEGs.mat$diffExpr %in% c("Signif Down", "Signif Up")), ]
medullo.DEGs.rank.mat <- cbind(medullo.DEGs.rank.mat, medullo.DEGs.rank.mat)
medullo.DEGs.rank.mat <- medullo.DEGs.rank.mat[,-1]
medullo.DEGs.rank.mat <- medullo.DEGs.rank.mat[,-3]
medullo.DEGs.high.rank.mat <- medullo.DEGs.rank.mat[abs(medullo.DEGs.rank.mat$log2FoldChange)>3, ]
write.table(medullo.DEGs.rank.mat, paste0(Table,"Table3.medullo.sifnificant.DEGs.fc.table.txt"), sep="\t", quote=F, row.names=F, col.names=F)
write.table(medullo.DEGs.high.rank.mat, paste0(Table,"Table3.medullo.sifnificant.DEGs.high.fc.table.txt"), sep="\t", quote=F, row.names=F, col.names=F)
write.csv(medullo.DEGs.high.rank.mat, paste0(Table,"Table3.medullo.sifnificant.DEGs.high.fc.table.csv"), row.names=F)
rhabdoid.DEGs.rank.mat <- rhabdoid.DEGs.mat[, colnames(rhabdoid.DEGs.mat) %in% c("log2FoldChange", "genes")][which(rhabdoid.DEGs.mat$diffExpr %in% c("Signif Down", "Signif Up")), ]
rhabdoid.DEGs.rank.mat <- cbind(rhabdoid.DEGs.rank.mat, rhabdoid.DEGs.rank.mat)
rhabdoid.DEGs.rank.mat <- rhabdoid.DEGs.rank.mat[,-1]
rhabdoid.DEGs.rank.mat <- rhabdoid.DEGs.rank.mat[,-3]
rhabdoid.DEGs.high.rank.mat <- rhabdoid.DEGs.rank.mat[abs(rhabdoid.DEGs.rank.mat$log2FoldChange)>3, ]
write.table(rhabdoid.DEGs.rank.mat, paste0(Table,"Table3.rhabdoid.sifnificant.DEGs.fc.table.txt"), sep="\t", quote=F, row.names=F, col.names=F)
supraten.DEGs.rank.mat <- supraten.DEGs.mat[, colnames(supraten.DEGs.mat) %in% c("log2FoldChange", "genes")][which(supraten.DEGs.mat$diffExpr %in% c("Signif Down", "Signif Up")), ]
supraten.DEGs.rank.mat <- cbind(supraten.DEGs.rank.mat, supraten.DEGs.rank.mat)
supraten.DEGs.rank.mat <- supraten.DEGs.rank.mat[,-1]
supraten.DEGs.rank.mat <- supraten.DEGs.rank.mat[,-3]
supraten.DEGs.high.rank.mat <- supraten.DEGs.rank.mat[abs(supraten.DEGs.rank.mat$log2FoldChange)>3, ]
write.table(supraten.DEGs.rank.mat, paste0(Table,"Table3.supraten.sifnificant.DEGs.fc.table.txt"), sep="\t", quote=F, row.names=F, col.names=F)



# PCA analysis
PCA <- function(input, total.group){
    PCA.res <- list()
    vsd <- vst(input, blind=FALSE)
    vsd$condition <- total.group
    pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
# PCA visualization
    p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size=3) + theme_classic() +
         xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed() + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=16))
    PCA.res[["PCA_vsd"]] <- vsd
    PCA.res[["PCA_Data"]] <- pcaData
    PCA.res[["Plot"]] <- p
    return(PCA.res)
}

medullo.DEGs.PCA <- PCA(medullo.DEGs$DESeq2, groups)
rhabdoid.DEGs.PCA <- PCA(rhabdoid.DEGs$DESeq2, groups)
supraten.DEGs.PCA <- PCA(supraten.DEGs$DESeq2, groups)

pdf(paste0(Fig, "Fig2.PCA.pdf"))
medullo.DEGs.PCA$Plot
rhabdoid.DEGs.PCA$Plot
supraten.DEGs.PCA$Plot
dev.off()


# prepare the meta-data (pathology information)
colData <- data.frame(row.names=colnames(count.mat)[3:24], pathology=groups)
colData.other <- data.frame(row.names=colnames(count.mat)[3:24], pathology=groups.other)
#DEGs.list <- unique(c(medullo.DEGs$DESeq2_Significant_Results$genes, rhabdoid.DEGs$DESeq2_Significant_Results$genes, supraten.DEGs$DESeq2_Significant_Results$genes))
DEGs.list <- unique(c(medullo.DEGs.top10, rhabdoid.DEGs.top10, supraten.DEGs.top10))
DEGs.normalized_counts <- medullo.DEGs$Normalized_Counts[medullo.DEGs$Normalized_Counts$genes %in% DEGs.list, ]
rownames(DEGs.normalized_counts) <- DEGs.normalized_counts$genes
DEGs.normalized_counts <- DEGs.normalized_counts[,-(1:2)]

# 병리학적 특성을 0과 1로 변환 (one-hot encoding)
pathology_matrix <- model.matrix(~ pathology - 1, data=colData.other)
colnames(pathology_matrix) <- levels(colData.other$pathology)

pathology_levels <- levels(colData.other$pathology)
DEGs.mean_expression <- sapply(pathology_levels, function(level) {
    rowMeans(DEGs.normalized_counts[, colData.other$pathology == level])
})

# NA 값을 제거하여 데이터 프레임으로 변환
DEGs.mean_expression_df <- as.data.frame(DEGs.mean_expression)
colnames(DEGs.mean_expression_df) <- pathology_levels

# NA 값 제거
DEGs.mean_expression_df <- na.omit(DEGs.mean_expression_df)

# 상관 행렬 계산
cor_matrix <- cor(DEGs.mean_expression_df)

# correlation matrix calculation
#cor_matrix <- cor(t(DEGs.normalized_counts), pathology_matrix)

# correlation matrix visualization
pdf(paste0(Fig, "Fig3.Correlation.DEGs.pdf"))
corrplot(cor_matrix, order="hclust", method="color", tl.col="black", tl.srt=90, is.corr=FALSE, 
         addCoef.col = "black", number.cex=1.5, col=colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F","#9E0142"))(20),
         title="Correlation between DEGs and Pathology")
dev.off()


# heatmap
major.significant.DEGs.list <- unique(c(medullo.DEGs.high.rank.mat$genes, rhabdoid.DEGs.high.rank.mat$gene, supraten.DEGs.high.rank.mat$genes))
filtered.DEGs.norm.count.mat <- medullo.DEGs$Normalized_Counts[medullo.DEGs$Normalized_Counts$genes %in% major.significant.DEGs.list, ]
rownames(filtered.DEGs.norm.count.mat) <- filtered.DEGs.norm.count.mat$genes
filtered.DEGs.norm.count.mat <- filtered.DEGs.norm.count.mat[,-(1:2)]
colnames(filtered.DEGs.norm.count.mat) <- clinical.data$ID


zscore <- function(input.mat){
    input.mat <- t(apply(input.mat, 1, scale))
    output <- na.omit(input.mat)
    return(output)
}

filtered.DEGs.norm.count.zscore.mat <- zscore(filtered.DEGs.norm.count.mat)

heatmap.plot <- function(input, clinical, title){
  clinical.info <- clinical
  pathology <- as.character(clinical.info$Pathology)
  pathology <- gsub("Pediatric_type_diffuse_low_grade_glioma|Subependymal_giant_cell_astocytoma|Posterior_fossa_ependymoma", "Others", pathology)
  who <- as.character(clinical.info$WHO)
  eor <- as.character(clinical.info$EOR)
  gender <- as.character(clinical.info$Gender)

  m1 <- as.matrix(input)

  hc <- dendsort(hclust(dist(m1), method="ward.D2"))
  cl.hc <- cutree(hc, k=4)
  split <- data.frame(cutree(hc, k=4))
  cluster <- as.character(cl.hc)

  ha <- HeatmapAnnotation(Pathology=pathology, WHO=who, EOR=eor, Gender=gender,
                          col=list(Pathology=c("Medulloblastoma"="#e76f51", "Pilocytic_astrocytoma"="#264653", "Atypical_teratoid_rhabdoid_tumor"="#2a9d8f", "Supratentorial_ependymoma"="#f4a261", "Pilocytic_astrocytoma"="#e9c46a", "Others"="#adb5bd"),
                                   WHO=c("1"="#ff99c8", "3"="#a9def9", "4"="#e4c1f9"),
                                   EOR=c("GTR"="#2ec4b6", "NTR"="#e71d36", "STR"="#ff9f1c"),
                                   Gender=c("M"="#b5e2fa", "F"="#eddea4")),
                          show_annotation_name=c(Pathology=FALSE, WHO=FALSE, EOR=FALSE, Gender=FALSE)
  )

  ha.row <- rowAnnotation(Cluster = paste0("C",cl.hc),
            #col = list(Cluster=c(C1="#f94144", C2="#f3722c", C3="#f8961e", C4="#f9844a", C5="#f9c74f", C6="#90be6d", C7="#43aa8b", C8="#4d908e", C9="#577590", C10="#277da1")),
            col = list(Cluster=c(C1="#f94144", C2="#f8961e", C3="#f9c74f", C4="#43aa8b")),
            show_annotation_name=c(Cluster=FALSE)
  )

  ht_list <- Heatmap(
                m1,
                name                   = "Scaled Expr",
                col                    = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
                bottom_annotation      = ha,
                left_annotation        = ha.row,
                use_raster             = TRUE,
                show_row_names         = FALSE,
                row_split              = split,
                row_title              = NULL,
                border                 = TRUE,
                show_column_names      = FALSE,
                clustering_method_rows = "ward.D2",
                cluster_rows           = TRUE,
                cluster_columns        = TRUE,
                column_title_gp        = gpar(fontsize = 8, fontface = "bold")
  )
  ht_list2 <- Heatmap(
                m1,
                name                   = "Scaled Expr",
                col                    = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
                bottom_annotation      = ha,
                left_annotation        = ha.row,
                use_raster             = TRUE,
                show_row_names         = FALSE,
                row_split              = split,
                row_title              = NULL,
                border                 = TRUE,
                show_column_names      = TRUE,
                column_names_gp        = gpar(fontzise = 6),
                clustering_method_rows = "ward.D2",
                cluster_rows           = TRUE,
                cluster_columns        = TRUE,
                column_title_gp        = gpar(fontsize = 8, fontface = "bold")
  )

  draw(ht_list, auto_adjust=FALSE,
       column_title=title,
       column_title_gp=gpar(fontsize=10, fontface="bold"),
       merge_legends=TRUE, heatmap_legend_side="right"
  )
  draw(ht_list2, auto_adjust=FALSE,
       column_title=title,
       column_title_gp=gpar(fontsize=10, fontface="bold"),
       merge_legends=TRUE, heatmap_legend_side="right"
  )
}
colnames(filtered.DEGs.norm.count.zscore.mat) <- clinical.data$ID
pdf(paste0(Fig,"Fig4.RNA.expr.landscape.heatmap.pdf"))
norm.count.heatmap <- heatmap.plot(filtered.DEGs.norm.count.zscore.mat, clinical.data, "STREAM Bulk RNA-seq")
dev.off()

medullo.norm.count <- filtered.DEGs.norm.count.mat[, colnames(filtered.DEGs.norm.count.mat) %in% clinical.data[clinical.data$Pathology=="Medulloblastoma",]$ID]
medullo.specific.DEGs.norm.count <- medullo.norm.count[rownames(medullo.norm.count) %in% medullo.DEGs.high.rank.mat$genes, ]
medullo.specific.DEGs.norm.count.zscore <- zscore(medullo.specific.DEGs.norm.count)
colnames(medullo.specific.DEGs.norm.count.zscore) <- colnames(medullo.specific.DEGs.norm.count)

specific.heatmap.plot <- function(input, clinical){
  res_list <- list()
  clinical.info <- clinical
  clinical.info <- clinical.info[clinical.info$ID %in% colnames(input), ]
  colnames(input) <- gsub("GINS-0029-", "", colnames(input))
  pathology <- as.character(clinical.info$Subpathology)
  who <- as.character(clinical.info$WHO)
  eor <- as.character(clinical.info$EOR)
  gender <- as.character(clinical.info$Gender)

  m1 <- as.matrix(input)

  hc <- dendsort(hclust(dist(m1), method="ward.D2"))
  cl.hc <- cutree(hc, k=4)
  split <- data.frame(cutree(hc, k=4))
  cluster <- as.character(cl.hc)

  ha <- HeatmapAnnotation(Pathology=pathology, WHO=who, EOR=eor, Gender=gender,
                          col=list(G.F=c("non-WNT/non-SHH"="#2a9d8f", "SHH-activated/TP53-mutant"="#e9c46a", "WNT-activated"="#e76f51"),
                                   WHO=c("1"="#ff99c8", "3"="#a9def9", "4"="#e4c1f9"),
                                   EOR=c("GTR"="#2ec4b6", "NTR"="#e71d36", "STR"="#ff9f1c"),
                                   Gender=c("M"="#b5e2fa", "F"="#eddea4")),
                         show_annotation_name=c(G.F=FALSE, WHO=FALSE, EOR=FALSE, Gender=FALSE)
  )

  ha.row <- rowAnnotation(Cluster = paste0("C",cl.hc),
            col = list(Cluster=c(C1="#f94144", C2="#f8961e", C3="#f9c74f", C4="#43aa8b")),
            show_annotation_name=c(Cluster=FALSE)
  )

  ht_list <- Heatmap(
                m1,
                name                   = "Scaled Expr",
                col                    = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
                bottom_annotation      = ha,
                left_annotation        = ha.row,
                use_raster             = TRUE,
                show_row_names         = FALSE,
                row_split              = split,
                row_title              = NULL,
                border                 = TRUE,
                show_column_names      = TRUE,
                column_names_gp        = gpar(fontzise = 6),
                clustering_method_rows = "ward.D2",
                cluster_rows           = TRUE,
                cluster_columns        = TRUE,
                column_title_gp        = gpar(fontsize = 6, fontface = "bold")
  )
  res_list[["cluster"]] <- split
  res_list[["plot"]] <- ht_list
  return(res_list)
}

medullo.heatmap.res <- specific.heatmap.plot(medullo.specific.DEGs.norm.count.zscore, clinical.data)

pdf(paste0(Fig,"Fig5.Medulloblastoma.specific.DEGs.heatmap.pdf"))
draw(medullo.heatmap.res$plot, auto_adjust=FALSE,
    column_title="Medulloblastoma RNA-seq heatmap",
    column_title_gp=gpar(fontsize=6, fontface="bold"),
    merge_legends=TRUE, heatmap_legend_side="right"
)
dev.off()
