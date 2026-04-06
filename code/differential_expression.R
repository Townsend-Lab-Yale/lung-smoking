# ============================================ #
# Setup
# ============================================ #

.libPaths(c("variants/.Rlibs", .libPaths()))
source("pkg_installation.R")
runtime_cran_packages <- unique(c(plotting_packages, "ashr", "ggsignif"))
runtime_bioc_packages <- unique(differential_expression_packages)
runtime_packages <- unique(c(runtime_cran_packages, runtime_bioc_packages))
install_missing_packages("variants", runtime_packages)

installed_runtime_packages <- rownames(installed.packages())
missing_runtime_packages <- setdiff(runtime_packages, installed_runtime_packages)
if (length(missing_runtime_packages) > 0) {
    if (is.null(getOption("repos")) || identical(getOption("repos")[["CRAN"]], "@CRAN@")) {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
    }

    cran_missing_packages <- intersect(missing_runtime_packages, runtime_cran_packages)
    if (length(cran_missing_packages) > 0) {
        install.packages(cran_missing_packages, lib = .libPaths()[1], dependencies = TRUE)
    }

    bioc_missing_packages <- intersect(setdiff(runtime_packages, rownames(installed.packages())), runtime_bioc_packages)
    if (length(bioc_missing_packages) > 0) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", lib = .libPaths()[1], dependencies = TRUE)
        }
        BiocManager::install(bioc_missing_packages, lib = .libPaths()[1], ask = FALSE, update = FALSE)
    }
}

library(data.table)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(fgsea)

source("de_functions.R")

location_data <- "../data"
location_gene_expression <- file.path(location_data, "gene_expression")
location_msigdb <- file.path(location_data, "msigdb_gene_sets")

dir.create(location_msigdb, recursive = TRUE, showWarnings = FALSE)
msigdb_release_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs"
msigdb_gmt_files <- c(
    "h.all.v2024.1.Hs.symbols.gmt",
    "c2.cp.v2024.1.Hs.symbols.gmt",
    "c5.go.v2024.1.Hs.symbols.gmt",
    "c6.all.v2024.1.Hs.symbols.gmt"
)
for (msigdb_gmt_file in msigdb_gmt_files) {
    local_msigdb_file <- file.path(location_msigdb, msigdb_gmt_file)
    if (!file.exists(local_msigdb_file)) {
        download.file(
            url = paste(msigdb_release_url, msigdb_gmt_file, sep = "/"),
            destfile = local_msigdb_file,
            mode = "wb"
        )
    }
}

# ============================================ #
# Data acquisition
# ============================================ #

transcriptome_cache_dir <- file.path("GDCdata", "TCGA-LUAD", "Transcriptome_Profiling", "Gene_Expression_Quantification")
mutation_cache_dir <- file.path("GDCdata", "TCGA-LUAD", "Simple_Nucleotide_Variation", "Masked_Somatic_Mutation")

query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)
transcriptome_cache_available <- dir.exists(transcriptome_cache_dir) &&
    length(list.files(transcriptome_cache_dir,
                      pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
                      recursive = TRUE,
                      full.names = TRUE)) > 0
if (!transcriptome_cache_available) {
    GDCdownload(query)
}
expr_data <- GDCprepare(query)

query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
)
mutation_cache_available <- dir.exists(mutation_cache_dir) &&
    length(list.files(mutation_cache_dir,
                      pattern = "\\.wxs\\.aliquot_ensemble_masked\\.maf\\.gz$",
                      recursive = TRUE,
                      full.names = TRUE)) > 0
if (!mutation_cache_available) {
    GDCdownload(query)
}
tcga_maf <- GDCprepare(query)

# ============================================ #
# Metadata and mutation/smoking annotation
# ============================================ #

# Extract counts and metadata
tcga_counts <- assay(expr_data)
tcga_metadata <- colData(expr_data)
tcga_gene_data <- rowData(expr_data)

# Obtain mutation status for EGFR, KRAS, and KEAP1 for each patient
tcga_maf <- as.data.table(tcga_maf)
tcga_maf_filtered <- tcga_maf[!Variant_Classification %in% c("Silent","IGR","Intron")]
tcga_gps <- produce_genes_per_sample(tcga_maf_filtered, c("EGFR","KRAS","KEAP1"))
any(!(tcga_gps[, mean(EGFR == "Mut"), by = case_id]$V1 %in% c(0, 1)))
any(!(tcga_gps[, mean(KRAS == "Mut"), by = case_id]$V1 %in% c(0, 1)))
any(!(tcga_gps[, mean(KEAP1 == "Mut"), by = case_id]$V1 %in% c(0, 1)))
tcga_gps[, .(percent_mutant = mean(KRAS == "Mut"), N = .N), by = case_id][!percent_mutant %in% c(0, 1)][order(percent_mutant)]
tcga_gps[, .(percent_mutant = mean(KEAP1 == "Mut"), N = .N), by = case_id][!percent_mutant %in% c(0, 1)][order(percent_mutant)]
tcga_gpp <- tcga_gps[, .(EGFR = ifelse(any(EGFR == "Mut"), "Mut", "WT"), 
                         KRAS = ifelse(any(KRAS == "Mut"), "Mut", "WT"),
                         KEAP1 = ifelse(any(KEAP1 == "Mut"), "Mut", "WT")), by = case_id]

tcga_clinical <- fread(file.path(location_data, "luad_tcga", "clinical.tsv")) # tcga_clinical contained the case_submitter_id needed to join with tcga_metadata
tcga_samples <- unique(tcga_clinical[, .(case_id, case_submitter_id)])
tcga_samples <- tcga_gpp[tcga_samples, on = "case_id"] # left join gpp onto tcga_samples

# Obtain smoking status for each patient
smoking_samples <- fread("../data/smoking_sample_ids.txt", header = F)[[1]]
nonsmoking_samples <- fread("../data/nonsmoking_sample_ids.txt", header = F)[[1]]
tcga_samples[case_id %in% smoking_samples, smoking_status := "smoking"]
tcga_samples[case_id %in% nonsmoking_samples, smoking_status := "nonsmoking"]

# Add mutation and smoking history information to metadata
if (all(is.na(match(tcga_metadata$patient, tcga_samples$case_submitter_id)))) {
    stop("No matching case_submitter_id found for all patients in tcga_metadata")
}
tcga_metadata[c("egfr_status", "kras_status", "keap1_status", "smoking_status")] <- tcga_samples[match(tcga_metadata$patient, tcga_samples$case_submitter_id), c("EGFR", "KRAS", "KEAP1", "smoking_status")]

# Set mutation status for normal samples to WT
tcga_metadata$egfr_status <- ifelse(tcga_metadata$sample_type == "Solid Tissue Normal", "WT", tcga_metadata$egfr_status)
tcga_metadata$kras_status <- ifelse(tcga_metadata$sample_type == "Solid Tissue Normal", "WT", tcga_metadata$kras_status)
tcga_metadata$keap1_status <- ifelse(tcga_metadata$sample_type == "Solid Tissue Normal", "WT", tcga_metadata$keap1_status)

# Factor all relevant metadata columns to set reference
tcga_metadata$egfr_status <- factor(tcga_metadata$egfr_status, levels = c("WT", "Mut"))
tcga_metadata$kras_status <- factor(tcga_metadata$kras_status, levels = c("WT", "Mut"))
tcga_metadata$keap1_status <- factor(tcga_metadata$keap1_status, levels = c("WT", "Mut"))
tcga_metadata$smoking_status <- factor(tcga_metadata$smoking_status, levels = c("nonsmoking", "smoking"))
tcga_metadata$prior_treatment <- factor(tcga_metadata$prior_treatment, levels = c("No", "Yes"))
tcga_metadata$sample_type <- factor(tcga_metadata$sample_type, levels = c("Solid Tissue Normal", "Primary Tumor"))

# ============================================ #
# Sample filtering
# ============================================ #

filtered_cases <- !is.na(tcga_metadata$smoking_status) & 
                    !is.na(tcga_metadata$egfr_status) & 
                    !is.na(tcga_metadata$kras_status) &
                    !is.na(tcga_metadata$keap1_status) &
                    !is.na(tcga_metadata$sample_type) &
                    tcga_metadata$preservation_method != "FFPE" &
                    !is.na(tcga_metadata$prior_treatment) &
                    tcga_metadata$prior_treatment == "No"

tcga_metadata_clean <- tcga_metadata[filtered_cases,]
tcga_counts_clean <- tcga_counts[,filtered_cases]

# Filter out genes with low counts across all samples
tcga_counts_clean <- tcga_counts_clean[rowSums(tcga_counts_clean) >= 10,]

dcast(as.data.table(dplyr::count(as.data.frame(tcga_metadata_clean), sample_type, smoking_status, egfr_status, keap1_status)), sample_type+smoking_status~egfr_status+keap1_status, value.var='n',fill=0)

# ============================================ #
# EGFR DESeq2 analysis
# ============================================ #

dds <- DESeqDataSetFromMatrix(countData = tcga_counts_clean, colData = tcga_metadata_clean, design = ~ sample_type + smoking_status + egfr_status + kras_status + smoking_status:egfr_status + smoking_status:kras_status)

options(repr.plot.width = 10, repr.plot.height = 10)
vsd <- vst(dds[,dds$sample_type=="Primary Tumor"], blind = TRUE)  
plotPCA(vsd,intgroup=c("smoking_status", "egfr_status")) + theme_bw()

dds <- DESeq(dds)

options(repr.plot.width = 8, repr.plot.height = 8)
plotDispEsts(dds)

# Results for EGFR-Smoking status interaction term
results <- results(dds, contrast = list("smoking_statussmoking.egfr_statusMut"))
results_lfc <- lfcShrink(dds, contrast = list("smoking_statussmoking.egfr_statusMut"), res = results, type = "ashr")
results_lfc$stat <- results[match(rownames(results_lfc), rownames(results)),"stat"]
rownames(results_lfc) <- tcga_gene_data[match(rownames(results_lfc), rownames(tcga_gene_data)),"gene_name"]
summary(results_lfc)

res = setDT(as.data.frame(results_lfc), keep.rownames = TRUE)[]
setnames(res,'rn','gene')

# Tumor vs Normal
results_tvn <- results(dds, contrast = list("sample_type_Primary.Tumor_vs_Solid.Tissue.Normal"))
results_tvn_lfc <- lfcShrink(dds, contrast = list("sample_type_Primary.Tumor_vs_Solid.Tissue.Normal"), res = results_tvn, type = "ashr")
results_tvn_lfc$stat <- results_tvn[match(rownames(results_tvn_lfc), rownames(results_tvn)),"stat"]
rownames(results_tvn_lfc) <- tcga_gene_data[match(rownames(results_tvn_lfc), rownames(tcga_gene_data)),"gene_name"]
# summary(results_tvn_lfc)

res_tvn = setDT(as.data.frame(results_tvn_lfc), keep.rownames = TRUE)[]
setnames(res_tvn,'rn','gene')

# Combining EGFR-Smoking & TvN results
alpha <- 0.05
tvn_tmp <- res_tvn[order(match(gene, res$gene))][,.(gene, log2FoldChange, padj)]
setnames(tvn_tmp, c("log2FoldChange", "padj"), c("log2FoldChange_tvn", "padj_tvn"))
res <- tvn_tmp[res,on="gene"]
res[,fc_match:=ifelse(padj_tvn<alpha,sign(log2FoldChange)==sign(log2FoldChange_tvn),NA)]

# ============================================ #
# EGFR GSEA
# ============================================ #

egfr_t_stats = results_lfc$stat
names(egfr_t_stats) = rownames(results_lfc)
egfr_t_stats = egfr_t_stats[order(egfr_t_stats,decreasing=TRUE)]
egfr_t_stats = egfr_t_stats[!is.na(names(egfr_t_stats))]
egfr_t_stats = egfr_t_stats[!duplicated(names(egfr_t_stats))]

pathways.hallmark = gmtPathways(file.path(location_msigdb, 'h.all.v2024.1.Hs.symbols.gmt'))
pathways.canonical = gmtPathways(file.path(location_msigdb, 'c2.cp.v2024.1.Hs.symbols.gmt'))
pathways.ontology = gmtPathways(file.path(location_msigdb, 'c5.go.v2024.1.Hs.symbols.gmt'))
pathways.oncogenic = gmtPathways(file.path(location_msigdb, 'c6.all.v2024.1.Hs.symbols.gmt'))

# Interesting genes will compile genes that drove signal of GSEA for pathways of interests (EMT, WNT, TGF-beta, ECM/Adhesion, PI3K-AKT)
interesting_genes <- data.table()

options(repr.plot.width=20, repr.plot.height = 6)
fgseaRes <- fgsea(pathways = pathways.hallmark, 
                  stats = egfr_t_stats,
                  minSize  = 10,
                  maxSize  = 500)
p1 <- ggplot(fgseaRes[padj<0.1], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
p2 <- plotEnrichment(pathways.hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, egfr_t_stats)
plot_grid(p1,p2,nrow=1, rel_widths = c(1, 0.5))

# Manually identifying gene sets of interest and extracting leading edge genes (genes that drive GSEA signal) for each of those gene sets.

leading_edges <- fgseaRes$leadingEdge
names(leading_edges) <- fgseaRes$pathway
gene_intersections <- table(unlist(leading_edges[c('HALLMARK_WNT_BETA_CATENIN_SIGNALING','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_MYC_TARGETS_V2','HALLMARK_HEDGEHOG_SIGNALING','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_APICAL_JUNCTION','HALLMARK_G2M_CHECKPOINT')]))
gene_intersections <- gene_intersections[gene_intersections > 1]
res[gene %in% names(gene_intersections) & padj<0.1][order(log2FoldChange)]

tmp <- leading_edges[c('HALLMARK_WNT_BETA_CATENIN_SIGNALING','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_MYC_TARGETS_V2','HALLMARK_HEDGEHOG_SIGNALING','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_APICAL_JUNCTION','HALLMARK_G2M_CHECKPOINT')]
tmp_df <- as.data.table(unlist(tmp),keep.rownames=T)
setnames(tmp_df,c("V1","V2"),c("pathway","gene"))
tmp_df[,pathway:=gsub('[0-9]+$','',pathway)]

interesting_genes <- rbind(interesting_genes,tmp_df)

fgseaRes <- fgsea(pathways = pathways.oncogenic, 
                  stats = egfr_t_stats,
                  minSize  = 10,
                  maxSize  = 500)

leading_edges <- fgseaRes$leadingEdge
names(leading_edges) <- fgseaRes$pathway
gene_intersections <- table(unlist(leading_edges[fgseaRes[padj<alpha,pathway]]))
gene_intersections <- gene_intersections[gene_intersections > 1]
res[gene %in% names(gene_intersections) & padj<alpha][order(log2FoldChange)][1:10]

pathways_of_interest <- c('CYCLIN_D1_UP.V1_UP','TGFB_UP.V1_DN','KRAS.LUNG.BREAST_UP.V1_DN','PDGF_ERK_DN.V1_DN','WNT_UP.V1_UP','PTEN_DN.V2_UP','AKT_UP.V1_DN')
tmp <- leading_edges[pathways_of_interest]
tmp_df <- as.data.table(unlist(tmp),keep.rownames=T)
setnames(tmp_df,c("V1","V2"),c("pathway","gene"))
tmp_df[,pathway:=gsub('[0-9]+$','',pathway)]

interesting_genes <- rbind(interesting_genes,tmp_df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fgseaRes <- fgsea(pathways = pathways.ontology, 
                  stats = egfr_t_stats,
                  minSize  = 10,
                  maxSize  = 500)

leading_edges <- fgseaRes$leadingEdge
names(leading_edges) <- fgseaRes$pathway
gene_intersections <- table(unlist(leading_edges[fgseaRes[padj<alpha,pathway]]))
gene_intersections <- gene_intersections[gene_intersections > 1]
res[gene %in% names(gene_intersections) & padj<alpha][order(log2FoldChange)][1:10]

pathways_of_interest <- fgseaRes[padj<alpha & stringr::str_detect(pathway,"(TGF)|(TRANSFORMING_GROWTH_FACTOR)|(WNT)|(CATENIN)|(EMT)|(MESENCHYMAL_TRANSITION)|(ADHESION)|(EXTRACELLULAR_MATRIX)|(PI3)|(EGFR)"),pathway]
tmp <- leading_edges[pathways_of_interest]
tmp_df <- as.data.table(unlist(tmp),keep.rownames=T)
setnames(tmp_df,c("V1","V2"),c("pathway","gene"))
tmp_df[,pathway:=gsub('[0-9]+$','',pathway)]

interesting_genes <- rbind(interesting_genes,tmp_df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fgseaRes <- fgsea(pathways = pathways.canonical, 
                  stats = egfr_t_stats,
                  minSize  = 10,
                  maxSize  = 500)

leading_edges <- fgseaRes$leadingEdge
names(leading_edges) <- fgseaRes$pathway
gene_intersections <- table(unlist(leading_edges[fgseaRes[padj<alpha,pathway]]))
gene_intersections <- gene_intersections[gene_intersections > 1]
res[gene %in% names(gene_intersections) & padj<alpha][order(log2FoldChange)][1:10]

pathways_of_interest <- fgseaRes[padj<alpha & stringr::str_detect(pathway,"(TGF)|(TRANSFORMING_GROWTH_FACTOR)|(WNT)|(CATENIN)|(HEDGEHOG)|(EMT)|(MESENCHYMAL_TRANSITION)|(ADHESION)|(ECM)|(COLLAGEN)|(FIBRES)|(EXTRACELLULAR_MATRIX)|(PI3)|(EGFR)"), pathway]
tmp <- leading_edges[pathways_of_interest]
tmp_df <- as.data.table(unlist(tmp),keep.rownames=T)
setnames(tmp_df,c("V1","V2"),c("pathway","gene"))
tmp_df[,pathway:=gsub('[0-9]+$','',pathway)]

interesting_genes <- rbind(interesting_genes,tmp_df)

interesting_genes <- unique(interesting_genes)
interesting_genes <- res[,.(gene,padj,log2FoldChange)][interesting_genes,on="gene"][order(gene)]

# ============================================ #
# EGFR postprocessing
# ============================================ #

interesting_genes_subset <- interesting_genes[padj<alpha*4]

tgf_beta <- interesting_genes_subset[stringr::str_detect(pathway,"(TGF)|(TRANSFORMING)"), .(gene=unique(gene),pathway="TGF_BETA")]
wnt <- interesting_genes_subset[stringr::str_detect(pathway,"(WNT)|(CATENIN)"), .(gene=unique(gene),pathway="WNT")]
hedgehog <- interesting_genes_subset[stringr::str_detect(pathway,"(HEDGEHOG)"), .(gene=unique(gene),pathway="HEDGEHOG")]
# wnt_and_hedgehog <- interesting_genes_subset[stringr::str_detect(pathway,"(WNT)|(CATENIN)|(HEDGEHOG)"), .(gene=unique(gene),pathway="WNT_AND_HH")]
emt <- interesting_genes_subset[stringr::str_detect(pathway,"(EMT)|(MESENCHYMAL_TRANSITION)"), .(gene=unique(gene),pathway="EMT")]
ecm <- interesting_genes_subset[stringr::str_detect(pathway,"(ADHESION)|(ECM)|(COLLAGEN)|(FIBRES)|(EXTRACELLULAR_MATRIX)"), .(gene=unique(gene),pathway="TM / ECM")]
pi3k <- interesting_genes_subset[stringr::str_detect(pathway,"(PI3)|(EGF)|(AKT)|(MEK)"), .(gene=unique(gene),pathway="PI3K")]

all <- rbind(tgf_beta, wnt, emt, ecm)

# Manually editing pathway assignments based on literature search
all <- rbind(all, data.table(
    gene = c("COL9A1","TM4SF4","S1PR3"),
    pathway = c("TM / ECM")
))
all <- rbind(all, data.table(
    gene = c("SMAD3","H19", "CAMK2B", "FZD10-AS1", "MYLK4", "PCDHA6", "HHIP", "PRSS2", 
             "IQGAP1", "SOX21", "S100A7", "TM4SF4", "TNFSF11", "FGF18", "S1PR3", "CD36"),
    pathway = c("EMT")
))
all <- rbind(all, data.table(
    gene = c("FZD10-AS1", "IQGAP1", "TRO", "LGALS4", "MUC6", "PCDHA6", "TENM3", "UNC5D",
             "TM4SF4", "PCDHGA3", "CD36"),
    pathway = c("WNT")
))
all <- rbind(all, data.table(
    gene = c("SMAD3", "HHIP", "PRG4", "FBLN5", "SOX21", "S100A7", "EMILIN3", "TNFSF11", "FGF18",
             "S1PR3"),
    pathway = c("TGF_BETA")
))
# all <- rbind(all, data.table(
#     gene = c("MYLK4", "F2", "F7", "CHAD"),
#     pathway = c("Other")
# ))

all[gene == "HEY2", pathway := c("TGF_BETA", "EMT")]
all[gene == "IL17RD", pathway := c("WNT", "EMT")]

all[,pathway:=factor(pathway, levels=c("TM / ECM","WNT","TGF_BETA","EMT","Other"))]

options(repr.plot.width = 4, repr.plot.height = 4)
pathway_mat <- dcast(all,gene~pathway,value.var="pathway",fun.aggregate=length)
# Transpose the data.table to cluster columns based on row values
pathway_mat_t <- t(as.matrix(pathway_mat[, -1, with = FALSE]))
colnames(pathway_mat_t) <- pathway_mat$gene

# Perform hierarchical clustering
dist_matrix <- dist(pathway_mat_t, method = "binary")
hclust_res <- hclust(dist_matrix, method = "ward.D2")

# Plot the dendrogram
plot(hclust_res, main = "Hierarchical Clustering of Columns", xlab = "", sub = "", cex = 0.9)

# Order columns based on clustering
ordered_columns <- rownames(pathway_mat_t)[hclust_res$order]
pathway_mat <- pathway_mat[, c("gene", ordered_columns), with = FALSE]

norm_counts_mat <- counts(dds, normalized=TRUE)
rownames(norm_counts_mat) <- tcga_gene_data[match(rownames(norm_counts_mat), rownames(tcga_gene_data)),"gene_name"]
norm_counts_mat <- as.data.table(norm_counts_mat, keep.rownames = TRUE)
setnames(norm_counts_mat,'rn','gene')
counts <- melt(norm_counts_mat, id.vars = "gene", variable.name = "barcode", value.name = "count")
counts <- as.data.table(colData(dds)[,c("barcode","smoking_status", "egfr_status", "kras_status")])[counts, on="barcode"]
counts$barcode <- NULL

# compute median counts for each group
median_counts = counts[,.(median_count=median(count)),by=.(smoking_status,egfr_status,gene)]
# reshape data
median_counts_x <- dcast(median_counts,gene+smoking_status~egfr_status,value.var="median_count")
# compute log2 fold change from EGFR WT to Mutant
median_counts_x[,log2fc:=log2(Mut)-log2(WT)]
# order genes by log2 fold change
median_counts_x$gene <- factor(median_counts_x$gene, levels=res[order(log2FoldChange,decreasing = T), gene])
# add info on change in expression from normal to tumor
median_counts_x <- res_tvn[,.(gene,tvn_change=ifelse(padj<alpha,ifelse(log2FoldChange>0,'+','-'),'none'))][median_counts_x,on="gene"]
# factor for ordering
median_counts_x$tvn_change <- factor(median_counts_x$tvn_change, levels=c('-','none','+')) 
# reshape to compare change in expression between smokers and nonsmokers
median_counts_x <- dcast(median_counts_x, gene~smoking_status, value.var="log2fc")[,.(gene,ratio_of_fc=2**(nonsmoking-smoking))][median_counts_x, on="gene"]
# add adjusted p-values for interaction term between smoking status and EGFR mutation status
median_counts_x <- res[,.(gene,padj)][median_counts_x,on="gene"]

# same with mean counts to check agreement (DESeq2 uses mean counts)
mean_counts = counts[,.(mean_count=mean(count)),by=.(smoking_status,egfr_status,gene)]
mean_counts_x <- dcast(mean_counts,gene+smoking_status~egfr_status,value.var="mean_count")
mean_counts_x[,log2fc:=log2(Mut)-log2(WT)]
mean_counts_x$gene <- factor(mean_counts_x$gene, levels=res[order(log2FoldChange,decreasing = T), gene])
mean_counts_x <- res[,.(gene,greater_in_ns=log2FoldChange<0,tvn_change=ifelse(padj_tvn<alpha,ifelse(log2FoldChange_tvn>0,'+','-'),'none'),padj)][mean_counts_x,on="gene"]
mean_counts_x$tvn_change <- factor(mean_counts_x$tvn_change, levels=c('-','none','+'))

smoking_status_abbrev <- c("nonsmoking" = "NS", "smoking" = "ES")

# median_counts_x$gene_name <- glue::glue("<i style='color:{ifelse(median_counts_x$padj<alpha,ifelse(median_counts_x$ratio_of_fc>1, 'red', 'blue'),'grey30')};'>{median_counts_filtered$gene}</i>")
median_counts_x$gene_name <- glue::glue("<i style='color:{ifelse(median_counts_x$ratio_of_fc>1, 'red', 'blue')};'>{median_counts_x$gene}</i>")
gene_name_order <- median_counts_x[smoking_status=="nonsmoking"][order(ratio_of_fc>1,tvn_change,log2fc),gene_name]
median_counts_x[,gene_name:=factor(gene_name, levels=gene_name_order)]
median_counts_x <- median_counts_x[order(gene_name)]
median_counts_x[,smoking_status:=factor(as.character(smoking_status_abbrev[smoking_status]),levels=c("ES","NS"))]

# mean_counts_x$gene_name <- glue::glue("<i style='color:{ifelse(mean_counts_x$padj<alpha,ifelse(mean_counts_x$greater_in_ns, 'red', 'blue'),'grey30')};'>{mean_counts_filtered$gene}</i>")
mean_counts_x$gene_name <- glue::glue("<i style='color:{ifelse(mean_counts_x$greater_in_ns, 'red', 'blue')};'>{mean_counts_x$gene}</i>")
gene_name_order <- mean_counts_x[smoking_status=="nonsmoking"][order(greater_in_ns,tvn_change,log2fc),gene_name]
mean_counts_x[,gene_name:=factor(gene_name, levels=gene_name_order)]
mean_counts_x <- mean_counts_x[order(gene_name)]
mean_counts_x[,smoking_status:=factor(as.character(smoking_status_abbrev[smoking_status]),levels=c("ES","NS"))]

tmp <- rbind(median_counts_x[,.(method="median",gene,greater_in_ns=ratio_of_fc>1,padj)], mean_counts_x[,.(method="mean",gene,greater_in_ns,padj)])
tmp <- unique(tmp)

genes_to_plot <- interesting_genes[padj<alpha*2 & pathway!="GOBP_LEUKOCYTE_CELL_CELL_ADHESION" & !(gene %in% c("MYPN","GRIN2B")),unique(gene)]
genes_to_plot <- c(genes_to_plot, "FZD10-AS1","SOX21")#, "PIK3C2G", "SMAD1")
# Only include genes for which the differential expression signal is consistent between median and mean counts
agreed_genes <- dcast(tmp, gene+padj~method,value.var="greater_in_ns")[mean==median & gene %in% genes_to_plot,gene]
agreed_genes <- agreed_genes[!agreed_genes %in% c("F2","F7","PROC","DEFB1","PENK","UNC5D","TLL2","HLX","DENND2A","BRMS1L")]

median_counts_filtered <- median_counts_x[gene %in% agreed_genes]
mean_counts_filtered <- mean_counts_x[gene %in% agreed_genes]

all_x <- rbind(all, data.table(gene=setdiff(agreed_genes,all$gene),pathway='Other'))

# ============================================ #
# EGFR figure assembly
# ============================================ #

options(repr.plot.width = 3.2, repr.plot.height = 10)
# genes_to_plot <- c(genes_to_plot, "SMAD3", "HHIP", "FZD10", "PRSS2","IQGAP1", "SOX21")
# agreed_genes <- c(agreed_genes, "SMAD3", "HHIP", "FZD10", "PRSS2","IQGAP1", "SOX21")

p1 = ggplot(median_counts_filtered, aes(x=smoking_status,y=gene_name)) + 
        geom_tile(fill=ifelse(median_counts_filtered$tvn_change=='+','red',ifelse(median_counts_filtered$tvn_change=='-','blue','gray')),alpha=0.1) +
        geom_point(aes(fill=log2fc, size=log2(Mut)),shape=21,color="black",alpha=1) + 
        geom_text(aes(label=fcase(
                smoking_status=="NS" & padj < alpha,   "*", 
                smoking_status=="NS" & padj < alpha*2, "."
                #smoking_status=="NS" & padj < alpha*3, "."
            ), vjust=ifelse(padj<alpha, 0.725, 0.1)), 
            hjust=0.5, size=4, nudge_x=-0.5) +
        scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0, breaks=c(-1.5,-0.5,0.5,1.3), labels=c(-1.5,-0.5,0.5,1.3)) + 
        scale_x_discrete(labels=c("NS"="NS-LUAD", "ES"="ES-LUAD")) +
        labs(x ="",
            y="", # "Gene",
            size = "log2(median expression) in<br>*EGFR* mutant samples",
            # fill = "Change in expression\nfrom Normal to Tumor",
            fill = "log2FC from<br>*EGFR* WT to mutant") +
        theme_minimal() + 
        theme(axis.text.y = ggtext::element_markdown(size=12, hjust=0.5),
                axis.text.x = element_text(size=12, angle=45, hjust=1),
                panel.grid.major = element_blank())
p2 = ggplot(all_x[gene %in% median_counts_filtered$gene][, gene := factor(gene, levels = unique(median_counts_filtered$gene))], aes(x=pathway, y=gene)) + 
        geom_tile(aes(fill=pathway),color="black",width=0.75,height=0.75) +
        scale_fill_brewer(palette="Pastel2") +
        scale_x_discrete(labels=c("TM / ECM"="TM / ECM","WNT"="  Wnt",
                                    "TGF_BETA"=expression(paste("  TGF-",beta)),
                                    "EMT"="  EMT")) + 
        labs(x="", y="") +
        theme_bw() + 
        theme(panel.grid = element_blank(), panel.border = element_blank(),
                plot.margin = margin(, , , 5, "mm"),
                axis.text.y = element_blank(), axis.ticks = element_blank(),
                axis.text.x = element_text(size=12, angle=45, hjust=1),
                legend.position = "none", legend.title =  ggtext::element_markdown(size=12, hjust=0.5))

plot_grid(p2, p1 + theme(legend.position="none"), nrow = 1, rel_widths = c(1, 1.35))

# ============================================ #
# KEAP1 GSVA analysis
# ============================================ #

library(GSVA)
library(GSEABase)

hallmark_sets <- getGmt(file.path(location_msigdb, 'h.all.v2024.1.Hs.symbols.gmt'),
                        geneIdType=SymbolIdentifier(),
                        collectionType = BroadCollection(category="h"))
canonical_sets <- getGmt(file.path(location_msigdb, 'c2.cp.v2024.1.Hs.symbols.gmt'),
                        geneIdType=SymbolIdentifier(),
                        collectionType = BroadCollection(category="c2"))
combined_sets <- GeneSetCollection(c(hallmark_sets, canonical_sets))

options(repr.plot.width = 10, repr.plot.height = 10)
dds <- DESeqDataSetFromMatrix(countData = tcga_counts_clean, colData = tcga_metadata_clean, design = ~ sample_type + smoking_status + keap1_status + smoking_status:keap1_status)

vsd <- vst(dds[,dds$sample_type=="Primary Tumor"], blind = TRUE)
plotPCA(vsd,intgroup=c("smoking_status", "keap1_status")) + theme_bw()

vsd_copy <- copy(vsd)
vsd_copy_gene_names <- tcga_gene_data[match(rownames(vsd_copy),rownames(tcga_gene_data)),"gene_name"]
vsd_copy <- vsd_copy[!is.na(vsd_copy_gene_names) & !duplicated(vsd_copy_gene_names), ]
rownames(vsd_copy) <- vsd_copy_gene_names[!is.na(vsd_copy_gene_names) & !duplicated(vsd_copy_gene_names)]
assayNames(vsd_copy) <- "counts"

keap1_par <- gsvaParam(exprData=vsd_copy[,vsd_copy$sample_type=="Primary Tumor"],geneSets=combined_sets,kcdf="Gaussian",assay='counts')
keap1_es <- gsva(keap1_par)

gsva_res <- as.data.table(t(assay(keap1_es)),keep.rownames = T)
setnames(gsva_res,'rn','barcode')
gsva_res <- as.data.table(tcga_metadata_clean)[,.(barcode,sample_type,smoking_status,keap1_status)][gsva_res,on="barcode"]

smoking_status_abbrev <- c("nonsmoking" = "NS-LUAD", "smoking" = "ES-LUAD")
gsva_res_to_plot <- gsva_res[sample_type=="Primary Tumor"]
gsva_res_to_plot[,smoking_status:=as.character(smoking_status_abbrev[smoking_status])]
gsva_res_to_plot[,group := paste0(smoking_status,'\n',keap1_status)]
gsva_res_to_plot[,group := factor(group, levels = c('NS-LUAD\nWT','NS-LUAD\nMut','ES-LUAD\nWT','ES-LUAD\nMut'))]

gsva_res_to_plot_2 <- copy(gsva_res_to_plot)
gsva_res_to_plot_2[,smoking_status:=fcase(
    smoking_status == "NS-LUAD", "NS-\nLUAD",
    smoking_status == "ES-LUAD", "ES-\nLUAD"
)]

options(repr.plot.width = 5, repr.plot.height = 6)
p = ggplot(gsva_res_to_plot_2,
        aes(x=group,
            y=KEGG_DNA_REPLICATION)) + 
    # geom_violin(scale="count",draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(aes(color=smoking_status),varwidth = T,outlier.shape=NA) +
    geom_jitter(aes(color=smoking_status),width=0.1) +
    geom_vline(xintercept=2.5,lty=2) +
    ggsignif::geom_signif(comparison = list(c('NS-LUAD\nWT','NS-LUAD\nMut'),
                                            c('ES-LUAD\nWT','ES-LUAD\nMut')),
                            map_signif_level = T,
                            textsize = 6) +
    scale_x_discrete(labels=rep(c("WT",'Mut'),2)) +
    labs(x="*KEAP1* status", y = "DNA Replication enrichment score") +
    theme_classic() +
    theme(axis.title.y = element_text(size=18), axis.text.y = element_text(size=12), axis.title.x = ggtext::element_markdown(size=18), axis.text.x = element_text(size=16),
            legend.position = c(0.35,0.75), legend.title=element_blank(), legend.text = element_text(size=14))
p
