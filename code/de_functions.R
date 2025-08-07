library(scales)

produce_genes_per_sample <- function(maf, gene_list) {
    maf <- maf[, .(case_id, Tumor_Sample_Barcode, Hugo_Symbol)]
    # some genes are mutated more than once but we only consider whether a gene is mutated or not
    maf <- unique(maf)
    # remove all genes that we do not want to consider (also removes samples whose mutations only occur in those genes)
    subset_maf <- maf[Hugo_Symbol %in% toupper(gene_list)]

    # adding the removed samples back
    subset_maf <- rbind(subset_maf,
        dplyr::setdiff(maf[, .(case_id, Tumor_Sample_Barcode)], subset_maf[, .(case_id, Tumor_Sample_Barcode)]),
        fill = T
    )

    #' Pivot to make a table where each row is a sample and the columns are genes
    #' with the value being a binary representation of whether the sample has the
    #' mutation or not
    samples_genes <- dcast(subset_maf, case_id + Tumor_Sample_Barcode ~ Hugo_Symbol,
        value.var = "Hugo_Symbol",
        fun.aggregate = function(x) ifelse(length(x) >= 1, "Mut", "WT")
    )
    samples_genes$`NA` <- NULL

    return(samples_genes)
}

get_smoker_nonsmoker_palette <- function() {
    c(
        "Ever-smoker" = hue_pal()(2)[1], "Smoker" = hue_pal()(2)[1], "smoking_plus" = hue_pal()(2)[1], "smoking" = hue_pal()(2)[1],
        "Never-smoker" = hue_pal()(2)[2], "nonsmoking_plus" = hue_pal()(2)[2], "nonsmoking" = hue_pal()(2)[2]
    )
}

get_counts <- function(gene_list, deseq_object, mapping_df = NULL, intgroup = NULL) {
    if (is.null(mapping_df)) {
        # All gene IDs are assumed to match those in the DESeq2 object
        filtered_gene_list <- gene_list[gene_list %in% rownames(deseq_object)]
        if (is.factor(filtered_gene_list)) filtered_gene_list <- filtered_gene_list[order(filtered_gene_list)]
        gene_ids <- filtered_gene_list
    } else {
        # Convert Hugo Symbol to respective gene ID using mapping dataframe
        filtered_gene_list <- gene_list[mapping_df[match(gene_list, mapping_df$gene_name), "gene_id"] %in% rownames(deseq_object)]
        if (is.factor(filtered_gene_list)) filtered_gene_list <- filtered_gene_list[order(filtered_gene_list)]
        gene_ids <- mapping_df[match(filtered_gene_list, mapping_df$gene_name), "gene_id"]
    }

    if (is.null(intgroup)) {
        # Extract all non-interaction-term variables from the design
        intgroup <- grep(":", strsplit(as.character(deseq_object@design)[2], " \\+ ")[[1]], invert = T, value = T)
        print(paste0("Setting intgroup to: ", paste0(intgroup, collapse = ", ")))
    }

    all_counts <- data.table()
    for (i in 1:length(filtered_gene_list)) {
        counts <- plotCounts(deseq_object,
            gene = gene_ids[i],
            intgroup = intgroup,
            returnData = TRUE
        )
        counts <- as.data.table(counts)[, gene := filtered_gene_list[i]]
        all_counts <- rbind(all_counts, counts)
    }
    all_counts$gene <- factor(all_counts$gene, levels = filtered_gene_list)
    return(all_counts)
}

get_counts_fast <- function(gene_list, deseq_object, mapping_df = NULL, intgroup = NULL, mc.cores = 8) {
    if (is.null(mapping_df)) {
        # All gene IDs are assumed to match those in the DESeq2 object
        filtered_gene_list <- gene_list[gene_list %in% rownames(deseq_object)]
        if (is.factor(filtered_gene_list)) filtered_gene_list <- filtered_gene_list[order(filtered_gene_list)]
        gene_ids <- filtered_gene_list
    } else {
        # Convert Hugo Symbol to respective gene ID using mapping dataframe
        filtered_gene_list <- gene_list[mapping_df[match(gene_list, mapping_df$gene_name), "gene_id"] %in% rownames(deseq_object)]
        if (is.factor(filtered_gene_list)) filtered_gene_list <- filtered_gene_list[order(filtered_gene_list)]
        gene_ids <- mapping_df[match(filtered_gene_list, mapping_df$gene_name), "gene_id"]
    }
    names(gene_ids) <- filtered_gene_list

    if (is.null(intgroup)) {
        # Extract all non-interaction-term variables from the design
        intgroup <- grep(":", strsplit(as.character(deseq_object@design)[2], " \\+ ")[[1]], invert = T, value = T)
        print(paste0("Setting intgroup to: ", paste0(intgroup, collapse = ", ")))
    }

    counts_list <- parallel::mclapply(gene_ids, plotCounts, dds=deseq_object, intgroup=intgroup, returnData=TRUE, mc.cores=mc.cores)
    all_counts <- rbindlist(counts_list,idcol=TRUE)
    setnames(all_counts, ".id", "gene")
    all_counts$gene <- factor(all_counts$gene, levels = filtered_gene_list)
    return(all_counts)
}

plot_counts <- function(counts, mean_or_median = "mean") {
    color_palette <- c(get_smoker_nonsmoker_palette(),
        "ES.+" = hue_pal()(2)[1], "ES.-" = hue_pal()(2)[1], "NS.+" = hue_pal()(2)[2], "NS.-" = hue_pal()(2)[2]
    )
    position_dodge <- 0.4

    count_stats <- counts[, .(median_count = median(count), mean_count = mean(count)),
        by = .(egfr_status, smoking_status, gene)
    ]

    p <- ggplot(counts, aes(x = egfr_status, y = count, color = smoking_status, group = smoking_status)) +
        geom_point(
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = position_dodge),
            alpha = 0.75
        ) +
        facet_grid(cols = vars(gene)) +
        scale_y_continuous(
            trans = "log2",
            breaks = scales::breaks_log(base = 2),
            labels = scales::label_log(base = 2)
        ) +
        labs(x = "EGFR status", y = "Normalized counts") +
        scale_color_manual(values = color_palette) +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(angle = 30, hjust = 1),
            axis.title.x = element_blank(),
            strip.text.x = element_text(size = 14),
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank()
        )

    if (mean_or_median == "mean") {
        p <- p +
            geom_line(
                data = count_stats,
                aes(x = egfr_status, y = mean_count, color = smoking_status, group = smoking_status),
                position = position_dodge(position_dodge)
            )
    } else if (mean_or_median == "median") {
        p <- p +
            geom_line(
                data = count_stats,
                aes(x = egfr_status, y = median_count, color = smoking_status, group = smoking_status),
                position = position_dodge(position_dodge)
            )
    } else {
        stop("mean_or_median must be either 'mean' or 'median'")
    }
    p
}