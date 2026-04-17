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

# Helpers for transcriptomics metadata preparation, filtering, matching, and DESeq2 reruns.

prepare_tcga_transcriptomics_metadata <- function(expr_data, tcga_maf, clinical_path, smoking_file, nonsmoking_file) {
    # Build counts, gene metadata, and sample metadata for the TCGA-LUAD transcriptomics analysis.
    counts <- assay(expr_data)
    gene_data <- rowData(expr_data)
    metadata <- as.data.frame(colData(expr_data))

    if (!"barcode" %in% colnames(metadata)) {
        metadata$barcode <- colnames(counts)
    }
    if (!"patient" %in% colnames(metadata)) {
        metadata$patient <- substr(metadata$barcode, 1, 12)
    }
    rownames(metadata) <- colnames(counts)

    tcga_maf <- as.data.table(tcga_maf)
    tcga_maf_filtered <- tcga_maf[!Variant_Classification %in% c("Silent", "IGR", "Intron")]
    tcga_gps <- produce_genes_per_sample(tcga_maf_filtered, c("EGFR", "KRAS", "KEAP1"))
    tcga_gpp <- tcga_gps[, .(
        EGFR = ifelse(any(EGFR == "Mut"), "Mut", "WT"),
        KRAS = ifelse(any(KRAS == "Mut"), "Mut", "WT"),
        KEAP1 = ifelse(any(KEAP1 == "Mut"), "Mut", "WT")
    ), by = case_id]

    tcga_clinical <- fread(clinical_path)
    clinical_cols <- c(
        "case_id", "case_submitter_id", "gender", "race", "ethnicity",
        "ajcc_pathologic_stage", "tumor_stage", "ajcc_clinical_stage"
    )
    tcga_clinical <- unique(tcga_clinical[, ..clinical_cols])

    tcga_samples <- unique(tcga_clinical[, .(case_id, case_submitter_id)])
    tcga_samples <- tcga_gpp[tcga_samples, on = "case_id"]

    smoking_samples <- fread(smoking_file, header = FALSE)[[1]]
    nonsmoking_samples <- fread(nonsmoking_file, header = FALSE)[[1]]
    tcga_samples[case_id %in% smoking_samples, smoking_status := "smoking"]
    tcga_samples[case_id %in% nonsmoking_samples, smoking_status := "nonsmoking"]

    if (all(is.na(match(metadata$patient, tcga_samples$case_submitter_id)))) {
        stop("No matching case_submitter_id found for all patients in tcga_metadata")
    }

    sample_idx <- match(metadata$patient, tcga_samples$case_submitter_id)
    metadata[c("egfr_status", "kras_status", "keap1_status", "smoking_status")] <- as.data.frame(
        tcga_samples[sample_idx, .(EGFR, KRAS, KEAP1, smoking_status)]
    )

    clinical_idx <- match(metadata$patient, tcga_clinical$case_submitter_id)
    metadata$case_id <- tcga_clinical$case_id[clinical_idx]
    metadata$sex <- tolower(tcga_clinical$gender[clinical_idx])
    metadata$race <- tolower(tcga_clinical$race[clinical_idx])
    metadata$ethnicity <- tolower(tcga_clinical$ethnicity[clinical_idx])

    metadata$stage_raw <- tcga_clinical$ajcc_pathologic_stage[clinical_idx]
    missing_stage <- is.na(metadata$stage_raw) | metadata$stage_raw %in% c("", "'--", "not reported")
    metadata$stage_raw[missing_stage] <- tcga_clinical$tumor_stage[clinical_idx][missing_stage]
    missing_stage <- is.na(metadata$stage_raw) | metadata$stage_raw %in% c("", "'--", "not reported")
    metadata$stage_raw[missing_stage] <- tcga_clinical$ajcc_clinical_stage[clinical_idx][missing_stage]

    stage_lower <- tolower(metadata$stage_raw)
    metadata$stage_group <- fcase(
        is.na(stage_lower) | stage_lower %in% c("", "'--", "not reported"), "Unknown",
        grepl("stage iv|^iv", stage_lower), "IV",
        grepl("stage iii|^iii", stage_lower), "III",
        grepl("stage ii|^ii", stage_lower), "II",
        grepl("stage i|^i", stage_lower), "I",
        default = "Unknown"
    )

    metadata$sex[is.na(metadata$sex) | metadata$sex %in% c("", "'--", "not reported")] <- "unknown"
    metadata$ancestry_proxy <- fcase(
        !is.na(metadata$ethnicity) & metadata$ethnicity == "hispanic or latino", "hispanic_or_latino",
        !is.na(metadata$race) & metadata$race == "white", "non_hispanic_white",
        !is.na(metadata$race) & metadata$race == "black or african american", "non_hispanic_black",
        !is.na(metadata$race) & metadata$race == "asian", "non_hispanic_asian",
        default = "other_or_unknown"
    )

    metadata$egfr_status <- ifelse(metadata$sample_type == "Solid Tissue Normal", "WT", metadata$egfr_status)
    metadata$kras_status <- ifelse(metadata$sample_type == "Solid Tissue Normal", "WT", metadata$kras_status)
    metadata$keap1_status <- ifelse(metadata$sample_type == "Solid Tissue Normal", "WT", metadata$keap1_status)

    metadata$egfr_status <- factor(metadata$egfr_status, levels = c("WT", "Mut"))
    metadata$kras_status <- factor(metadata$kras_status, levels = c("WT", "Mut"))
    metadata$keap1_status <- factor(metadata$keap1_status, levels = c("WT", "Mut"))
    metadata$smoking_status <- factor(metadata$smoking_status, levels = c("nonsmoking", "smoking"))
    metadata$prior_treatment <- factor(metadata$prior_treatment, levels = c("No", "Yes"))
    metadata$sample_type <- factor(metadata$sample_type, levels = c("Solid Tissue Normal", "Primary Tumor"))
    metadata$sex <- factor(metadata$sex, levels = c("female", "male", "unknown"))
    metadata$stage_group <- factor(metadata$stage_group, levels = c("I", "II", "III", "IV", "Unknown"))
    metadata$ancestry_proxy <- factor(
        metadata$ancestry_proxy,
        levels = c("non_hispanic_white", "non_hispanic_black", "non_hispanic_asian", "hispanic_or_latino", "other_or_unknown")
    )

    list(counts = counts, gene_data = gene_data, metadata = metadata)
}

filter_transcriptomics_samples <- function(counts, gene_data, metadata, required_covariates,
                                           exclude_ffpe = TRUE, exclude_prior_treatment = TRUE,
                                           min_total_counts = 10) {
    # Apply the sample-level exclusions and low-count gene filter used by the transcriptomics analysis.
    keep <- rep(TRUE, nrow(metadata))
    filter_log <- data.table(step = "start", n_samples = nrow(metadata), n_genes = nrow(counts))

    if (length(required_covariates) > 0) {
        keep <- keep & stats::complete.cases(metadata[, required_covariates, drop = FALSE])
        filter_log <- rbind(filter_log, data.table(step = "required_covariates", n_samples = sum(keep), n_genes = nrow(counts)))
    }

    if (exclude_ffpe) {
        keep <- keep & !is.na(metadata$preservation_method) & metadata$preservation_method != "FFPE"
        filter_log <- rbind(filter_log, data.table(step = "exclude_ffpe", n_samples = sum(keep), n_genes = nrow(counts)))
    }

    if (exclude_prior_treatment) {
        keep <- keep & !is.na(metadata$prior_treatment) & metadata$prior_treatment == "No"
        filter_log <- rbind(filter_log, data.table(step = "exclude_prior_treatment", n_samples = sum(keep), n_genes = nrow(counts)))
    }

    metadata <- metadata[keep, , drop = FALSE]
    counts <- counts[, keep, drop = FALSE]

    gene_keep <- rowSums(counts) >= min_total_counts
    counts <- counts[gene_keep, , drop = FALSE]
    gene_data <- gene_data[gene_keep, ]
    filter_log <- rbind(filter_log, data.table(step = "gene_filter", n_samples = nrow(metadata), n_genes = nrow(counts)))

    subgroup_counts <- dcast(
        as.data.table(dplyr::count(as.data.frame(metadata), sample_type, smoking_status, egfr_status, keap1_status)),
        sample_type + smoking_status ~ egfr_status + keap1_status,
        value.var = "n",
        fill = 0
    )

    list(counts = counts, gene_data = gene_data, metadata = metadata, subgroup_counts = subgroup_counts, filter_log = filter_log)
}

match_smoking_groups <- function(metadata, focal_status_col = NULL, priority_covariates,
                                 summary_covariates = NULL, id_col = "barcode",
                                 smoking_col = "smoking_status",
                                 sample_type_col = "sample_type", seed = 1) {
    # Match smoking groups within sample type and an optional focal genotype stratum.
    # Inputs/outputs: metadata with one row per sample plus matching columns; returns selected IDs, matched metadata, balance summaries, and a stratum log.
    # Assumptions: nonsmoking samples define the reference cohort and priority_covariates are ordered from highest to lowest matching priority.
    stratum_cols <- c(sample_type_col, focal_status_col)
    if (is.null(summary_covariates)) {
        summary_covariates <- priority_covariates
    }
    match_cols <- unique(c(id_col, smoking_col, stratum_cols, priority_covariates, summary_covariates))
    metadata <- as.data.table(copy(as.data.frame(metadata)))
    metadata <- metadata[stats::complete.cases(metadata[, ..match_cols])]
    metadata <- metadata[do.call(order, c(lapply(c(stratum_cols, smoking_col, id_col), function(col) metadata[[col]]), list(na.last = TRUE)))]

    summarize_balance <- function(dt, covariates) {
        if (length(covariates) == 0) {
            return(data.table())
        }

        rbindlist(lapply(covariates, function(covariate) {
            if (is.numeric(dt[[covariate]])) {
                out <- dt[, .(n = .N, mean = mean(get(covariate)), sd = stats::sd(get(covariate))), by = smoking_col]
                setnames(out, smoking_col, "smoking_status")
                out[, `:=`(covariate = covariate, value = NA_character_, proportion = NA_real_)]
                return(out[, .(smoking_status, covariate, value, n, proportion, mean, sd)])
            }

            out <- dt[, .N, by = c(smoking_col, covariate)]
            setnames(out, c(smoking_col, covariate, "N"), c("smoking_status", "value", "n"))
            out[, proportion := n / sum(n), by = smoking_status]
            out[, `:=`(covariate = covariate, mean = NA_real_, sd = NA_real_)]
            out[, .(smoking_status, covariate, value, n, proportion, mean, sd)]
        }), fill = TRUE)
    }

    priority_covariates <- unique(priority_covariates)
    summary_covariates <- unique(c(stratum_cols, summary_covariates))
    balance_before <- summarize_balance(metadata, summary_covariates)

    relax_order <- rev(priority_covariates)
    priority_weights <- 2^((length(priority_covariates) - 1):0)
    ref_dt <- metadata[get(smoking_col) == "nonsmoking"]
    comp_dt <- metadata[get(smoking_col) == "smoking"]
    strata <- unique(ref_dt[, ..stratum_cols])

    selected_ref <- data.table()
    selected_comp <- data.table()
    matching_log <- data.table()

    set.seed(seed)

    for (i in seq_len(nrow(strata))) {
        stratum <- strata[i]
        ref_stratum <- ref_dt[stratum, on = c(sample_type_col, focal_status_col), nomatch = 0]
        comp_stratum <- comp_dt[stratum, on = c(sample_type_col, focal_status_col), nomatch = 0]
        if (nrow(comp_stratum) < nrow(ref_stratum)) {
            stop(sprintf(
                "Insufficient smoking samples for matching in stratum %s",
                paste(paste(stratum_cols, as.character(unlist(stratum)), sep = "="), collapse = ", ")
            ))
        }
        available <- copy(comp_stratum)
        matched_ref <- data.table()
        matched_comp <- data.table()

        for (j in seq_len(nrow(ref_stratum))) {
            if (nrow(available) == 0) {
                break
            }

            ref_row <- ref_stratum[j]
            candidate <- copy(available)

            for (relaxed_n in 0:length(relax_order)) {
                exact_covariates <- setdiff(priority_covariates, relax_order[seq_len(relaxed_n)])
                candidate <- copy(available)

                if (length(exact_covariates) > 0) {
                    for (covariate in exact_covariates) {
                        ref_value <- ref_row[[covariate]]
                        if (is.na(ref_value)) {
                            candidate <- candidate[is.na(get(covariate))]
                        } else {
                            candidate <- candidate[get(covariate) == ref_value]
                        }
                    }
                }

                if (nrow(candidate) > 0) {
                    break
                }
            }

            candidate[, match_score := 0L]
            for (covariate_i in seq_along(priority_covariates)) {
                covariate <- priority_covariates[covariate_i]
                weight <- priority_weights[covariate_i]
                ref_value <- ref_row[[covariate]]
                if (is.na(ref_value)) {
                    candidate[, match_score := match_score + weight * as.integer(is.na(get(covariate)))]
                } else {
                    candidate[, match_score := match_score + weight * as.integer(!is.na(get(covariate)) & get(covariate) == ref_value)]
                }
            }

            chosen <- candidate[order(-match_score, get(id_col))][1]
            matched_ref <- rbind(matched_ref, ref_row, fill = TRUE)
            matched_comp <- rbind(matched_comp, chosen[, !"match_score"], fill = TRUE)
            available <- available[get(id_col) != chosen[[id_col]]]
        }

        selected_ref <- rbind(selected_ref, matched_ref, fill = TRUE)
        selected_comp <- rbind(selected_comp, matched_comp, fill = TRUE)
        matching_log_row <- data.table(
            sample_type = stratum[[sample_type_col]],
            nonsmoking_n = nrow(ref_stratum),
            smoking_n = nrow(comp_stratum),
            matched_n = nrow(matched_ref)
        )
        if (!is.null(focal_status_col)) {
            matching_log_row[, focal_status := stratum[[focal_status_col]]]
        }
        matching_log <- rbind(matching_log, matching_log_row, fill = TRUE)
    }

    matched_metadata <- rbind(selected_ref, selected_comp, fill = TRUE)
    matched_metadata <- matched_metadata[do.call(order, c(lapply(c(stratum_cols, smoking_col, id_col), function(col) matched_metadata[[col]]), list(na.last = TRUE)))]
    balance_after <- summarize_balance(matched_metadata, summary_covariates)

    list(
        selected_ids = list(
            nonsmoking = selected_ref[[id_col]],
            smoking = selected_comp[[id_col]]
        ),
        matched_metadata = matched_metadata,
        balance_before = balance_before,
        balance_after = balance_after,
        matching_log = matching_log
    )
}

run_transcriptomics_deseq <- function(counts, gene_data, metadata, design_terms,
                                      primary_contrast = NULL, tumor_normal_contrast = NULL,
                                      pca_subset = NULL, pca_intgroup = NULL, shrink_type = "ashr") {
    # Build the DESeq2 object, VST object, and optionally the requested result tables for one transcriptomics design.
    metadata <- metadata[colnames(counts), , drop = FALSE]
    modeled_cols <- unique(unlist(strsplit(design_terms, ":", fixed = TRUE)))

    if (!all(modeled_cols %in% colnames(metadata))) {
        stop("Missing design columns in metadata")
    }

    for (modeled_col in modeled_cols) {
        if (length(unique(stats::na.omit(metadata[[modeled_col]]))) < 2) {
            stop(paste("Design column has fewer than two levels:", modeled_col))
        }
    }

    design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design_formula)

    if (is.null(pca_subset)) {
        pca_subset <- rep(TRUE, ncol(dds))
    }
    vsd <- vst(dds[, pca_subset], blind = TRUE)

    if (!is.null(primary_contrast) || !is.null(tumor_normal_contrast)) {
        dds <- DESeq(dds)
    }

    primary_results <- NULL
    primary_results_lfc <- NULL
    primary_table <- NULL
    tumor_normal_table <- NULL

    if (!is.null(primary_contrast)) {
        primary_results <- results(dds, contrast = primary_contrast)
        primary_results_lfc <- lfcShrink(dds, contrast = primary_contrast, res = primary_results, type = shrink_type)
        primary_results_lfc$stat <- primary_results[match(rownames(primary_results_lfc), rownames(primary_results)), "stat"]
        rownames(primary_results_lfc) <- gene_data[match(rownames(primary_results_lfc), rownames(gene_data)), "gene_name"]
        primary_table <- setDT(as.data.frame(primary_results_lfc), keep.rownames = TRUE)[]
        setnames(primary_table, "rn", "gene")
    }

    if (!is.null(tumor_normal_contrast)) {
        tumor_normal_results <- results(dds, contrast = tumor_normal_contrast)
        tumor_normal_results_lfc <- lfcShrink(dds, contrast = tumor_normal_contrast, res = tumor_normal_results, type = shrink_type)
        tumor_normal_results_lfc$stat <- tumor_normal_results[match(rownames(tumor_normal_results_lfc), rownames(tumor_normal_results)), "stat"]
        rownames(tumor_normal_results_lfc) <- gene_data[match(rownames(tumor_normal_results_lfc), rownames(gene_data)), "gene_name"]
        tumor_normal_table <- setDT(as.data.frame(tumor_normal_results_lfc), keep.rownames = TRUE)[]
        setnames(tumor_normal_table, "rn", "gene")
    }

    list(
        dds = dds,
        vsd = vsd,
        design_formula = design_formula,
        primary_results = primary_results,
        primary_results_lfc = primary_results_lfc,
        primary_table = primary_table,
        tumor_normal_table = tumor_normal_table
    )
}
