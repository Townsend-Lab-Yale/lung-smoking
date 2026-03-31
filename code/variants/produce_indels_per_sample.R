library(data.table)

produce_indels_per_sample = function(cesa_maf,
                                     preloaded_maf,
                                     gene_list,
                                     samples_template,
                                     save_results = TRUE) {
  #' Build a sample-by-gene table for frameshift indels.
  #'
  #' Inputs
  #' ------
  #' cesa_maf:
  #'   cancereffectsizeR MAF-like table containing gene assignments and
  #'   variant identifiers.
  #' preloaded_maf:
  #'   preloaded MAF table containing `Variant_Classification`.
  #' gene_list:
  #'   vector of genes to retain as columns.
  #' samples_template:
  #'   `genes_per_sample`-style table used to enforce the same sample universe.
  #'
  #' Outputs
  #' -------
  #' A binary table with the same sample rows as `samples_template` and one
  #' column per gene. A gene is marked as 1 if the sample has a
  #' `Frame_Shift_Ins` or `Frame_Shift_Del` in that gene.
  #'
  #' Assumptions
  #' -----------
  #' - `variant_id` links `cesa_maf` and `preloaded_maf`.
  #' - Single-gene assignments are stored in `top_gene`, or can be recovered
  #'   from `genes` when there is only one mapped gene.

  cesa_maf = copy(cesa_maf)
  cesa_maf = cesa_maf[lengths(genes) == 1 & !is.na(genes) & is.na(top_gene),
                      top_gene := as.character(genes)]

  template_dt = as.data.table(copy(samples_template))
  gene_columns = setdiff(colnames(template_dt), c("Sample ID", "Source", "Panel"))

  frameshift_variants = preloaded_maf[
    is.na(problem)
    & Variant_Classification %in% c("Frame_Shift_Ins", "Frame_Shift_Del")
    & !is.na(variant_id),
    .(Unique_Patient_Identifier, variant_id)
  ]

  indel_maf = merge(
    cesa_maf,
    frameshift_variants,
    by = c("Unique_Patient_Identifier", "variant_id"),
    all = FALSE
  )

  indel_maf = indel_maf[
    top_gene %in% toupper(gene_list),
    .(`Sample ID` = Unique_Patient_Identifier, Gene = toupper(top_gene))
  ]
  indel_maf = unique(indel_maf)

  if (nrow(indel_maf) > 0) {
    indel_maf = merge(
      indel_maf,
      template_dt[, .(`Sample ID`, Source, Panel)],
      by = "Sample ID",
      all.x = TRUE,
      sort = FALSE
    )
    indel_maf = indel_maf[, .(`Sample ID`, Source, Panel, Gene)]

    missing_samples = fsetdiff(
      template_dt[, .(`Sample ID`, Source, Panel)],
      unique(indel_maf[, .(`Sample ID`, Source, Panel)])
    )
    subset_maf = rbind(indel_maf, missing_samples, fill = TRUE)
    subset_maf[is.na(subset_maf)] <- ""

    samples_indels = dcast(subset_maf,
                           `Sample ID` + Source + Panel ~ Gene,
                           fill = 0)
  } else {
    samples_indels = template_dt[, .(`Sample ID`, Source, Panel)]
  }

  if ("V1" %in% colnames(samples_indels)) {
    samples_indels$V1 = NULL
  }

  for (gene in gene_columns) {
    if (!(gene %in% colnames(samples_indels))) {
      samples_indels[, (gene) := 0]
    }
  }

  samples_indels = samples_indels[, c("Sample ID", "Source", "Panel", gene_columns), with = FALSE]

  tmp = samples_indels[, ..gene_columns]
  tmp[tmp != 0] <- 1
  samples_indels = cbind(samples_indels[, .(`Sample ID`, Source, Panel)], tmp)

  samples_indels = merge(
    template_dt[, .(`Sample ID`, Source, Panel)],
    samples_indels,
    by = c("Sample ID", "Source", "Panel"),
    all.x = TRUE,
    sort = FALSE
  )

  for (gene in gene_columns) {
    set(samples_indels,
        which(is.na(samples_indels[[gene]])),
        gene,
        0)
  }

  if (save_results) {
    fwrite(samples_indels, paste0(location_output, "indels_per_sample.txt"))
  }

  return(samples_indels)
}
