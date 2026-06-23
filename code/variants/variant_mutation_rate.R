library(data.table)
library(cancereffectsizeR)
library(Biostrings)
library(stringr)


# pick gene
# get rate
# collect all variants of interest for gene, reverse complementing as necessary
# find freq of that variant type (trinuc_mut) for the gene
# find #times that context is observed
# var_mut_rate = gene_rate * freq / #times-context-observed

variant_mutation_rate = function(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, ref_genome_version, variants=NULL, maf_df=NULL, samples=NULL){
  #' Purpose: assign mutation-rate mass to each observed SNV in one gene.
  #' Inputs:
  #'   gene                       HGNC symbol used for reporting/context.
  #'   gene_syn_mutation_rate     scalar synonymous mutation rate for the gene/cohort.
  #'   trinuc_proportion_matrix   named vector/table of cohort trinucleotide mutation proportions.
  #'   ref_genome_version         "hg19" or "hg38".
  #'   variants                   character vector of variant_ids formatted chr:pos_ref>alt.
  #'   maf_df, samples            optional source for variants when `variants` is NULL.
  #' Output: one row per unique variant_id with trinucleotide context and `mut_rate`.
  #' Assumptions:
  #'   The input variants define the complete variant universe for the requested
  #'   rate calculation. Subset-invariant rates therefore require callers to pass
  #'   the full default gene universe here, then subset by summing the returned rows.
  gene_x = gene

  if(ref_genome_version == "hg19"){
    reference_genome = BSgenome.Hsapiens.UCSC.hg19}
  else if(ref_genome_version == "hg38"){
    reference_genome = BSgenome.Hsapiens.UCSC.hg38}
  else{stop('Reference genome version must be one of `hg19` or `hg38`')}

  # print(gene_x)
  
  # certain variants have an NA variant_id, but none of these are in genes in our 1290 genes of interest.
  
  if(is.null(variants)){
    print('No variant list provided. Using only the variants observed within the gene from the MAF dataframe provided.')
    if(is.null(maf_df)){stop("MAF dataframe must be provided if variant list is not provided")}
    if(is.null(samples)){stop("List of samples must be provided if variant list is not provided.
                              If you want to include all samples, use cesa$samples$Unique_Patient_Identifier")}
    variants = unique(maf_df[Unique_Patient_Identifier %in% samples & top_gene == gene_x & 
                               variant_type == 'snv' & !is.na(variant_id), 
                             variant_id])
    if(length(variants) == 0){warning('There are no variants present in the dataset corresponding to this gene'); return(NULL)}
  } else if(length(variants) == 0){warning('The list of variants provided is empty.'); return(NULL)}
  
  unique_variants = unique(variants)
  
  tmp1 = str_split_fixed(unique_variants, pattern = ':', n=2)
  tmp2 = str_split_fixed(tmp1[,2], pattern= '_', n=2)
  variant_mut_rate_df = as.data.table(cbind(unique_variants, tmp1[,1],tmp2))
  colnames(variant_mut_rate_df) = c('variant_id','chr','pos','mut')
  variant_mut_rate_df[,pos := strtoi(pos)]
  
  #get trinuc context for variants
  
  #' 
  #' ASSUMPTION: all sequences taken from POSITIVE strand
  #' 
  variant_mut_rate_df[, trinuc := as.vector(as.character(getSeq(reference_genome, names = str_c('chr', chr), start = pos-1, width=3, strand = '+')))]
  
  #reverse complementing variants as necessary
  middle = str_c(substr(variant_mut_rate_df$trinuc,2,2),collapse="")
  revcomp_ind = c(str_locate_all(middle,'G|A')[[1]][,1])
  variant_mut_rate_df[revcomp_ind, ':=' (trinuc = as.character(reverseComplement(DNAStringSet(variant_mut_rate_df[revcomp_ind, trinuc]))), 
                              mut = paste0(as.character(reverseComplement(DNAStringSet(variant_mut_rate_df[revcomp_ind, substr(mut,1,1)]))),
                                           '>',
                                           as.character(reverseComplement(DNAStringSet(variant_mut_rate_df[revcomp_ind, substr(mut,3,3)])))))]
  
  #generating COSMIC type format for trinuc mutations
  variant_mut_rate_df[,trinuc_mut := paste0(substr(trinuc,1,1),'[', mut,']',substr(trinuc,3,3))]
  
  #gets trinuc composition for each variant
  variant_mut_rate_df[,trinuc_mut_prop := as.vector(trinuc_proportion_matrix[trinuc_mut])]
  
  variant_mut_rate_df[,gene_rate := gene_syn_mutation_rate]
  
  tri_nt_contexts = table(variant_mut_rate_df$trinuc)
  variant_mut_rate_df[, context_freq := as.numeric(tri_nt_contexts[trinuc])]
  
  #' #' This only works with the singular table returned by compute tri-nt contexts, 
  #' #' if we were to vectorize this, we would need to make it dataframe-workable
  #' tri_nt_contexts = compute_trinucleotide_contexts(gene_x)
  #' variant_mut_rate_df[,context_freq := tri_nt_contexts[trinuc]] # should use data.table notation
  
  #calculate mutation rate for variants
  variant_mut_rate_df[,mut_rate := gene_rate * trinuc_mut_prop / context_freq]
  
  return(variant_mut_rate_df)
}

gene_mutation_rate = function(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, ref_genome_version, variants=NULL, maf_df=NULL, samples=NULL){
  #' Purpose: sum variant-level mutation rates for one gene.
  #' Inputs/outputs/assumptions: same as variant_mutation_rate(); returns one scalar.
  variant_rates = variant_mutation_rate(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, ref_genome_version, variants, maf_df, samples)
  if (is.null(variant_rates)) {
    return(0)
  }
  return(sum(variant_rates$mut_rate))
  
  
  # variant_mut_rate_df = variant_mutation_rate(gene, gene_syn_mutation_rate, trinuc_proportion_matrix, variants, maf_df, samples)
  # if(!is.null(variant_mut_rate_df)){
  #   return(sum(variant_mut_rate_df$mut_rate))}
  # else{return(NA)}
}

compute_fixed_variant_mutation_rates = function(variants_per_gene,
                                                pan_data_gene_rates,
                                                smoking_gene_rates,
                                                nonsmoking_gene_rates,
                                                pan_trinuc,
                                                smoking_trinuc,
                                                nonsmoking_trinuc,
                                                ref_genome_version) {
  #' Purpose: compute one invariant mutation rate per observed variant_id.
  #' Inputs:
  #'   variants_per_gene       default, unfiltered universe with columns top_gene, variant_id.
  #'   *_gene_rates            per-cohort gene synonymous rates with columns gene, rate_grp_1.
  #'   *_trinuc                per-cohort trinucleotide mutation proportions.
  #'   ref_genome_version      "hg19" or "hg38".
  #' Output:
  #'   data.table keyed by top_gene and variant_id with pan/smoking/nonsmoking
  #'   variant rates. Dimensions are one row per unique (top_gene, variant_id).
  #' Assumptions:
  #'   Variant rates are defined relative to the full default variant universe for
  #'   each parent gene. Task-specific filters change only which fixed rows are
  #'   summed, not the rate assigned to any retained variant.

  required_cols = c("top_gene", "variant_id")
  missing_cols = setdiff(required_cols, colnames(variants_per_gene))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "variants_per_gene is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  variant_universe = unique(
    variants_per_gene[
      !is.na(top_gene) & !is.na(variant_id),
      .(top_gene, variant_id)
    ]
  )

  if (nrow(variant_universe) == 0) {
    return(data.table(
      top_gene = character(),
      variant_id = character(),
      chr = character(),
      pos = integer(),
      mut = character(),
      trinuc = character(),
      trinuc_mut = character(),
      context_freq = numeric(),
      pan_data_trinuc_mut_prop = numeric(),
      smoking_trinuc_mut_prop = numeric(),
      nonsmoking_trinuc_mut_prop = numeric(),
      pan_data = numeric(),
      smoking = numeric(),
      nonsmoking = numeric()
    ))
  }

  out_list = vector("list", length(unique(variant_universe$top_gene)))
  genes = sort(unique(variant_universe$top_gene))

  for (i in seq_along(genes)) {
    gene_x = genes[i]
    gene_variants = variant_universe[top_gene == gene_x, variant_id]

    pan_rate = pan_data_gene_rates[gene == gene_x, rate_grp_1]
    smoking_rate = smoking_gene_rates[gene == gene_x, rate_grp_1]
    nonsmoking_rate = nonsmoking_gene_rates[gene == gene_x, rate_grp_1]
    if (length(pan_rate) != 1 || length(smoking_rate) != 1 ||
        length(nonsmoking_rate) != 1) {
      stop(sprintf(
        "Expected exactly one synonymous mutation rate per cohort for %s; got pan=%d, smoking=%d, nonsmoking=%d",
        gene_x, length(pan_rate), length(smoking_rate), length(nonsmoking_rate)
      ))
    }

    pan_rates = variant_mutation_rate(
      gene_x, pan_rate, pan_trinuc, ref_genome_version, variants = gene_variants
    )[, .(
      top_gene = gene_x,
      variant_id,
      chr,
      pos,
      mut,
      trinuc,
      trinuc_mut,
      context_freq,
      pan_data_trinuc_mut_prop = trinuc_mut_prop,
      pan_data = mut_rate
    )]

    smoking_rates = variant_mutation_rate(
      gene_x, smoking_rate, smoking_trinuc, ref_genome_version, variants = gene_variants
    )[, .(
      variant_id,
      smoking_trinuc_mut_prop = trinuc_mut_prop,
      smoking = mut_rate
    )]

    nonsmoking_rates = variant_mutation_rate(
      gene_x, nonsmoking_rate, nonsmoking_trinuc, ref_genome_version, variants = gene_variants
    )[, .(
      variant_id,
      nonsmoking_trinuc_mut_prop = trinuc_mut_prop,
      nonsmoking = mut_rate
    )]

    gene_rates = merge(pan_rates, smoking_rates, by = "variant_id", all = FALSE)
    gene_rates = merge(gene_rates, nonsmoking_rates, by = "variant_id", all = FALSE)
    setcolorder(
      gene_rates,
      c(
        "top_gene",
        "variant_id",
        "chr",
        "pos",
        "mut",
        "trinuc",
        "trinuc_mut",
        "context_freq",
        "pan_data_trinuc_mut_prop",
        "smoking_trinuc_mut_prop",
        "nonsmoking_trinuc_mut_prop",
        "pan_data",
        "smoking",
        "nonsmoking"
      )
    )
    out_list[[i]] = gene_rates
  }

  variant_rate_table = rbindlist(out_list, use.names = TRUE)
  missing_rates = variant_rate_table[
    is.na(pan_data) | is.na(smoking) | is.na(nonsmoking)
  ]
  if (nrow(missing_rates) > 0) {
    shown = head(missing_rates[, .(top_gene, variant_id)], 20)
    stop(sprintf(
      "Fixed variant-rate table contains missing rates for %d variants. First affected rows: %s",
      nrow(missing_rates),
      paste(sprintf("(%s, %s)", shown$top_gene, shown$variant_id), collapse = "; ")
    ))
  }
  variant_rate_table
}

sum_fixed_variant_mutation_rates = function(variant_rate_table,
                                            variants_per_gene,
                                            genes) {
  #' Purpose: sum fixed per-variant rates over an analysis-specific variant set.
  #' Inputs:
  #'   variant_rate_table  output of compute_fixed_variant_mutation_rates().
  #'   variants_per_gene   analysis universe with columns top_gene, variant_id.
  #'   genes               ordered vector of genes to report.
  #' Output:
  #'   data.table with columns gene, pan_data, smoking, nonsmoking.
  #' Assumptions:
  #'   Each retained (top_gene, variant_id) must exist in the fixed rate table.
  #'   A gene with no retained variants receives rate 0 for all cohorts.

  required_cols = c("top_gene", "variant_id")
  missing_cols = setdiff(required_cols, colnames(variants_per_gene))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "variants_per_gene is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  rate_cols = c("top_gene", "variant_id", "pan_data", "smoking", "nonsmoking")
  missing_rate_cols = setdiff(rate_cols, colnames(variant_rate_table))
  if (length(missing_rate_cols) > 0) {
    stop(sprintf(
      "variant_rate_table is missing required columns: %s",
      paste(missing_rate_cols, collapse = ", ")
    ))
  }

  gene_order = data.table(gene = as.character(genes), gene_order = seq_along(genes))
  membership = unique(
    variants_per_gene[
      !is.na(top_gene) & !is.na(variant_id),
      .(top_gene, variant_id)
    ]
  )

  if (nrow(membership) > 0) {
    merged = merge(
      membership,
      variant_rate_table[, ..rate_cols],
      by = c("top_gene", "variant_id"),
      all.x = TRUE,
      sort = FALSE
    )
    missing_rates = merged[
      is.na(pan_data) | is.na(smoking) | is.na(nonsmoking)
    ]
    if (nrow(missing_rates) > 0) {
      shown = head(missing_rates[, .(top_gene, variant_id)], 20)
      stop(sprintf(
        "Analysis variant universe contains %d variants without fixed rates. First affected rows: %s",
        nrow(missing_rates),
        paste(sprintf("(%s, %s)", shown$top_gene, shown$variant_id), collapse = "; ")
      ))
    }

    summed = merged[, .(
      pan_data = sum(pan_data),
      smoking = sum(smoking),
      nonsmoking = sum(nonsmoking)
    ), by = .(gene = top_gene)]
  } else {
    summed = data.table(
      gene = character(),
      pan_data = numeric(),
      smoking = numeric(),
      nonsmoking = numeric()
    )
  }

  out = merge(gene_order, summed, by = "gene", all.x = TRUE, sort = FALSE)
  out[is.na(pan_data), pan_data := 0]
  out[is.na(smoking), smoking := 0]
  out[is.na(nonsmoking), nonsmoking := 0]
  setorder(out, gene_order)
  out[, gene_order := NULL]
  out
}
