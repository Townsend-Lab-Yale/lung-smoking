library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)
library(dplyr)
library(tidyr)
library(glue)
library(repr)
library(cowplot)
library(stringr)

location_output = "../../output/"
location_cesR_output = glue("{location_output}/output_for_transfer/genes/cesR/")
location_variant_output = glue("{location_output}/output_for_transfer/genes/variant/")

# M = 1

load_M1_results = function(){
    # Average fluxes and selection intensities across genetic backgrounds (M=1 model)
    M1_gammas = fread("M1_gene_gammas.csv", header = T)
    M1_fluxes = fread("M1_gene_fluxes.csv", header = T)

    # mutation frequencies for each gene (note that it doesn't matter which mutation rate method is used)
    samples_per_combo_pd = fread(glue("{location_cesR_output}/M1/samples_per_combination_pan_data.csv"))
    samples_per_combo_pd$key = "pan_data"
    samples_per_combo_s = fread(glue("{location_cesR_output}/M1/samples_per_combination_smoking.csv"))
    samples_per_combo_s$key = "smoking"
    samples_per_combo_ns = fread(glue("{location_cesR_output}/M1/samples_per_combination_nonsmoking.csv"))
    samples_per_combo_ns$key = "nonsmoking"
    samples_per_combo_sp = fread(glue("{location_cesR_output}/M1/samples_per_combination_smoking_plus.csv"))
    samples_per_combo_sp$key = "smoking_plus"
    samples_per_combo_nsp = fread(glue("{location_cesR_output}/M1/samples_per_combination_nonsmoking_plus.csv"))
    samples_per_combo_nsp$key = "nonsmoking_plus"
    samples_per_combo = bind_rows(samples_per_combo_pd, samples_per_combo_s, samples_per_combo_ns, samples_per_combo_sp, samples_per_combo_nsp)
    colnames(samples_per_combo) = c("gene","0","1","key")
    samples_per_combo = 
        samples_per_combo %>%
        mutate(gene = gsub("[\\(\\)',]","",gene),
                freq = `1` / (`0` + `1`))

    # load mutation rates
    mus = fread("mutation_rates.csv", header = T)
    tmp = mus %>% filter(key %in% c("smoking","nonsmoking"))
    tmp = tmp %>% mutate(key = glue("{key}_plus")) 
    mus = mus %>% bind_rows(tmp)

    # combine all information (gammas, fluxes, mus, mutation frequencies)
    M1_results = 
        M1_gammas %>%
        left_join(M1_fluxes, by=c("method","key","gene")) %>%
        left_join(mus, by=c("method","key","gene")) %>%
        rename(mu = rate, mu_ci_low = rate_ci_low, mu_ci_high = rate_ci_high) %>%
        left_join(samples_per_combo %>% select(key, gene, freq), by=c("key","gene"))

    return(M1_results)
}


# Plot all information
plot_M1_results = function(df, dataset_key, mu_method, var_to_plot, show_freq_legend=TRUE, show_x_axis_title=TRUE){

    plotting_df = df %>% filter(key == dataset_key, method == mu_method)

    if(var_to_plot == "freq") {
        plot = plotting_df %>%
                ggplot(aes(x=reorder(gene, -freq), y=freq*100)) + 
                geom_point(aes(fill = log(gamma_mle)), pch = 21, color = "black", size = 4) + 
                scale_fill_viridis_c() +
                labs(y = "Frequency (%)", title = "Prevalence of mutations in lung adenocarcinoma", fill = "Log(Selection Intensity)") + 
                scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
                theme_minimal_grid() +
                
                theme(axis.title.x = element_blank(), 
                        axis.title.y = element_text(size = 18),
                        axis.text.x = element_text(size = 14),
                        plot.title = element_text(size = 20, hjust = 0.5),
                        panel.grid.major.x = element_line(size=0.2),
                        legend.position = c(0.6, 0.8),
                        legend.background = element_rect(fill = "white"))
    } else {
        if (var_to_plot == "gamma") {
            plot = plotting_df %>% ggplot(aes(x=reorder(gene, -gamma_mle), y=gamma_mle)) + 
                geom_errorbar(aes(ymin = gamma_ci_low, ymax = gamma_ci_high), width=0) +
                labs(x = "Gene", y = "Selection Intensity", color="Mutation frequency")
        } else if (var_to_plot == "flux") {
            plot = plotting_df %>% ggplot(aes(x=reorder(gene, -gamma_mle), y=flux_mle)) + 
                geom_errorbar(aes(ymin = flux_ci_low, ymax = flux_ci_high), width=0) +
                labs(x = "Gene", y = "Fixation Rate", color="Mutation frequency")
        } else if (var_to_plot == "mu") {
            plot = plotting_df %>% ggplot(aes(x=reorder(gene, -gamma_mle), y=mu)) + 
                geom_errorbar(aes(ymin = mu_ci_low, ymax = mu_ci_high), width=0) +
                labs(x = "Gene", y = "Mutation Rate", color="Mutation frequency")
        } else {stop("var_to_plot must be gamma (selection intensity), flux (mutation acquisition rate), mu (mutation rate), or freq (mutation frequency)")}

        plot = plot + 
            geom_point(aes(size=freq, color=freq)) + 
            scale_color_viridis_c() +
            scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
            theme_classic() +
            theme(
                axis.title.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 16),
                axis.text.x = element_text(angle = 90),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, hjust=0.5),
                #legend.position = c(0.9,0.8)
                legend.position = "bottom"
            ) + 
            coord_flip()
        
        if(!show_x_axis_title){plot = plot + theme(axis.title.x = element_blank())}
        if(!show_freq_legend){plot = plot + theme(legend.position="none")}
    }

    
    return(plot)
}

get_genes_with_gxe_effects = function(df, mu_method){
    gxe_effects = 
        df %>%
            filter(method == mu_method,
                    key %in% c("smoking_plus","nonsmoking_plus")) %>%
            pivot_wider(
                names_from = key,
                id_cols = c(method, gene),
                values_from  = c(gamma_mle, gamma_ci_low, gamma_ci_high)
            ) %>%
            filter((gamma_ci_low_nonsmoking_plus > gamma_ci_high_smoking_plus) | (gamma_ci_low_smoking_plus > gamma_ci_high_nonsmoking_plus)) %>%
            pull(gene)
    
    return(gxe_effects)
}

plot_GxE_results = function(df, mu_method){
    ranked_genes = df %>% filter(key == "pan_data", method==mu_method) %>% arrange(desc(gamma_mle)) %>% pull(gene)
    
    gxe_effects = get_genes_with_gxe_effects(df, mu_method)

    plotting_df = 
        df %>%
            filter(key %in% c("smoking_plus","nonsmoking_plus")) %>%
            filter(method == mu_method) %>%
            mutate(signif = ifelse(gene %in% gxe_effects, "Significant", "Not significant"),
                    signif_mark = ifelse(signif == "Significant", "*", ""),
                    gene = factor(gene, levels = ranked_genes))

    raw_gamma_plot = 
        plotting_df %>%
            mutate(key = ifelse(key == "smoking_plus", "Ever-Smokers", ifelse(key == "nonsmoking_plus","Never-Smokers",""))) %>%
            ggplot(aes(x = gene, y = gamma_mle, group=key)) +
            geom_col(aes(fill = key), alpha = 0.7, position="dodge", width = 0.75) +
            geom_errorbar(aes(ymin = gamma_ci_low, ymax = gamma_ci_high), position="dodge", width = 0.75) +
            # geom_point(data = df %>% filter(key=="pan_data", method == mu_method), aes(x=reorder(gene, -gamma_mle),y=gamma_mle)) + 
            # geom_errorbar(data = df %>% filter(key=="pan_data", method == mu_method), aes(ymin = gamma_ci_low, ymax = gamma_ci_high), width=0.5) +

            labs(x = "Gene", y = "Selection Intensity (Gamma)", title = "Average selection intensity for SNVs between ever- and never-smoker LUAD", fill = "Smoker Status") +
            theme_classic() +
            theme(
                axis.title.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 16),
                axis.text.x = element_text(angle = 90),
                axis.text = element_text(size = 14),
                legend.title = element_text(size = 15),
                legend.text = element_text(size = 13),
                plot.title = element_text(size = 20, hjust=0.5)
            ) + 
            geom_text(aes(label=signif_mark, x=gene), y=log10(0), size = 10)
    
    comparison_df = plotting_df %>%
        pivot_wider(
            names_from = key,
            id_cols = c(gene, signif),
            values_from = c(gamma_mle, flux_mle, mu, freq)
        ) %>%
        mutate(ratio = gamma_mle_smoking_plus/gamma_mle_nonsmoking_plus,
                ratio = ifelse(ratio < 1, 1/ratio, ratio)) %>%
        mutate(which_greater = ifelse(gamma_mle_smoking_plus > gamma_mle_nonsmoking_plus, "Greater in Ever-Smokers", 
                                ifelse(!(is.na(gamma_mle_smoking_plus) | is.na(gamma_mle_nonsmoking_plus)), "Greater in Never-Smokers", "null")))

    ratio_plot = 
        comparison_df %>%
        ggplot(aes(x=gene)) + 
            geom_point(aes(y=ratio, alpha=signif), size=5) + 
            scale_alpha_manual(values = c(0.5, 1)) + 
            geom_hline(yintercept=1, color="black", lty=2) +
            facet_wrap(~which_greater) + 
            labs(x = "Gene", y = "Ratio of selection intensities",title="Ratios of selection intensities (using cesR mutation rates)",alpha="Significant Difference") +
            theme_minimal_vgrid() +
            theme(
                axis.title.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 16),
                axis.text.x = element_text(angle = 90),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, hjust=0.5),
                strip.text = element_text(size = 16),
                panel.border = element_rect(colour = "black", fill = NA, size = 1)
            ) 

    
    combined_plot = plot_grid(raw_gamma_plot, ratio_plot, nrow=2, labels="AUTO", align="v", axis="l", label_size = 20)

    return(list(df = comparison_df, plot = combined_plot))
}


# M = 3

find_mutated_gene = function(mutation_indices){
    numeric_indices = as.integer(stringr::str_extract_all(mutation_indices,"\\d")[[1]])
    len = length(numeric_indices)
    # find mutated gene by finding the index where the first genotype (0,0,0) is different from the second (0,0,1)
    mutated_index = which((numeric_indices[(len/2 + 1):len] - numeric_indices[1:(len/2)]) == 1)
    return(mutated_index)
}

get_Mk_samples_per_combination = function(k){
    # consider switching to loading all availables samples per combo files in future
    pd_spc = fread(glue("{location_cesR_output}/M{k}/samples_per_combination_pan_data.csv"))
    s_plus_spc = fread(glue("{location_cesR_output}/M{k}/samples_per_combination_smoking_plus.csv"))
    ns_plus_spc = fread(glue("{location_cesR_output}/M{k}/samples_per_combination_nonsmoking_plus.csv"))

    pd_spc = pd_spc %>% 
                    pivot_longer(cols=starts_with("("),
                                    names_to = "state",
                                    values_to = "count") %>%
                    mutate(state = gsub("[()]","",state)) %>% 
                    mutate(key = "pan_data")

    s_plus_spc = s_plus_spc %>% 
                    pivot_longer(cols=starts_with("("),
                                    names_to = "state",
                                    values_to = "count") %>%
                    mutate(state = gsub("[()]","",state)) %>% 
                    mutate(key = "smoking_plus")

    ns_plus_spc = ns_plus_spc %>% 
                    pivot_longer(cols=starts_with("("),
                                    names_to = "state",
                                    values_to = "count") %>%
                    mutate(state = gsub("[()]","",state)) %>%
                    mutate(key = "nonsmoking_plus")

    samples_per_combination = bind_rows(pd_spc, s_plus_spc, ns_plus_spc)

    samples_per_combination = samples_per_combination %>%
        mutate(gene_set = gsub(" ","_",
                            gsub("['\\)\\(\\,]","",gene_combination)),
                gene_combination = NULL) 

    return(samples_per_combination)
}


load_M3_results = function(){
    samples_per_combination = get_Mk_samples_per_combination(k=3)
    lower_bound = 1e-2 # Note that the exact value of the lower bound has no practical significance, at least with these results

    gammas = fread("M3_gene_gammas.csv")
    gammas_df = gammas %>%
        # Put a lower bound on gamma because distinctions between strengths of negative selections are impractical to accurately infer
        mutate(gamma_ci_low = ifelse(gamma_ci_low < lower_bound, lower_bound, gamma_ci_low),
                gamma_ci_high = ifelse(gamma_ci_high < lower_bound, lower_bound, gamma_ci_high),
                gamma_mle = ifelse(gamma_mle < lower_bound, lower_bound, gamma_mle)) %>%

        # Information on mutated genes
        mutate(gene_set = paste0(first_gene, "_", second_gene, "_", third_gene)) %>%
        mutate(from = stringr::str_extract(mutation, "\\d, \\d, \\d(?=\\), )"),
            to = stringr::str_extract(mutation, "(?<=, \\()\\d, \\d, \\d"),
            mutated_ind = as.character(unlist(lapply(mutation, find_mutated_gene)))) %>%
        mutate(mutated_gene = ifelse(mutated_ind == 1, first_gene, ifelse(mutated_ind == 2, second_gene, third_gene))) %>%

        # Information on number of samples in `from` and `to` states
        left_join(samples_per_combination, by=c("key"="key", "gene_set"="gene_set", 
                                            "from"="state")) %>%
            dplyr::rename(from_count = count) %>%
        left_join(samples_per_combination, by=c("key"="key", "gene_set"="gene_set", 
                                            "to"="state")) %>%
        dplyr::rename(to_count = count) %>%

        # Add column for initial genotype 
        mutate(tmp_from = lapply(strsplit(from, ", "), function(x) as.logical(as.integer(x))), 
                tmp_set = strsplit(gene_set, "_"),
                from_gt = lapply(1:nrow(.), function(i) tmp_set[[i]][tmp_from[[i]]]),
                tmp_from = NULL, tmp_set = NULL) %>%
        mutate(from_gt = ifelse(from == "0, 0, 0", "WT", from_gt)) %>%
        mutate(from_gt = lapply(from_gt, function(x) x[order(x)]))

    return(gammas_df)
}

# Test for epistatic interactions from WT genotype
# will return info on all interactions, significant or not

calc_WT_epistasis_ratios = function(data){
    WT_row = data[from_gt == "WT"]
    mutant_rows = data[from_gt != "WT"]
    mutated_gene = WT_row$mutated_gene

    all_comparisons = data.table()
    for(i in 1:nrow(mutant_rows)){
        comparison = data.table(epistatic_gt = paste(sort(mutant_rows[[i, 'from_gt']]),collapse="_"),
                                mutated_gene = mutated_gene,
                                ratio = mutant_rows[[i, 'gamma_mle']] / WT_row$gamma_mle, 
                                signif = (mutant_rows[[i, 'gamma_ci_low']] > WT_row$gamma_ci_high) | (WT_row$gamma_ci_low > mutant_rows[[i, 'gamma_ci_high']]))
        all_comparisons = rbind(all_comparisons, comparison)
    }

    return(all_comparisons)
}

calc_WT_epistasis_ratios_for_all_data = function(gammas_df){
    gammas_df_list = 
        gammas_df %>% 
        select(key, gene_set, mutated_gene, from, from_gt, gamma_mle, gamma_ci_low, gamma_ci_high) %>%
        split(gammas_df %>% mutate(id = paste0(gene_set,":",mutated_gene)) %>% pull(id))

    all_epistatic_comparisons = rbindlist(lapply(gammas_df_list, calc_WT_epistasis_ratios))
    all_epistatic_comparisons = all_epistatic_comparisons %>%
        group_by(epistatic_gt, mutated_gene) %>%
        summarize(median_ratio = median(ratio), 
                    sd_ratio = sqrt(var(ratio)),
                    rel_sd = sd_ratio / median_ratio, 
                    n=n(), 
                    signif = ifelse(sum(signif) > length(signif)/2, 'significant', 'non-significant')) %>%
        ungroup() %>%
        mutate(epistasis_sign = ifelse(median_ratio > 1, "synergistic", "antagonistic"))


    return(all_epistatic_comparisons)
}

# checks for epistatic interactions from any genotype, not just WT
# only returns significant interactions

check_ci_overlap = function(data, test_gge) {
  # order from low to high
  #   this minimizes the number of comparisons that need to be checked
  #   since if a row's upper bound is less than the next upper bound,
  #   then the it is impossible for the lower bound of the row to be greater
  #   than the next row's upper bound.
  data = data[order(data$gamma_ci_high), ]
  
  non_overlaps = data.table()
  
  n = nrow(data)
  
  # iterate through each row, checking if the next rows' lower 
  # bounds are greater than the current row's upper bound
  for (i in 1:(n - 1)) {
    var1 <- data[i]
    candidates <- data[(i + 1):n, ]
    candidates <- candidates[candidates$gamma_ci_low > var1$gamma_ci_high, ]
    
    if(nrow(candidates) > 0) {
      for(j in 1:nrow(candidates)) {
        row = data.table(lower = var1$from, 
                        upper = candidates$from[j], 
                        ratio = candidates$gamma_mle[j]/var1$gamma_mle,
                        lower_gt = var1$from_gt, 
                        upper_gt = candidates$from_gt[j],
                        under_neg_sel = var1$gamma_mle < 1)
        if(test_gge) {
            row$lower_key = var1$key
            row$upper_key = candidates$key[j]
        }
        non_overlaps = rbind(non_overlaps, row)
      }
    }
  }
  if(nrow(non_overlaps) > 0) {
    non_overlaps$mutated_gene = var1$mutated_gene
    non_overlaps$gene_set = var1$gene_set
  }

  return(non_overlaps)
}

get_epistasis_results = function(gammas_df, test_gge = FALSE){
    gammas_df_list = 
        gammas_df %>% 
        select(key, gene_set, mutated_gene, from, from_gt, gamma_mle, gamma_ci_low, gamma_ci_high) %>%
        split(gammas_df %>% mutate(id = paste0(gene_set,":",mutated_gene)) %>% pull(id))

    if(test_gge){
        gammas_df_list = gammas_df_list[unlist(lapply(gammas_df_list, function(df) nrow(df) == 8))]
    }

    epistasis_results = rbindlist(lapply(gammas_df_list, function(df) check_ci_overlap(df, test_gge) ))

    return(epistasis_results)
}

summarize_epistasis_results = function(epistasis_results, exclude_synergy_with_TP53 = TRUE){
    summarized_results = epistasis_results %>%
        group_by(lower_gt, upper_gt, mutated_gene) %>%
        summarize(median_ratio = round(median(ratio),3),
                sd_ratio = round(sqrt(var(ratio)),3), 
                rel_sd = round(sqrt(var(ratio))/median(ratio),3), 
                n = n()) %>%
        ungroup()
    
    return(summarized_results)
}

# Plot epistasic interactions

plot_epistatic_interaction = function(mg, gt1, gt2, df, condensed = FALSE, edge_color="black"){
    gt1 = gt1[order(gt1)]
    gt1_collapsed = paste(gt1, collapse="_")

    gt2 = gt2[order(gt2)]
    gt2_collapsed = paste(gt2, collapse="_")

    necessary_genes = c(gt1, gt2)
    necessary_genes = necessary_genes[necessary_genes != "WT"]

    plotting_df = df %>%
        filter(mutated_gene == mg) %>%
        rowwise() %>%
        mutate(from_gt = paste(from_gt, collapse="_")) %>%
        filter(all(unlist(lapply(necessary_genes, function(gene) str_detect(gene_set,gene))))) %>%
        ungroup() %>%
        filter(from_gt %in% c(gt1_collapsed , gt2_collapsed)) %>%
        pivot_wider(names_from = from_gt, names_glue = "{.value}_{from_gt}", id_cols = gene_set, values_from = c(gamma_ci_low, gamma_mle, gamma_ci_high))


    if(condensed){
        if(gt1 == "WT"){
            xlabel = glue("SI: WT --> {mg}")
        } else {
            xlabel = glue("SI: {paste(gt1,collapse=' + ')} --> {paste(c(gt1, mg),collapse=' + ')}")
        }
        ylabel = glue("SI: {paste(gt2,collapse=' + ')} --> {paste(c(gt2, mg),collapse=' + ')}")
        title = glue("{mg} | {necessary_genes}")
    } else{
        if(gt1 == "WT"){
            xlabel = glue("Selection intensity for WT --> {mg}")
        } else {
            xlabel = glue("Selection intensity for {paste(gt1,collapse=' + ')} --> {paste(c(gt1, mg),collapse=' + ')}")
        }
        ylabel = glue("Selection intensity for {paste(gt2,collapse=' + ')} --> {paste(c(gt2, mg),collapse=' + ')}")
        title = glue("Epistatic interaction of {mg} with {paste(necessary_genes,collapse=', ')}")
    }

    p = plotting_df %>%
        ggplot(aes(x = plotting_df %>% select(starts_with("gamma_mle") & ends_with(gt1_collapsed)) %>% pull, 
                y = plotting_df %>% select(starts_with("gamma_mle") & ends_with(gt2_collapsed)) %>% pull)) + 
            geom_point(
            size = 3, shape="o", color=edge_color, alpha=1) + 
            geom_errorbar(aes(ymin = plotting_df %>% select(starts_with("gamma_ci_low") & ends_with(gt2_collapsed)) %>% pull, 
                            ymax = plotting_df %>% select(starts_with("gamma_ci_high") & ends_with(gt2_collapsed)) %>% pull), width=0, size=0.2, alpha=0.5) +
            geom_errorbarh(aes(xmin = plotting_df %>% select(starts_with("gamma_ci_low") & ends_with(gt1_collapsed)) %>% pull, 
                            xmax = plotting_df %>% select(starts_with("gamma_ci_high") & ends_with(gt1_collapsed)) %>% pull), height=0, size=0.2, alpha=0.5) +
            labs(x=xlabel,
                y=ylabel,
                title = title) +
            geom_abline(slope = 1, intercept = 0, lty=2) +
            theme_classic() + 
            theme(axis.title = element_text(size = 20), 
                    plot.title = element_text(size = 20, hjust = 0.5)) +
            geom_point(aes(x=0,y=0),color="white") + 
            scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

    return(p)
}

plot_all_epistatic_interactions_for_one_gene = function(mg, gt_base, df, condensed = FALSE, log_scale = FALSE, add_labels = TRUE){
    gt_base = gt_base[order(gt_base)]
    gt_base_collapsed = paste(gt_base, collapse="_")
    necessary_genes = gt_base[gt_base != "WT"]

    tmp_df = df %>%
    filter(mutated_gene == mg) %>%
    rowwise() %>%
    mutate(from_gt = paste(from_gt, collapse="_"),
            from_base = ifelse(from_gt == gt_base, "base", "other")) %>%
    filter(all(unlist(lapply(necessary_genes, function(gene) str_detect(gene_set,gene))))) %>%
    ungroup() %>%
    pivot_wider(names_from = from_base, names_glue = "{.value}_{from_base}", id_cols = c(gene_set, from_gt), values_from = c(gamma_ci_low, gamma_mle, gamma_ci_high))

    wt_gammas = tmp_df %>%
        filter(from_gt == gt_base_collapsed) %>% 
        select(gene_set, ends_with("base"))

    mut_gammas = tmp_df %>%
        filter(from_gt != gt_base_collapsed) %>% 
        select(gene_set, from_gt, ends_with("other"))

    plotting_df = mut_gammas %>% left_join(wt_gammas, by = "gene_set")


    if(condensed){
        if(gt_base == "WT"){
            xlabel = glue("SI: WT --> {mg}")
        } else {
            xlabel = glue("SI: {paste(gt_base,collapse=' + ')} --> {paste(c(gt_base, mg),collapse=' + ')}")
        }
        ylabel = glue("SI: [Mutant Genotype] --> [Mutant Genotype] +  {mg}")
        title = glue("{mg} | {necessary_genes}")
    } else{
        xlabel = glue("Selection intensity for {mg} mutation from {paste(gt_base,collapse=' + ')}")
        ylabel = glue("Selection intensity for {mg} mutation from Mutant Genotype")
        title = glue("Epistatic interactions of {mg} with {paste(necessary_genes,collapse=', ')}")
    }

    labels = plotting_df %>%
        group_by(from_gt) %>%
        summarize(median_ratio_base = median(gamma_mle_base), median_ratio_other = median(gamma_mle_other))

    p = plotting_df %>%
        # filter(gamma_mle_other > 1 & gamma_mle_base > 1) %>%
        ggplot(aes(x = gamma_mle_base, y = gamma_mle_other)) +
            geom_point(aes(color = from_gt), size = 3, shape = "o", alpha=1) + 
            geom_errorbar(aes(color = from_gt, ymin = gamma_ci_low_other, ymax = gamma_ci_high_other), width=0, size=0.2, alpha=0.5) +
            geom_errorbarh(aes(color = from_gt, xmin = gamma_ci_low_base, xmax = gamma_ci_high_base), height=0, size=0.2, alpha=0.5) +
            #geom_label(aes(label = from_gt, x=median(gamma_mle_base), y = median(gamma_mle_other))) +        
            labs(x = xlabel, y = ylabel, color = "Mutant Genotype") +
            geom_abline(slope = 1, intercept = 0, lty=2) +
            theme_classic() +
            theme(axis.title = element_text(size = 20), 
                        plot.title = element_text(size = 20, hjust = 0.5)) +
                geom_point(aes(x=0,y=0),color="white")
    
    if(add_labels){p = p + ggrepel::geom_label_repel(data = labels, aes(label = from_gt, x=median_ratio_base, y = median_ratio_other), size = 5, max.overlaps = 100)}
    if(log_scale){p = p + scale_x_log10() + scale_y_log10()}

    return(p)
}

pairwise_comparisons = function(data){
    WT_row = data[from_gt == "WT"]
    mutant_rows = data[from_gt != "WT"]
    mutated_gene = WT_row$mutated_gene

    all_comparisons = data.table()
    for(i in 1:nrow(mutant_rows)){
        epi_gt = paste(sort(mutant_rows[[i, 'from_gt']]),collapse="_")
        tested_combo_ = paste0(epi_gt,'_',mutated_gene)
        ratio_ = mutant_rows[[i, 'gamma_mle']] / WT_row$gamma_mle
        signif_ = (mutant_rows[[i, 'gamma_ci_low']] > WT_row$gamma_ci_high) | (WT_row$gamma_ci_low > mutant_rows[[i, 'gamma_ci_high']])

        new_WT_row = data.table(gene_set = WT_row$gene_set,
                                tested_combo = tested_combo_,
                                epistatic_gt = "WT",
                                mutated_gene = mutated_gene,
                                gamma_ci_low = WT_row$gamma_ci_low,
                                gamma_mle = WT_row$gamma_mle,
                                gamma_ci_high = WT_row$gamma_ci_high,
                                to_count = WT_row$to_count,
                                ratio = ratio_, 
                                signif = signif_)
        new_WT_row[,key:=WT_row$key]
        new_epi_row = data.table(gene_set = WT_row$gene_set,
                                tested_combo = tested_combo_,
                                epistatic_gt = epi_gt,
                                mutated_gene = mutated_gene,
                                gamma_ci_low = mutant_rows[[i, 'gamma_ci_low']],
                                gamma_mle = mutant_rows[[i, 'gamma_mle']],
                                gamma_ci_high = mutant_rows[[i, 'gamma_ci_high']],
                                to_count = mutant_rows[[i, 'to_count']],
                                ratio = ratio_, 
                                signif = signif_)
        new_epi_row[,key:=mutant_rows[[i, 'key']]]
        all_comparisons = rbind(all_comparisons, new_WT_row, new_epi_row)
    }

    return(all_comparisons)
}

plot_interactions = function(data, interactions_to_plot=NULL, custom_order=FALSE, n_interactions=10, synergy_or_antagonism = "synergy", include_higher_order=FALSE,log_scale=FALSE,spread=2/3,title="default", add_annotations=TRUE){
    gammas_df_list = 
        data %>% 
        select(key, gene_set, from_gt, mutated_gene, gamma_ci_low, gamma_mle, gamma_ci_high, to_count) %>%
        split(data %>% mutate(id = paste0(gene_set,":",mutated_gene)) %>% pull(id))
    
    interaction_df = rbindlist(lapply(gammas_df_list, pairwise_comparisons))
    if(is.null(interactions_to_plot)){
        interactions_to_plot = interaction_df %>% 
                                filter(if(!include_higher_order) str_count(tested_combo,'_')<2 else TRUE) %>% 
                                arrange(desc(ratio)) %>%
                                pull(tested_combo) %>% unique()
        if(synergy_or_antagonism %in% c("s","synergy")) {interactions_to_plot = interactions_to_plot %>% head(n_interactions)}
        else if(synergy_or_antagonism %in% c("a","antagonism")){interactions_to_plot = interactions_to_plot %>% tail(n_interactions)}
        else if(synergy_or_antagonism %in% c("b","both")){interactions_to_plot = c(interactions_to_plot %>% head(n_interactions/2),  interactions_to_plot %>% tail(n_interactions/2))}
        else{stop("`synergy_or_antagonism` can only take on character values `synergy` or `antagonism` or `both`.")}
    }                        
    interaction_df = interaction_df %>% filter(tested_combo %in% interactions_to_plot)
    interaction_df = interaction_df %>%
                        group_by(tested_combo, epistatic_gt) %>%
                        mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE),
                                combo_name = sapply(tested_combo, 
                                                function(tested_combo){tmp = str_split(tested_combo,'_')[[1]]; 
                                                                        glue("{tmp[length(tmp)]} [{tmp[-length(tmp)]}]")}),
                                combo_name = factor(combo_name, levels = sapply(interactions_to_plot, 
                                                function(combo){tmp = str_split(combo,'_')[[1]]; 
                                                                        glue("{tmp[length(tmp)]} [{tmp[-length(tmp)]}]")})),
                                median_ratio = median(ratio)) 
    interaction_df = interaction_df %>% ungroup() %>% 
        mutate(order = ifelse(epistatic_gt!="WT",median(gamma_mle),NA), .by = c(tested_combo, epistatic_gt)) %>% 
        mutate(order = ifelse(epistatic_gt=="WT",median(order,na.rm=T),order), .by=tested_combo) %>%
        rowwise() %>% mutate(other_genes = lapply(gene_set, function(x){genes = str_split(gene_set,'_')[[1]]; genes[!(genes %in% c(epistatic_gt, mutated_gene))]})) %>% ungroup() 

    alpha_palette = c("TRUE" = 1, "FALSE" = 0.1)

    if(log_scale) {interaction_df = interaction_df %>% mutate(across(starts_with("gamma"), log10))} 

    plot = interaction_df %>%
        ggplot(aes(x=if(custom_order){combo_name}else{reorder(combo_name, order)}, y=gamma_mle, size=to_count, alpha=signif)) + 
            geom_errorbar(data=interaction_df %>% filter(epistatic_gt=="WT"),
                        aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt=="WT") %>% pull(nudge_dist)),
                        width=0,
                        linewidth=0.2) + 
            geom_point(data=interaction_df %>% filter(epistatic_gt=="WT"),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt=="WT") %>% pull(nudge_dist)),
                        shape=21, color="black",fill="maroon") + 
            geom_errorbar(data=interaction_df %>% filter(epistatic_gt!="WT"),
                        aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist)),
                        width=0,
                        linewidth=0.2) + 
            geom_point(data=interaction_df %>% filter(epistatic_gt!="WT"), 
                        aes(fill=ratio),
                        shape=21, color="black",
                        position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist))) + 
            labs(x="Epistatic Pair\nMutated Gene [Genetic Context]", y=if(log_scale){expression(paste(log[10],"(Scaled Selection Coefficient)"))}else{"Scaled Selection Coefficient"}, title=if(title=="default"){"Epistatic interactions in LUAD"}else{title},
                    size="Sample count") +
            scale_alpha_manual(values = alpha_palette, name="Significant Difference in Selection") +
            scale_fill_viridis_c(name="Ratio of Selection") +
            theme_classic() +
            theme(plot.title = element_text(size = 24, hjust=0.5),
                    axis.title = element_text(size = 20),
                    axis.text = element_text(size = 16),
                    axis.ticks.x = element_blank(),
                    legend.position = c(0.8,0.2),
                    legend.direction="horizontal",
                    legend.key.size = unit(1.5, 'cm'),
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size = 20),
                    panel.grid.major.y = element_line(color="gray",linewidth=0.75, linetype=3)) +
            guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5),
                    size = guide_legend(title.position="top", title.hjust = 0.5),
                    alpha = "none")+#guide_legend(title.position="top", title.hjust = 0.5)) + 
            coord_flip()
        if(add_annotations){
            plot = plot + 
                    ggrepel::geom_text_repel(data=interaction_df %>% filter(epistatic_gt!="WT"), 
                                            aes(label=other_genes),
                                            position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist)),
                                            size = 3,
                                            max.overlaps = 20)
        }
        if(!log_scale){
            plot = plot + geom_hline(yintercept = 1, color="grey", lty=2)
        }

    return(list("plot" = plot, "df" = interaction_df))
}


plot_gge_interactions = function(data, n_interactions=10, synergy_or_antagonism = "synergy", include_higher_order=FALSE,log_scale=FALSE,spread=2/3,title="default"){
    gammas_df_list = 
        data %>% 
        select(gene_set, from_gt, mutated_gene, gamma_ci_low, gamma_mle, gamma_ci_high, to_count) %>%
        split(data %>% mutate(id = paste0(gene_set,":",mutated_gene)) %>% pull(id))
    
    interaction_df = rbindlist(lapply(gammas_df_list, pairwise_comparisons))
    interactions_to_plot = interaction_df %>% 
                            filter(if(!include_higher_order) str_count(tested_combo,'_')<2 else TRUE) %>% 
                            arrange(desc(ratio)) %>%
                            pull(tested_combo) %>% unique()
    if(synergy_or_antagonism %in% c("s","synergy")) {interactions_to_plot = interactions_to_plot %>% head(n_interactions)}
    else if(synergy_or_antagonism %in% c("a","antagonism")){interactions_to_plot = interactions_to_plot %>% tail(n_interactions)}
    else if(synergy_or_antagonism %in% c("b","both")){interactions_to_plot = c(interactions_to_plot %>% head(n_interactions/2),  interactions_to_plot %>% tail(n_interactions/2))}
    else{stop("`synergy_or_antagonism` can only take on character values `synergy` or `antagonism` or `both`.")}
                            
    interaction_df = interaction_df %>% filter(tested_combo %in% interactions_to_plot)
    interaction_df = interaction_df %>%
                        group_by(tested_combo, epistatic_gt) %>%
                        mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE),
                                combo_name = sapply(tested_combo, 
                                                function(tested_combo){tmp = str_split(tested_combo,'_')[[1]]; 
                                                                        glue("{tmp[length(tmp)]} [{tmp[-length(tmp)]}]")}),
                                median_ratio = median(ratio)) 
    interaction_df = interaction_df %>% ungroup() %>% 
        mutate(order = ifelse(epistatic_gt!="WT",median(gamma_mle),NA), .by = c(tested_combo, epistatic_gt)) %>% 
        mutate(order = ifelse(epistatic_gt=="WT",median(order,na.rm=T),order), .by=tested_combo) %>%
        rowwise() %>% mutate(other_genes = lapply(gene_set, function(x){genes = str_split(gene_set,'_')[[1]]; genes[!(genes %in% c(epistatic_gt, mutated_gene))]})) %>% ungroup() 

    alpha_palette = c("TRUE" = 1, "FALSE" = 0.1)

    if(log_scale) {interaction_df = interaction_df %>% mutate(across(starts_with("gamma"), log10))} 

    plot = interaction_df %>%
        ggplot(aes(x=reorder(combo_name, order), y=gamma_mle, size=to_count, alpha=signif)) + 
            geom_errorbar(data=interaction_df %>% filter(epistatic_gt=="WT"),
                        aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt=="WT") %>% pull(nudge_dist)),
                        width=0,
                        linewidth=0.2) + 
            geom_point(data=interaction_df %>% filter(epistatic_gt=="WT"),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt=="WT") %>% pull(nudge_dist)),
                        shape=21, color="black",fill="maroon") + 
            geom_errorbar(data=interaction_df %>% filter(epistatic_gt!="WT"),
                        aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist)),
                        width=0,
                        linewidth=0.2) + 
            geom_point(data=interaction_df %>% filter(epistatic_gt!="WT"), 
                        aes(fill=ratio),
                        shape=21, color="black",
                        position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist))) + 
            ggrepel::geom_text_repel(data=interaction_df %>% filter(epistatic_gt!="WT"), 
                        aes(label=other_genes),
                        position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist)),
                        size = 3,
                        max.overlaps = 20) +
            labs(x="Epistatic Pair\nMutated Gene [Genetic Context]", y=if(log_scale){expression(paste(log[10],"(Scaled Selection Coefficient)"))}else{"Scaled Selection Coefficient"}, title=if(title=="default"){"Epistatic interactions in LUAD"}else{title},
                    size="Sample count") +
            scale_alpha_manual(values = alpha_palette, name="Significant Difference in Selection") +
            scale_fill_viridis_c(name="Ratio of Selection") +
            theme_classic() +
            theme(plot.title = element_text(size = 24, hjust=0.5),
                    axis.title = element_text(size = 20),
                    axis.text = element_text(size = 16),
                    axis.ticks.x = element_blank(),
                    legend.position = c(0.8,0.2),
                    legend.direction="horizontal",
                    legend.key.size = unit(1.5, 'cm'),
                    legend.text = element_text(size = 16),
                    legend.title = element_text(size = 20),
                    panel.grid.major.y = element_line(color="gray",linewidth=0.75, linetype=3)) +
            guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5),
                    size = guide_legend(title.position="top", title.hjust = 0.5),
                    alpha = "none")+#guide_legend(title.position="top", title.hjust = 0.5)) + 
            coord_flip()
        if(!log_scale){
            plot = plot + geom_hline(yintercept = 1, color="grey", lty=2)
        }

    return(list("plot" = plot, "df" = interaction_df))
}



add_default_theme = function(p){
    p = p + 
        theme(axis.title = element_text(size = 18), 
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 20, hjust = 0.5))
}



# Miscellaneous functions


load_M1_all_results = function(){
    # Average fluxes and selection intensities across genetic backgrounds (M=1 model)
    M1_all_gammas = fread("M1_all_gene_gammas.csv", header = T)
    M1_all_fluxes = fread("M1_all_gene_fluxes.csv", header = T)

    # mutation frequencies for each gene (note that it doesn't matter which mutation rate method is used)
    samples_per_combo_pd = fread(glue("{location_cesR_output}/M1_all/samples_per_combination_pan_data.csv"))
    samples_per_combo_pd$key = "pan_data"
    samples_per_combo_sp = fread(glue("{location_cesR_output}/M1_all/samples_per_combination_smoking_plus.csv"))
    samples_per_combo_sp$key = "smoking_plus"
    samples_per_combo_nsp = fread(glue("{location_cesR_output}/M1_all/samples_per_combination_nonsmoking_plus.csv"))
    samples_per_combo_nsp$key = "nonsmoking_plus"
    samples_per_combo = bind_rows(samples_per_combo_pd, samples_per_combo_sp, samples_per_combo_nsp)
    colnames(samples_per_combo) = c("gene","0","1","key")
    samples_per_combo = 
        samples_per_combo %>%
        mutate(gene = gsub("[\\(\\)',]","",gene),
                freq = `1` / (`0` + `1`))

    # load mutation rates
    mus = fread("mutation_rates.csv", header = T)
    tmp = mus %>% filter(key %in% c("smoking","nonsmoking"))
    tmp = tmp %>% mutate(key = glue("{key}_plus")) 
    mus = mus %>% bind_rows(tmp)

    # combine all information (gammas, fluxes, mus, mutation frequencies)
    M1_results = 
        M1_all_gammas %>%
        left_join(M1_all_fluxes, by=c("method","key","gene")) %>%
        left_join(mus, by=c("method","key","gene")) %>%
        rename(mu = rate, mu_ci_low = rate_ci_low, mu_ci_high = rate_ci_high) %>%
        left_join(samples_per_combo %>% select(key, gene, freq), by=c("key","gene"))

    return(M1_results)
}

load_M2_results = function(){
    get_M2_samples_per_combination = function(){
        s_plus_spc = fread(glue("{location_cesR_output}/M2/samples_per_combination_smoking_plus.csv"))

        s_plus_spc = s_plus_spc %>% 
                        pivot_longer(cols=starts_with("("),
                                        names_to = "state",
                                        values_to = "count") %>%
                        mutate(state = gsub("[()]","",state)) %>% 
                        mutate(key = "smoking_plus")

        samples_per_combination = bind_rows(s_plus_spc)

        samples_per_combination = samples_per_combination %>%
            mutate(gene_set = gsub(" ","_",
                                gsub("['\\)\\(\\,]","",gene_combination)),
                    gene_combination = NULL) 

        return(samples_per_combination)
    }
    samples_per_combination = get_M2_samples_per_combination()

    gammas = fread("M2_gene_gammas.csv")
    gammas_df = gammas %>%
        # Put a lower bound on gamma because distinctions between strengths of negative selections are impractical to accurately infer
        mutate(gamma_ci_low = ifelse(gamma_ci_low < 1e-1, 1e-1, gamma_ci_low),
                gamma_ci_high = ifelse(gamma_ci_high < 1e-1, 1e-1, gamma_ci_high),
                gamma_mle = ifelse(gamma_mle < 1e-1, 1e-1, gamma_mle)) %>%

        # Information on mutated genes
        mutate(gene_set = paste0(first_gene, "_", second_gene)) %>%
        mutate(from = stringr::str_extract(mutation, "\\d, \\d(?=\\), )"),
            to = stringr::str_extract(mutation, "(?<=, \\()\\d, \\d"),
            mutated_ind = as.character(unlist(lapply(mutation, find_mutated_gene)))) %>%
        mutate(mutated_gene = ifelse(mutated_ind == 1, first_gene, second_gene)) %>%

        # Information on number of samples in `from` and `to` states
        left_join(samples_per_combination, by=c("key"="key", "gene_set"="gene_set", 
                                            "from"="state")) %>%
            dplyr::rename(from_count = count) %>%
        left_join(samples_per_combination, by=c("key"="key", "gene_set"="gene_set", 
                                            "to"="state")) %>%
        dplyr::rename(to_count = count) %>%

        # Add column for initial genotype 
        mutate(tmp_from = lapply(strsplit(from, ", "), function(x) as.logical(as.integer(x))), 
                tmp_set = strsplit(gene_set, "_"),
                from_gt = lapply(1:nrow(.), function(i) tmp_set[[i]][tmp_from[[i]]]),
                tmp_from = NULL, tmp_set = NULL) %>%
        mutate(from_gt = ifelse(from == "0, 0", "WT", from_gt)) %>%
        mutate(from_gt = lapply(from_gt, function(x) x[order(x)]))

    return(gammas_df)
}