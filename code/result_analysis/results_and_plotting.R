library(repr)

library(data.table)
library(tidyr)
library(dplyr)

library(ggplot2)
library(ggrepel)
library(ggbreak)
library(scales)
library(cowplot)

library(glue)
library(stringr)


# location_output = "../../output/"
# location_variant_output = glue("{location_output}/output_for_transfer/genes/cesR/")
# location_variant_output = glue("{location_output}/output_for_transfer/genes/variant/")
location_variant_output = glue("variant_results/")
location_cesR_output = glue("cesR_results/")

# M = 1

load_M1_results = function(mu_method){
    if(mu_method == "cesR"){location_output = location_cesR_output}
    else if (mu_method == "variant") {location_output = location_variant_output}

    # Average fluxes and selection intensities across genetic backgrounds (M=1 model)
    M1_gammas = fread(glue(location_output,'M1_gene_gammas.csv', header = T))
    M1_fluxes = fread(glue(location_output,'M1_gene_fluxes.csv', header = T))

    # mutation frequencies for each gene (note that it doesn't matter which mutation rate method is used)
    samples_per_combo = fread(glue(location_output,'M1_samples_per_combination.csv'), header=T)
    samples_per_combo = samples_per_combo %>% mutate(freq = `(1,)` / (`(0,)` + `(1,)`))

    # load mutation rates
    mus = fread(glue(location_output,'mutation_rates.csv'), header = T)

    # combine all information (gammas, fluxes, mus, mutation frequencies)
    M1_results = 
        M1_gammas %>%
        left_join(M1_fluxes, by=c("key","gene")) %>%
        left_join(mus, by=c("method","key","gene")) %>%
        left_join(samples_per_combo %>% select(key, gene, freq), by=c("key","gene"))

    return(M1_results)
}


# Plot all information
plot_M1_results = function(df, dataset_key, mu_method, var_to_plot, show_freq_legend=TRUE, show_genes=TRUE, include_mu_error=FALSE){

    plotting_df = df %>% filter(key == dataset_key, method == mu_method)

    # if(var_to_plot == "freq") {
    #     plot = plotting_df %>%
    #             ggplot(aes(x=reorder(gene, -freq), y=freq*100)) + 
    #             geom_point(aes(fill = log(gamma_mle)), pch = 21, color = "black", size = 4) + 
    #             scale_fill_viridis_c() +
    #             labs(y = "Frequency (%)", title = "Prevalence of mutations in lung adenocarcinoma", fill = "Log(Selection Intensity)") + 
    #             scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    #             theme_minimal_grid() +
                
    #             theme(axis.title.x = element_blank(), 
    #                     axis.title.y = element_text(size = 18),
    #                     axis.text.x = element_text(size = 14),
    #                     plot.title = element_text(size = 20, hjust = 0.5),
    #                     panel.grid.major.x = element_line(size=0.2),
    #                     legend.position = c(0.6, 0.8),
    #                     legend.background = element_rect(fill = "white"))
    # } else {
        if (var_to_plot == "selection") {
            plot = plotting_df %>% #mutate(gamma_mle = gamma_mle / 10^5, gamma_ci_low = gamma_ci_low / 10^5, gamma_ci_high = gamma_ci_high / 10^5) %>% 
                ggplot(aes(x=reorder(gene, gamma_mle), y=gamma_mle)) +
                geom_errorbar(aes(ymin = gamma_ci_low, ymax = gamma_ci_high), width=0) +
                labs(x = "Gene", y = "Scaled selection coefficient")
        } else if (var_to_plot == "fixation") {
            plot = plotting_df %>% ggplot(aes(x=reorder(gene, gamma_mle), y=flux_mle)) +
                geom_errorbar(aes(ymin = flux_ci_low, ymax = flux_ci_high), width=0) +
                labs(x = "Gene", y = "Fixation rate")
        } else if (var_to_plot == "mutation") {
            plot = plotting_df %>% # mutate(mu = mu * 10^6, mu_ci_low = mu_ci_low * 10^6, mu_ci_high = mu_ci_high * 10^6) %>% 
                ggplot(aes(x=reorder(gene, gamma_mle), y=mu)) +
                labs(x = "Gene", y="Mutation rate") 
            if (include_mu_error) {plot = plot + geom_errorbar(aes(ymin = mu_ci_low, ymax = mu_ci_high), width=0)}
        } else if (var_to_plot == "frequency") {
            plot = plotting_df %>% mutate(freq = 100*freq) %>% ggplot(aes(x=reorder(gene, gamma_mle), y=freq)) +
                labs(x = "Gene", y = "Prevalence (%)")
        } else {stop("var_to_plot must be selection (selection intensity), fixation (mutation acquisition rate), mutation (mutation rate), or frequency (mutation frequency)")}

        plot = plot + 
            geom_point(aes(size=freq*100, color=freq*100)) + 
            scale_color_viridis_c(labels = ~paste0(.x, "%")) +
            scale_size_continuous(labels = ~paste0(.x, "%")) + 
            scale_y_continuous(labels = ifelse(var_to_plot %in% c("fixation","frequency"), function(x)format(x, scientific=F), fancy_scientific)) +
            guides(color=guide_legend(title="Prevalence"), size = guide_legend(title="Prevalence")) +
            theme_classic() +
            theme(
                axis.title.y = element_text(size = 18),
                axis.title.x = element_text(size = 18),
                axis.ticks.y = element_blank(),
                axis.text.y = element_text(size = 16),
                #axis.text.x = element_text(angle = 90),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, hjust=0.5),
                panel.grid.major.y = element_line(color="gray",linewidth = 0.5, linetype=3),
                #legend.position = c(0.9,0.8)
                legend.position = "bottom",
            ) + 
            coord_flip()
        
        if(!show_genes){plot = plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank())}
        if(!show_freq_legend){plot = plot + theme(legend.position="none")}
    # }

    
    return(plot)
}

get_genes_with_gxe_effects = function(df, mu_method, include_panel_data){
    if(include_panel_data){keys = c("smoking_plus","nonsmoking_plus")}
    else{keys = c("smoking","nonsmoking")}

    gxe_effects = 
        df %>%
            filter(method == mu_method,
                    key %in% keys) %>%
            pivot_wider(
                names_from = key,
                id_cols = c(method, gene),
                values_from  = c(gamma_mle, gamma_ci_low, gamma_ci_high)
            ) 

    if(include_panel_data){
         gxe_effects =  gxe_effects %>%
            filter((gamma_ci_low_nonsmoking_plus > gamma_ci_high_smoking_plus) | (gamma_ci_low_smoking_plus > gamma_ci_high_nonsmoking_plus))
    } else{gxe_effects =  gxe_effects %>% filter((gamma_ci_low_nonsmoking > gamma_ci_high_smoking) | (gamma_ci_low_smoking > gamma_ci_high_nonsmoking))}
    
    gxe_effects = gxe_effects %>% pull(gene)
    
    return(gxe_effects)
}

get_smoker_nonsmoker_palette = function(){c("Ever-smoker" = hue_pal()(2)[1], "Smoker" = hue_pal()(2)[1], "Never-smoker" = hue_pal()(2)[2])}

plot_GxE_results = function(df, mu_method, ratio_plot=TRUE, include_panel_data=TRUE){
    
    gxe_effects = get_genes_with_gxe_effects(df, mu_method, include_panel_data)

    if(include_panel_data){keys = c("smoking_plus","nonsmoking_plus")}
    else{keys = c("smoking","nonsmoking")}

    plotting_df = 
        df %>%
            filter(key %in% keys) %>%
            filter(method == mu_method) %>%
            mutate(signif = ifelse(gene %in% gxe_effects, "Significant", "Not significant"),
                    signif_mark = ifelse(signif == "Significant", "*", ""))

    ranked_genes = plotting_df %>% 
        pivot_wider(
            names_from = key,
            id_cols = gene,
            values_from = gamma_mle
        ) %>%
        rowwise() %>% mutate(diff_ns = ifelse(include_panel_data, nonsmoking_plus - smoking_plus, nonsmoking-smoking), max_sel = ifelse(include_panel_data, max(nonsmoking_plus, smoking_plus), max(nonsmoking, smoking))) %>%
        arrange(diff_ns<0, ifelse(diff_ns>0, desc(max_sel), max_sel)) %>%
        pull(gene)

    plotting_df = plotting_df %>% mutate(gene = factor(gene, levels = ranked_genes))

    raw_gamma_plot = 
        plotting_df %>%
            mutate(key = ifelse(str_detect(key,"nonsmoking"), "Never-smoker", ifelse(str_detect(key,"smoking"),"Ever-smoker","")),
                    key = factor(key, levels = c("Never-smoker","Ever-smoker"))) %>%
            ggplot(aes(x = gene, y = gamma_mle, group=key)) +
            geom_col(aes(fill = key), alpha = 0.7, position="dodge", width = 0.75) +
            geom_errorbar(aes(ymin = gamma_ci_low, ymax = gamma_ci_high), position=position_dodge(width=0.75), width = 0) +
            # geom_point(data = df %>% filter(key=="pan_data", method == mu_method), aes(x=reorder(gene, -gamma_mle),y=gamma_mle)) + 
            # geom_errorbar(data = df %>% filter(key=="pan_data", method == mu_method), aes(ymin = gamma_ci_low, ymax = gamma_ci_high), width=0.5) +
            scale_fill_manual(values = get_smoker_nonsmoker_palette()) +
            scale_y_continuous(labels = fancy_scientific) +
            labs(x = "Gene", y = "Scaled selection coefficient", title = "Selection intensity for SNVs in ever-smoker and never-smoker LUAD", fill = "Smoker status") +
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

    if(!ratio_plot){return(raw_gamma_plot)}

    comparison_df = plotting_df %>%
        pivot_wider(
            names_from = key,
            id_cols = c(gene, signif),
            values_from = c(gamma_mle, flux_mle, mu, freq)
        )
    if(include_panel_data){
        comparison_df = comparison_df %>%
                        mutate(ratio = gamma_mle_smoking_plus/gamma_mle_nonsmoking_plus,
                                ratio = ifelse(ratio < 1, 1/ratio, ratio)) %>%
                        mutate(which_greater = ifelse(gamma_mle_smoking_plus > gamma_mle_nonsmoking_plus, "Greater in Ever-Smokers", 
                                                ifelse(!(is.na(gamma_mle_smoking_plus) | is.na(gamma_mle_nonsmoking_plus)), "Greater in Never-Smokers", "null")))
    } else{
        comparison_df = comparison_df %>%
                        mutate(ratio = gamma_mle_smoking/gamma_mle_nonsmoking,
                                ratio = ifelse(ratio < 1, 1/ratio, ratio)) %>%
                        mutate(which_greater = ifelse(gamma_mle_smoking > gamma_mle_nonsmoking, "Greater in Ever-Smokers", 
                                                ifelse(!(is.na(gamma_mle_smoking) | is.na(gamma_mle_nonsmoking)), "Greater in Never-Smokers", "null")))
    }

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
    samples_per_combination = fread(glue('{location_variant_output}/M{k}_samples_per_combination.csv', header = T))
    if(k==2){samples_per_combination = samples_per_combination %>% mutate(gene_set = glue("{first_gene}_{second_gene}"))}
    else if(k==3){samples_per_combination = samples_per_combination %>% mutate(gene_set = glue("{first_gene}_{second_gene}_{third_gene}"))}
    else{stop("k must be 2 or 3")}

    return(samples_per_combination)
}

# @param lower_bound: lower bound on gamma values
# Note that the exact value of the lower bound has no practical significance, at least with these results
load_M3_results = function(mu_method, lower_bound = 1e-2){
    if(mu_method == "cesR"){location_output = location_cesR_output}
    else if (mu_method == "variant") {location_output = location_variant_output}

    samples_per_combination = get_Mk_samples_per_combination(k=3)
    colnames(samples_per_combination) = trimws(colnames(samples_per_combination),"both","[\\(\\)]")
    samples_per_combination = samples_per_combination %>% 
                                select(key, gene_set, starts_with(c("0","1"))) %>%
                                pivot_longer(starts_with(c("0","1")), names_to = "state", values_to = "count")

    gammas = fread(glue(location_output,'M3_gene_gammas.csv', header = T))
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
    data = as.data.table(data)
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

get_interaction_df = function(data){
    gammas_df_list = 
        data %>% 
        select(key, gene_set, from_gt, mutated_gene, gamma_ci_low, gamma_mle, gamma_ci_high, to_count) %>%
        split(data %>% mutate(id = paste0(gene_set,":",mutated_gene)) %>% pull(id))
    interaction_df = rbindlist(lapply(gammas_df_list, pairwise_comparisons))
    interaction_df = interaction_df %>%
                        group_by(tested_combo, epistatic_gt) %>%
                        mutate(combo_name = sapply(tested_combo, 
                                                function(combo){tmp = str_split(combo,'_')[[1]]; 
                                                                        glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")}),
                                median_ratio = median(ratio))
                                
    interaction_df = interaction_df %>% ungroup() %>% 
        mutate(order = ifelse(epistatic_gt!="WT",median(gamma_mle),NA), .by = c(tested_combo, epistatic_gt)) %>% 
        mutate(order = ifelse(epistatic_gt=="WT",median(order,na.rm=T),order), .by=tested_combo) %>%
        rowwise() %>% mutate(other_genes = lapply(gene_set, 
                        function(x){genes = str_split(gene_set,'_')[[1]]; 
                                    genes[!(genes %in% c(str_split(epistatic_gt,'_')[[1]], mutated_gene))]})) %>% ungroup()                                
    return(interaction_df)
}

process_interaction_df = function(df, interactions_to_plot=NULL, n_interactions=10, synergy_or_antagonism = "both", include_higher_order=FALSE, spread=2/3){
    interaction_df = df
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
                                combo_name = factor(combo_name, levels = sapply(interactions_to_plot, 
                                                function(combo){tmp = str_split(combo,'_')[[1]]; 
                                                                        glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")})))
    return(interaction_df)
}

plot_interactions_from_df = function(df, custom_order=FALSE, log_scale=FALSE,title="default", add_annotations=TRUE, minimal_labels=TRUE){
    alpha_palette = c("TRUE" = 1, "FALSE" = 0.1)
    interaction_df = df

    if(log_scale) {interaction_df = interaction_df %>% mutate(across(starts_with("gamma"), log10))} 

    max_gamma = interaction_df %>% pull(gamma_ci_high) %>% max(.,na.rm=T)

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
            labs(x="Epistatic pair\nMutated gene [genetic context]", y=if(log_scale){expression(paste(log[10],"(Scaled selection coefficient)"))}else{"Scaled selection coefficient"}, title=if(title=="default"){"Epistatic interactions in LUAD"}else{title},
                    size="Sample count") +
            scale_alpha_manual(values = alpha_palette, name="Significant Difference in Selection") +
            scale_fill_viridis_c(name="Ratio of Selection") +
            scale_y_continuous(labels = fancy_scientific) +
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
                                            size = 5,
                                            max.overlaps = ifelse(minimal_labels,6,20))
        }
        if(!log_scale){
            plot = plot + geom_hline(yintercept = 1, color="grey", lty=2)
        }

    return(plot)
}

plot_interactions = function(df, interactions_to_plot=NULL, n_interactions=10, synergy_or_antagonism = "both", include_higher_order=FALSE, custom_order=FALSE, log_scale=FALSE,spread=2/3,title="default", add_annotations=TRUE, minimal_labels=TRUE){
    processed_interaction_df = process_interaction_df(df, interactions_to_plot, n_interactions, synergy_or_antagonism, include_higher_order, spread)
    plot_interactions_from_df(processed_interaction_df, custom_order, log_scale, title, add_annotations, minimal_labels)
}

plot_interactions_by_env = function(env1_df, env2_df, interactions_to_plot, labels=c("Environment 1","Environment 2"), log_scale=FALSE,spread=2/3, add_annotations=TRUE, minimal_labels=TRUE){
    env1_df = env1_df %>% filter(tested_combo %in% interactions_to_plot)
    env2_df = env2_df %>% filter(tested_combo %in% interactions_to_plot)
    env1_missing = setdiff(env2_df %>% pull(tested_combo), env1_df %>% pull(tested_combo))
    env2_missing = setdiff(env1_df %>% pull(tested_combo), env2_df %>% pull(tested_combo))

    if(length(env1_missing) > 0){
        env1_rows_to_add = data.table(combo_name = sapply(env1_missing, function(combo){tmp = str_split(combo,'_')[[1]]; 
                                               glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")}),
                                 epistatic_gt = rep("WT",length(env1_missing)))
        env1_df = env1_df %>% bind_rows(env1_rows_to_add)
    }
    if(length(env2_missing) > 0){
        env2_rows_to_add = data.table(combo_name = sapply(env2_missing, function(combo){tmp = str_split(combo,'_')[[1]]; 
                                               glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")}),
                                 epistatic_gt = rep("WT",length(env2_missing)))
        env2_df = env2_df %>% bind_rows(env2_rows_to_add)
    }

    interactions_to_plot = env1_df %>% 
        filter(ifelse(median_ratio>1, epistatic_gt != "WT", epistatic_gt == "WT")) %>%
        group_by(tested_combo) %>%
        summarize(mean_max_gamma = mean(gamma_mle), median_ratio = median(median_ratio)) %>%
        arrange(median_ratio<1, ifelse(median_ratio>1,desc(mean_max_gamma), mean_max_gamma)) %>%
        pull(tested_combo) %>%
        rev

    env1_df = env1_df %>%
                group_by(tested_combo, epistatic_gt) %>%
                mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE),
                        combo_name = factor(combo_name, levels = sapply(interactions_to_plot, 
                                        function(combo){tmp = str_split(combo,'_')[[1]]; 
                                                                glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")}))) %>%
                ungroup()
    env2_df = env2_df %>%
                group_by(tested_combo, epistatic_gt) %>%
                mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE),
                        combo_name = factor(combo_name, levels = sapply(interactions_to_plot, 
                                        function(combo){tmp = str_split(combo,'_')[[1]]; 
                                                                glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")}))) %>%
                ungroup()
    
    fill_range = range(c(env1_df %>% pull(ratio), env2_df %>% pull(ratio)), na.rm=T)
    fill_scale = scale_fill_viridis_c(limits = fill_range, name="Ratio of Selection")
    # size_range = range(c(env1_df %>% pull(to_count), env2_df %>% pull(to_count)), na.rm=T)
    # size_scale = scale_size(range = size_range/c(size_range[1], size_range[2]/6))

    # to_count_table = bind_rows(env1_df %>% mutate(key="env1") %>% select(key, to_count), env2_df %>% mutate(key="env2") %>% select(key, to_count)) %>% 
    #                 mutate(to_count = rescale(to_count, size_range))
    # env1_df = env1_df %>% mutate(to_count = to_count_table %>% filter(key=="env1") %>% pull(to_count))
    # env2_df = env2_df %>% mutate(to_count = to_count_table %>% filter(key=="env2") %>% pull(to_count))

    gamma_ceiling_env1 = (env1_df %>% pull(gamma_ci_high) %>% max(.,na.rm=T) / 1e5) %>% ceiling * 1e5
    gamma_mle_ceiling_env1 = (env1_df %>% pull(gamma_mle) %>% max(.,na.rm=T) / 1e5) %>% ceiling * 1e5
    gamma_ceiling_env2 = (env2_df %>% pull(gamma_ci_high) %>% max(.,na.rm=T) / 1e5) %>% ceiling * 1e5
    gamma_mle_ceiling_env2 = (env2_df %>% pull(gamma_mle) %>% max(.,na.rm=T) / 1e5) %>% ceiling * 1e5

    # TODO: let order be custom, env1, or env2
    env1_plot = plot_interactions_from_df(env1_df, custom_order=TRUE, log_scale, title=labels[1], add_annotations, minimal_labels) +
                    fill_scale + 
                    guides(fill = "none") + 
                    theme(axis.text.y=element_text(size=18, hjust=0.5), axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
                        axis.ticks.x = element_line(colour = "grey50", linewidth=2)) +
                    scale_y_continuous(labels = fancy_scientific, limits = c(0, gamma_mle_ceiling_env2), breaks = seq(0, gamma_mle_ceiling_env2, by = 1e5), 
                                expand=expand_scale(mult = c(0, 0),add = c(0, 2e4)))
                    
    env1_legend = get_legend(env1_plot)
    env1_plot = env1_plot + theme(legend.position="none")

    env2_plot = plot_interactions_from_df(env2_df, custom_order=TRUE, log_scale, title=labels[2], add_annotations, minimal_labels) +
                    fill_scale + 
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                        legend.box = "horizontal") + 
                    scale_y_break(c(gamma_mle_ceiling_env2 + 1e5, gamma_ceiling_env2-2e5), space=2, expand=FALSE) +
                    scale_y_reverse(labels = fancy_scientific, limits = c (gamma_ceiling_env2,0), breaks = seq(0, gamma_ceiling_env2, by = 1e5),
                                    expand = c(0, 0)) + 
                    theme(
                        axis.text.x.top = element_blank(),
                        axis.line.x.top = element_blank(),
                        axis.ticks.x.top = element_blank(),
                        axis.ticks.x.bottom = element_line(colour = "grey50", linewidth=2)
                    )

    env2_legend = get_legend(env2_plot)
    env2_plot = env2_plot + theme(legend.position="none")

    aplot::plot_list(env2_plot, env1_plot, env1_legend, env2_legend, nrow=2, heights = c(1,0.05))
}


plot_gge_interaction = function(mg, gt, smoking_df, nonsmoking_df){
    smoking_df %>% bind_rows(nonsmoking_df, .id="key") %>% filter(tested_combo == paste0(gt,"_",mg)) %>%
    mutate(key = ifelse(key==1, "Smoker","Never-smoker")) %>%
    mutate(nudge_dist = scale((1:n())/n()**(1/.6), scale=FALSE), .by = c(key, epistatic_gt),
            x_label = ifelse(epistatic_gt == "WT", paste0(gt,"\nWT"), paste0(gt,"-\nmutant")),
            x_label = ifelse(key == "Smoker", paste0(x_label," "), x_label),
            x_label = factor(x_label, levels = c(paste0(gt,"\nWT "), paste0(gt,"-\nmutant "),paste0(gt,"\nWT"), paste0(gt,"-\nmutant")))) %>% {#levels=c("Ever-smoker WT",paste0("Ever-smoker ",gt),"Never-smoker WT", paste0("Never-smoker ",gt)))) %>% {
        ggplot(.,aes(x=x_label, y = gamma_mle)) + 
            geom_errorbar(aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                            width=0,linewidth=0.3,
                            position=position_nudge(.$nudge_dist)) + 
            geom_point(aes(fill=key, size=to_count),
                            color="black",
                            shape=21, 
                            position=position_nudge(.$nudge_dist)) + 
            geom_vline(xintercept = 2.5, lty=2,col="grey") +
            scale_size_continuous(range=c(2,8)) +
            scale_fill_manual(values = get_smoker_nonsmoker_palette()) +
            scale_y_continuous(labels=function(x)x/1e5) +
            # scale_size_continuous() +
            #scale_size_binned(n.breaks = 2) +
            labs(y="Scaled selection coefficient", title= paste("Selection for",mg,"mutations"), size="Sample count") +
            guides(size=guide_legend(nrow=1), fill=guide_legend(override.aes = list(size=6))) +
            theme_classic() +
            theme(axis.title.x = element_blank(),
                plot.title = element_text(size = 20, hjust=0.5),
                axis.text = element_text(size = 16), 
                axis.title.y = element_text(size = 20),
                legend.text = element_text(size=14))
    }
}

# With modification from https://groups.google.com/g/ggplot2/c/a_xhMoQyxZ4
fancy_scientific = function(l) {
    l = format(l, scientific = TRUE)
    l = gsub("^0(.*)e(.*)", "0", l) # Just show 0 when is 0e...
    l = gsub("^(.*)e", "'\\1'e", l)
    l = gsub("\\+","",l) # modification to remove unnecessary pluses
    l = gsub("e", "%*%10^", l)
    return(parse(text=l))
}

log_labels = function(l) {
    l = format(l, scientific = TRUE)
    l = gsub("\\+","",l)
    l = gsub("^1e", "10^", l)
    return(parse(text=l))
}

# plot_interactions = function(df, interactions_to_plot=NULL, custom_order=FALSE, log_scale=FALSE,spread=2/3,title="default", add_annotations=TRUE, minimal_labels=TRUE){
#     interaction_df = df
#     if(is.null(interactions_to_plot)){
#         interactions_to_plot = interaction_df %>% 
#                                 filter(if(!include_higher_order) str_count(tested_combo,'_')<2 else TRUE) %>% 
#                                 arrange(desc(ratio)) %>%
#                                 pull(tested_combo) %>% unique()
#         if(synergy_or_antagonism %in% c("s","synergy")) {interactions_to_plot = interactions_to_plot %>% head(n_interactions)}
#         else if(synergy_or_antagonism %in% c("a","antagonism")){interactions_to_plot = interactions_to_plot %>% tail(n_interactions)}
#         else if(synergy_or_antagonism %in% c("b","both")){interactions_to_plot = c(interactions_to_plot %>% head(n_interactions/2),  interactions_to_plot %>% tail(n_interactions/2))}
#         else{stop("`synergy_or_antagonism` can only take on character values `synergy` or `antagonism` or `both`.")}
#     }                        
#     interaction_df = interaction_df %>% filter(tested_combo %in% interactions_to_plot)
#     interaction_df = interaction_df %>%
#                         group_by(tested_combo, epistatic_gt) %>%
#                         mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE),
#                                 combo_name = factor(combo_name, levels = sapply(interactions_to_plot, 
#                                                 function(combo){tmp = str_split(combo,'_')[[1]]; 
#                                                                         glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")})))

#     alpha_palette = c("TRUE" = 1, "FALSE" = 0.1)

#     if(log_scale) {interaction_df = interaction_df %>% mutate(across(starts_with("gamma"), log10))} 

#     plot = interaction_df %>%
#         ggplot(aes(x=if(custom_order){combo_name}else{reorder(combo_name, order)}, y=gamma_mle, size=to_count, alpha=signif)) + 
#             geom_errorbar(data=interaction_df %>% filter(epistatic_gt=="WT"),
#                         aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
#                         position=position_nudge(interaction_df %>% filter(epistatic_gt=="WT") %>% pull(nudge_dist)),
#                         width=0,
#                         linewidth=0.2) + 
#             geom_point(data=interaction_df %>% filter(epistatic_gt=="WT"),
#                         position=position_nudge(interaction_df %>% filter(epistatic_gt=="WT") %>% pull(nudge_dist)),
#                         shape=21, color="black",fill="maroon") + 
#             geom_errorbar(data=interaction_df %>% filter(epistatic_gt!="WT"),
#                         aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
#                         position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist)),
#                         width=0,
#                         linewidth=0.2) + 
#             geom_point(data=interaction_df %>% filter(epistatic_gt!="WT"), 
#                         aes(fill=ratio),
#                         shape=21, color="black",
#                         position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist))) + 
#             labs(x="Epistatic Pair\nMutated Gene [Genetic Context]", y=if(log_scale){expression(paste(log[10],"(Scaled Selection Coefficient)"))}else{"Scaled Selection Coefficient"}, title=if(title=="default"){"Epistatic interactions in LUAD"}else{title},
#                     size="Sample count") +
#             scale_alpha_manual(values = alpha_palette, name="Significant Difference in Selection") +
#             scale_fill_viridis_c(name="Ratio of Selection") +
#             theme_classic() +
#             theme(plot.title = element_text(size = 24, hjust=0.5),
#                     axis.title = element_text(size = 20),
#                     axis.text = element_text(size = 16),
#                     axis.ticks.x = element_blank(),
#                     legend.position = c(0.8,0.2),
#                     legend.direction="horizontal",
#                     legend.key.size = unit(1.5, 'cm'),
#                     legend.text = element_text(size = 16),
#                     legend.title = element_text(size = 20),
#                     panel.grid.major.y = element_line(color="gray",linewidth=0.75, linetype=3)) +
#             guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5),
#                     size = guide_legend(title.position="top", title.hjust = 0.5),
#                     alpha = "none")+#guide_legend(title.position="top", title.hjust = 0.5)) + 
#             coord_flip()
#         if(add_annotations){
#             plot = plot + 
#                     ggrepel::geom_text_repel(data=interaction_df %>% filter(epistatic_gt!="WT"), 
#                                             aes(label=other_genes),
#                                             position=position_nudge(interaction_df %>% filter(epistatic_gt!="WT") %>% pull(nudge_dist)),
#                                             size = 3,
#                                             max.overlaps = ifelse(minimal_labels,6,20))
#         }
#         if(!log_scale){
#             plot = plot + geom_hline(yintercept = 1, color="grey", lty=2)
#         }

#     return(plot)
# }

add_default_theme = function(p){
    p = p + 
        theme(axis.title = element_text(size = 18), 
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 20, hjust = 0.5))
}

sym_diff = function(s1,s2){unique(c(setdiff(s1,s2), setdiff(s2,s1)))}

t_test = function(mu1, mu2, sd1, sd2, n1, n2, mu0 = 0, detailed=F){
    t_stat = (mu1 - mu2 - mu0) / sqrt((sd1^2/n1) + (sd2^2/n2))
    df = ((sd1^2/n1 + sd2^2/n2)^2) / ((sd1^4/(n1^2*(n1 - 1))) + (sd2^4/(n2^2*(n2 - 1))))
    p_val = 2 * pt(-abs(t_stat), df)

    if(!detailed){
        return(p_val)
    } else{
        return(list(t_stat, p_val))
    }
}

head_and_tail = function(df, n=6, n_head=n, n_tail=n){rbind(head(df, n_head),tail(df, n_tail))}


# Miscellaneous functions

load_M1_all_results = function(mu_method){
    if(mu_method == "cesR"){location_output = location_cesR_output}
    else if (mu_method == "variant") {location_output = location_variant_output}
    # Average fluxes and selection intensities across genetic backgrounds (M=1 model)
    M1_gammas = fread(glue(location_output,'M1_all_gene_gammas.csv', header = T))
    M1_fluxes = fread(glue(location_output,'M1_all_gene_fluxes.csv', header = T))

    # mutation frequencies for each gene (note that it doesn't matter which mutation rate method is used)
    samples_per_combo = fread(glue(location_output,'M1_all_samples_per_combination.csv'), header=T)
    samples_per_combo = samples_per_combo %>% mutate(freq = `(1,)` / (`(0,)` + `(1,)`))

    # load mutation rates
    mus = fread(glue(location_output,'mutation_rates.csv'), header = T)

    # combine all information (gammas, fluxes, mus, mutation frequencies)
    M1_results = 
        M1_gammas %>%
        left_join(M1_fluxes, by=c("key","gene")) %>%
        left_join(mus, by=c("method","key","gene")) %>%
        left_join(samples_per_combo %>% select(key, gene, freq), by=c("key","gene"))

    return(M1_results)
}

load_M2_results = function(mu_method, lower_bound = 1e-2){
    if(mu_method == "cesR"){location_output = location_cesR_output}
    else if (mu_method == "variant") {location_output = location_variant_output}

    samples_per_combination = get_Mk_samples_per_combination(k=2)
    colnames(samples_per_combination) = trimws(colnames(samples_per_combination),"both","[\\(\\)]")
    samples_per_combination = samples_per_combination %>% 
                                select(key, gene_set, starts_with(c("0","1"))) %>%
                                pivot_longer(starts_with(c("0","1")), names_to = "state", values_to = "count")

    gammas = fread(glue(location_output,'M2_gene_gammas.csv', header = T))
    gammas_df = gammas %>%
        # Put a lower bound on gamma because distinctions between strengths of negative selections are impractical to accurately infer
        mutate(gamma_ci_low = ifelse(gamma_ci_low < lower_bound, lower_bound, gamma_ci_low),
                gamma_ci_high = ifelse(gamma_ci_high < lower_bound, lower_bound, gamma_ci_high),
                gamma_mle = ifelse(gamma_mle < lower_bound, lower_bound, gamma_mle)) %>%

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