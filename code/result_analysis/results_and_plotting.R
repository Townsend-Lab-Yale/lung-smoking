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


location_variant_output = "variant_results/"
location_cesR_output = "cesR_results/"

subset_genes = c("TP53",    "KRAS",  "EGFR",  "BRAF",   "CTNNB1",
                "KEAP1",   "STK11", "ATM",   "PIK3CA", "RBM10",
                "SMARCA4", "SMAD4", "ALK",   "ARID1A", "APC",
                "MET",     "RB1",   "SETD2", "BRCA2",  "MGA",
                "GNAS")

smoking_comparison_keys = c("smoking_plus","nonsmoking_plus")

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
            scale_color_viridis_c(labels = ~paste0(.x, "%"), breaks=c(5,10,15,20,30,40), limits=c(0,41)) +
            scale_size_continuous(labels = ~paste0(.x, "%"), breaks=c(5,10,15,20,30,40), limits=c(0,41)) + 
            #scale_y_continuous(labels = ifelse(var_to_plot %in% c("fixation","frequency"), function(x)format(x, scientific=F), scientific_expr)) +
            guides(color=guide_legend(title="Prevalence", nrow=1), size = guide_legend(title="Prevalence",nrow=1)) +
            theme_classic() +
            theme(
                axis.title.y = element_text(size = 16),
                axis.title.x = element_text(size = 16),
                axis.ticks.y = element_blank(),
                axis.text.y = element_text(size = 14, face="italic"),
                #axis.text.x = element_text(angle = 90),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, hjust=0.5),
                panel.grid.major.y = element_line(color="gray",linewidth = 0.5, linetype=3),
                #legend.position = c(0.9,0.8)
                legend.position = "right",
                legend.title = element_text(hjust=0.5, vjust=-0.5)
            ) + 
            coord_flip()
        
        if(!show_genes){plot = plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank())}
        if(!show_freq_legend){plot = plot + theme(legend.position="none")}
    # }

    
    return(plot)
}

# M = 1 GxE effects
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

get_smoker_nonsmoker_palette = function(){c("Ever-smoker" = hue_pal()(2)[1], "Smoker" = hue_pal()(2)[1], "smoking_plus" = hue_pal()(2)[1], "smoking" = hue_pal()(2)[1], 
                                            "Never-smoker" = hue_pal()(2)[2], "nonsmoking_plus" = hue_pal()(2)[2], "nonsmoking" = hue_pal()(2)[2])}

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
            scale_y_continuous(labels = scientific_expr) +
            labs(x = "Gene", y = "Scaled selection coefficient", title = "Selection intensity for SNVs in ever-smoker and never-smoker LUAD", fill = "Smoker status") +
            theme_classic() +
            theme(
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 16),
                axis.text.x = element_text(angle = 45, hjust=0.95, face="italic"),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 14),
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
                axis.text.x = element_text(angle = 90, face="italic"),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 18, hjust=0.5),
                strip.text = element_text(size = 16),
                panel.border = element_rect(colour = "black", fill = NA, size = 1)
            ) 

    
    combined_plot = plot_grid(raw_gamma_plot, ratio_plot, nrow=2, labels="AUTO", align="v", axis="l", label_size = 20)

    return(list(df = comparison_df, plot = combined_plot))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# M = 3

# Load results
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Test for epistatic interactions from WT genotype
# will return info on all interactions, significant or not
# Only used for epistasis heatmaps and networks

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Scatterplot representation of epistatic interactions, showing selection from one genotype on one axis 
# and selection from the other genotype on the other axis

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# G x G (Pairwise)
## Both pairwise and higher-order epistasis analyses begin with getting the interaction dataframe

# Welch's t-test for whether epistatic interactions are significant, assuming normality of 
# scaled selection coefficient estimates
epistasis_t_test = function(wt_gamma_mle, mut_gamma_mle, wt_gamma_ci_low, wt_gamma_ci_high, 
                            mut_gamma_ci_low, mut_gamma_ci_high, wt_n, mut_n){
    # Assigning standard deviation values based on direction of comparison due to asymmetric 
    # bounds
    if(mut_gamma_mle > wt_gamma_mle){
        mut_sd = (mut_gamma_mle - mut_gamma_ci_low)
        wt_sd = (wt_gamma_ci_high - wt_gamma_mle)
    } else {
        mut_sd = (mut_gamma_ci_high - mut_gamma_mle)
        wt_sd = (wt_gamma_mle - wt_gamma_ci_low)}
    
    z_stat=qnorm(0.975) ## z score value for 95% confidence interval
    wt_sd = (wt_gamma_ci_high - wt_gamma_mle)*sqrt(wt_n)/z_stat
    mut_sd = (mut_gamma_ci_high - mut_gamma_mle)*sqrt(mut_n)/z_stat

    # print(paste0("WT std dev: ",wt_sd))
    # print(paste0("Mut std dev: ",mut_sd))

    return(t_test(wt_gamma_mle, mut_gamma_mle, wt_sd, mut_sd, wt_n, mut_n))
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
        # Assigning p-values via Welch's t-test
        # p_val_ = epistasis_t_test(WT_row$gamma_mle, mutant_rows[[i, 'gamma_mle']], WT_row$gamma_ci_low, WT_row$gamma_ci_high, mutant_rows[[i, 'gamma_ci_low']], mutant_rows[[i, 'gamma_ci_high']], WT_row$from_count, mutant_rows[[i, 'from_count']])

        new_WT_row = data.table(gene_set = WT_row$gene_set,
                                tested_combo = tested_combo_,
                                epistatic_gt = "WT",
                                mutated_gene = mutated_gene,
                                gamma_ci_low = WT_row$gamma_ci_low,
                                gamma_mle = WT_row$gamma_mle,
                                gamma_ci_high = WT_row$gamma_ci_high,
                                from_count = WT_row$from_count,
                                to_count = WT_row$to_count,
                                ratio = ratio_, 
                                signif = signif_)
                                #p_val = p_val_)
        new_WT_row[,key:=WT_row$key]
        new_epi_row = data.table(gene_set = WT_row$gene_set,
                                tested_combo = tested_combo_,
                                epistatic_gt = epi_gt,
                                mutated_gene = mutated_gene,
                                gamma_ci_low = mutant_rows[[i, 'gamma_ci_low']],
                                gamma_mle = mutant_rows[[i, 'gamma_mle']],
                                gamma_ci_high = mutant_rows[[i, 'gamma_ci_high']],
                                from_count = mutant_rows[[i, 'from_count']],
                                to_count = mutant_rows[[i, 'to_count']],
                                ratio = ratio_, 
                                signif = signif_)
                                #p_val = p_val_)
        new_epi_row[,key:=mutant_rows[[i, 'key']]]
        all_comparisons = rbind(all_comparisons, new_WT_row, new_epi_row)
    }

    return(all_comparisons)
}

get_interaction_df = function(data){
    gammas_df_list = 
        data %>% 
        select(key, gene_set, from_gt, mutated_gene, gamma_ci_low, gamma_mle, gamma_ci_high, from_count, to_count) %>%
        split(data %>% mutate(id = paste0(gene_set,":",mutated_gene)) %>% pull(id))
    rbindlist(lapply(gammas_df_list, pairwise_comparisons)) %>%
            group_by(tested_combo, epistatic_gt) %>%
            mutate(combo_name = sapply(tested_combo, 
                                    function(combo){tmp = str_split(combo,'_')[[1]]; 
                                                            glue("{tmp[length(tmp)]} [{paste(tmp[-length(tmp)],collapse='+')}]")})) %>% 
            ungroup() %>% 
            mutate(order = ifelse(epistatic_gt!="WT",median(gamma_mle),NA), .by = c(tested_combo, epistatic_gt)) %>% 
            mutate(order = ifelse(epistatic_gt=="WT",median(order,na.rm=T),order), .by=tested_combo) %>%
            rowwise() %>%
            mutate(other_genes = ifelse(epistatic_gt == "WT",str_remove(gene_set,mutated_gene),
                                        ifelse(str_detect(epistatic_gt,'_'),'',
                                            str_remove_all(gene_set,str_c(c(mutated_gene,epistatic_gt),collapse='|'))))) %>%
            ungroup() %>% 
            mutate(other_genes = gsub('(__|_)','+',trimws(other_genes,which="both",whitespace='_')))
}

# Plotting selection of synergistic and/or antagonistic pairwise interactions

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
            scale_y_continuous(labels = scientific_expr) +
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

plot_all_pairwise_epistatic_effects = function(interaction_df, baseline_selection_df, genes=NULL,ceiling=NULL){
    alpha_palette = c("TRUE" = 1, "FALSE" = 0.2)
    fill_palette = c("WT" = "red", "Mutant genotype\n([gene] mutated)" = "darkgray")

    tmp = interaction_df %>%
        filter(signif) %>%
        filter(if(!is.null(genes)) {mutated_gene %in% genes} else TRUE) %>%
        group_by(mutated_gene, epistatic_gt == "WT") %>%
            arrange(gamma_mle, .by_group=TRUE) %>%
            mutate(nudge_dist = scale((1:n())/n()**(1/0.8), scale=FALSE),
                        label_nudge_dist = ifelse(nudge_dist>0,nudge_dist+0.15,nudge_dist-0.15))  %>%
        ungroup()

    tmp = tmp %>% rowwise() %>% mutate(label_y_position = ifelse(gamma_mle<1,0,min(gamma_ci_high+5e4,ceiling))) %>% ungroup()

    tmp = tmp %>%
        filter(epistatic_gt != "WT") %>%
        bind_rows(baseline_selection_df %>% 
            filter(if(!is.null(genes)) {gene %in% genes} else TRUE) %>%
            mutate(mutated_gene=gene, epistatic_gt="WT", signif=TRUE, nudge_dist=0) %>% 
            left_join(tmp %>% 
                        filter(epistatic_gt == "WT") %>% 
                        group_by(mutated_gene) %>% 
                        summarize(from_count = mean(from_count), to_count=mean(to_count)), by = "mutated_gene")) %>% 
        mutate(context = ifelse(epistatic_gt == "WT","WT","Mutant genotype\n([gene] mutated)"),
                color = ifelse(epistatic_gt=="WT","WT",ifelse(ratio>1,"Synergism","Antagonism")))

    gamma_order = tmp %>% filter(context=="WT") %>% arrange(gamma_mle) %>% pull(mutated_gene)
    tmp = tmp %>% mutate(mutated_gene = factor(mutated_gene, levels = gamma_order))

    if(!is.null(ceiling)){
        tmp = tmp %>% mutate(gamma_ci_high = ifelse(gamma_ci_high > ceiling, ceiling, gamma_ci_high))
    }

    plot = tmp %>%
        ggplot(aes(x=mutated_gene, y=gamma_mle, alpha=signif)) + 
                geom_errorbar(
                            aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                            width=0,
                            linewidth=0.2,
                            position=position_nudge(tmp %>% pull(nudge_dist))) + 
                geom_point(
                            aes(fill=context),
                            size=3,
                            shape=21, color="black",
                            position=position_nudge(tmp %>% pull(nudge_dist))) + 
                geom_text(data=tmp %>% filter(context!="WT", signif), 
                                                # aes(label=epistatic_gt),
                                                # position=position_nudge(tmp %>% filter(context!="WT", signif) %>% pull(label_nudge_dist)),
                                                aes(label=paste0('[',epistatic_gt,']'), y=0.1),
                                                position=position_nudge(x=tmp %>% filter(context!="WT", signif) %>% pull(nudge_dist), 
                                                                        y=tmp %>% filter(context!="WT", signif) %>% pull(label_y_position)),
                                                size = 5,
                                                check_overlap = TRUE,
                                                fontface="italic") +
                labs(x="Mutation under selection", y="Scaled selection coefficient for mutation in somatic genotype",
                        size="Number of samples\nwith both mutations") +
                scale_alpha_manual(values = alpha_palette, name="Significant difference in selection") +
                scale_fill_manual(values = fill_palette, name = "Genetic context") +
                scale_y_continuous(labels = scientific_expr) +
                theme_classic() +
                theme(plot.title = element_text(size = 24, hjust=0.5),
                        axis.title = element_text(size = 24),
                        axis.text.y = element_text(size = 20, face="italic"),
                        axis.text.x = element_text(size = 20),
                        legend.position.inside = c(0.8,0.2),
                        legend.key.size = unit(1.5, 'cm'),
                        legend.text = element_text(size = 20),
                        legend.title = element_text(size = 20),
                        panel.grid.major.y = element_line(color="gray",linewidth=0.75, linetype=3)) +
                guides(fill = guide_legend(title.position="top", title.hjust = 0.5, order=1, override.aes = list(size=6)),
                        size = guide_legend(title.position="top", title.hjust = 0.5, order=2),
                        alpha = "none")+#guide_legend(title.position="top", title.hjust = 0.5)) + 
                coord_flip()
    return(plot)
}

plot_interactions = function(df, interactions_to_plot=NULL, n_interactions=10, synergy_or_antagonism = "both", include_higher_order=FALSE, custom_order=FALSE, log_scale=FALSE,spread=2/3,title="default", add_annotations=TRUE, minimal_labels=TRUE){
    processed_interaction_df = process_interaction_df(df, interactions_to_plot, n_interactions, synergy_or_antagonism, include_higher_order, spread)
    plot_interactions_from_df(processed_interaction_df, custom_order, log_scale, title, add_annotations, minimal_labels)
}

# Plotting multiple comparisons of pairwise epistatic interactions between two environments

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
                    scale_y_continuous(labels = scientific_expr, limits = c(0, gamma_mle_ceiling_env2), breaks = seq(0, gamma_mle_ceiling_env2, by = 1e5), 
                                expand=expand_scale(mult = c(0, 0),add = c(0, 2e4)))
                    
    env1_legend = get_legend(env1_plot)
    env1_plot = env1_plot + theme(legend.position="none")

    env2_plot = plot_interactions_from_df(env2_df, custom_order=TRUE, log_scale, title=labels[2], add_annotations, minimal_labels) +
                    fill_scale + 
                    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                        legend.box = "horizontal") + 
                    scale_y_break(c(gamma_mle_ceiling_env2 + 1e5, gamma_ceiling_env2-2e5), space=2, expand=FALSE) +
                    scale_y_reverse(labels = scientific_expr, limits = c (gamma_ceiling_env2,0), breaks = seq(0, gamma_ceiling_env2, by = 1e5),
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

# Plotting singular comparisons of pairwise epistasis between two environments

plot_gge_interaction = function(mg, gt, smoking_df, nonsmoking_df){
    y_scale_factor=1e-5

    smoking_df %>% bind_rows(nonsmoking_df, .id="key") %>% filter(tested_combo == paste0(gt,"_",mg)) %>%
    mutate(key = ifelse(key==1, "Ever-smoker","Never-smoker")) %>%
    mutate(x_label = ifelse(epistatic_gt == "WT", "WT", "Mut"),
           x_label = ifelse(key == "Ever-smoker", paste0(x_label," "), x_label),
           x_label = factor(x_label, levels = c("WT ","Mut ","WT","Mut"))) %>% {
    # mutate(nudge_dist = scale((1:n())/n()**(1/.6), scale=FALSE), .by = c(key, epistatic_gt),
    #         x_label = ifelse(epistatic_gt == "WT", paste0(gt,"\nWT"), paste0(gt,"-\nmutant")),
    #         x_label = ifelse(key == "Ever-smoker", paste0(x_label," "), x_label),
    #         x_label = factor(x_label, levels = c(paste0(gt,"\nWT "), paste0(gt,"-\nmutant "),paste0(gt,"\nWT"), paste0(gt,"-\nmutant")))) %>% {
        ggplot(.,aes(x=x_label, y = gamma_mle)) + 
            geom_errorbar(aes(ymin=gamma_ci_low,ymax=gamma_ci_high),
                            width=0,linewidth=0.3) + 
            geom_point(aes(fill=key, size=to_count),
                            color="black",
                            shape=21) + 
            geom_vline(xintercept = 2.5, lty=2,col="grey") +
            scale_size_continuous(range=c(2,8)) +
            scale_fill_manual(values = get_smoker_nonsmoker_palette()) +
            expand_limits(y = 0) +
            scale_y_continuous(labels=function(y)y*y_scale_factor) +
            # scale_size_continuous() +
            #scale_size_binned(n.breaks = 2) +
            labs(y=TeX(r'(Scaled selection coefficient $(\times 10^{5})$)'), x=paste0('*',gt,'* genotype'),title= paste0('Selection for *',mg,'* mutations'), size="") +#"Sample count") +
            guides(size=guide_legend(nrow=1), fill="none") +#guide_legend(override.aes = list(size=6))) +
            theme_classic() +
            theme(
                plot.title = ggtext::element_markdown(size = 22, hjust=0.5),
                axis.text = element_text(size=18), 
                axis.title.y = element_text(size = 22),
                axis.title.x = ggtext::element_markdown(size = 22),
                legend.text = element_text(size=14),
                legend.position = "bottom")
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# G x G x G

find_extra_effects = function(mg, gt, df, spread=0.8){
    if(df %>% filter(mutated_gene == mg, str_detect(epistatic_gt, gt)) %>% nrow > 0){
        return(
        df %>% 
        filter(mutated_gene == mg, str_detect(epistatic_gt, gt)) %>%
        group_by(gene_set) %>%
        mutate(extra_effect_strict = max(gamma_ci_low) > min(gamma_ci_high)) %>%
        filter(str_detect(epistatic_gt,"_")) %>%
        
        mutate(pairwise_combo = paste0(mg,' [',gt,']')) %>%
        arrange(gene_set) %>%
        group_by(pairwise_combo) %>%
        mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE), 
                label_nudge_dist = ifelse(nudge_dist>0,nudge_dist+0.125,nudge_dist-0.125)) %>%
        mutate(genes_to_label = gsub("_","",gsub(gt,"",epistatic_gt))))
    }
}

# @param epistatic_pairs: pair of genes in the form of [prior mutation, subsequent mutation]
#                           additive and higher-order interactions will be evaluated against
#                           the pairwise interactions between the two genes of the pair
get_multi_gene_effects = function(epistatic_pairs, M3_df, M2_df){
    mutant_genes = str_split_i(epistatic_pairs,'_',2)
    epi_genes = str_split_i(epistatic_pairs,'_',1)

    plotting_df = rbindlist(lapply(1:length(epistatic_pairs), 
                    function(i){find_extra_effects(mutant_genes[i], epi_genes[i], M3_df)}))

    plotting_df = plotting_df %>%
        bind_rows(M2_df %>% filter(tested_combo %in% epistatic_pairs) %>% mutate(pairwise_combo = combo_name, nudge_dist = 0)) %>%
        mutate(group = ifelse(epistatic_gt=="WT","WT", ifelse(str_detect(epistatic_gt, "_"), "double","single"))) %>%
        arrange(pairwise_combo)

    return(plotting_df)
}


find_extra_effects_2 = function(mg, gt, df, spread=0.8){
    if(df %>% filter(mutated_gene == mg, str_detect(epistatic_gt, gt)) %>% nrow > 0){
        return(
        df %>% 
        filter(mutated_gene == mg, str_detect(epistatic_gt, gt)) %>%
        group_by(gene_set) %>%
        mutate(extra_effect_strict = max(gamma_ci_low) > min(gamma_ci_high)) %>%
        ungroup() %>%
        mutate(pairwise_combo = paste0(mg,' [',gt,']'),
               genes_to_label = gsub("_","",gsub(gt,"",epistatic_gt))))
    }
}

# @param epistatic_pairs: pair of genes in the form of [prior mutation, subsequent mutation]
#                           additive and higher-order interactions will be evaluated against
#                           the pairwise interactions between the two genes of the pair
get_multi_gene_effects_2 = function(epistatic_pairs, M3_df){
    mutant_genes = str_split_i(epistatic_pairs,'_',2)
    epi_genes = str_split_i(epistatic_pairs,'_',1)

    plotting_df = rbindlist(lapply(1:length(epistatic_pairs), 
                    function(i){find_extra_effects_2(mutant_genes[i], epi_genes[i], M3_df)}))

    #plotting_df = plotting_df %>% bind_rows(M3_df %>% filter(epistatic_gt=="WT") %>% mutate(pairwise_combo = combo_name))# %>% distinct(pick(gene_set, mutated_gene),.keep_all = TRUE))

    return(plotting_df)
}



# create_plotting_df = function(mg, gt, df, spread=0.8){
#     # tmp is a list of all two genes to third combinations where the mutated gene is the third
#     #   and the pairwise interactor is in the initial two
#     tmp = df %>% 
#         filter(str_count(tested_combo,"_")>1) %>% # two genes to third (110 -> 111)
#         filter(mutated_gene == mg) %>%
#         filter(str_detect(epistatic_gt, gt)) %>%
#         pull(tested_combo)
    
#     if(length(tmp)==0){return(data.frame())}

#     # returns df of only results from two genes to third
#     return(
#         df %>%
#         filter(tested_combo %in% tmp, epistatic_gt != "WT") %>%
#         mutate(pairwise_combo = paste0(mg,' [',gt,']')) %>%
#         arrange(gene_set) %>%
#         group_by(pairwise_combo) %>%
#         mutate(nudge_dist = scale((1:n())/n()**(1/spread), scale=FALSE), 
#                 label_nudge_dist = ifelse(nudge_dist>0,nudge_dist+0.125,nudge_dist-0.125)) %>%
#         mutate(genes_to_label = gsub("_","",gsub(gt,"",epistatic_gt)))
#     )
# }

# find_extra_effects = function(data){
#     pairwise_row = data %>% filter(group == "single")
#     extra_gene_rows = data %>% filter(group == "double")

#     extra_gene_rows = 
#         extra_gene_rows %>%
#         mutate(difference_from_pairwise = gamma_mle - pairwise_row$gamma_mle,
#                 ratio_to_pairwise = gamma_mle/pairwise_row$gamma_mle, 
#                 # confidence interval for double mutation does not overlap with MLE for pairwise interaction
#                 extra_effect = gamma_ci_high < pairwise_row$gamma_mle |
#                                 gamma_ci_low > pairwise_row$gamma_mle,
#                 extra_effect_strict = gamma_ci_high < pairwise_row$gamma_ci_low | gamma_ci_low > pairwise_row$gamma_ci_high,
#                 signif = any(signif,extra_effect_strict))
    
#     return(data %>% filter(group != "double") %>% bind_rows(extra_gene_rows))
# }

# # @param epistatic_pairs: pair of genes in the form of [prior mutation, subsequent mutation]
# #                           additive and higher-order interactions will be evaluated against
# #                           the pairwise interactions between the two genes of the pair
# get_multi_gene_effects = function(epistatic_pairs, M3_df, M2_df){
#     mutant_genes = str_split_i(epistatic_pairs,'_',2)
#     epi_genes = str_split_i(epistatic_pairs,'_',1)

#     plotting_df = rbindlist(lapply(1:length(epistatic_pairs), 
#                     function(i){create_plotting_df(mutant_genes[i], epi_genes[i], M3_df)}))

#     plotting_df = plotting_df %>%
#         bind_rows(M2_df %>% filter(tested_combo %in% epistatic_pairs) %>% mutate(pairwise_combo = combo_name, nudge_dist = 0)) %>%
#         mutate(group = ifelse(epistatic_gt=="WT","WT", ifelse(str_detect(epistatic_gt, "_"), "double","single"))) %>%
#         arrange(pairwise_combo)

#     gene_set_list = split(plotting_df, plotting_df$pairwise_combo)
#     # currently find extra effects is not working becuase is comparing M=3 results to M=2 results
#     plotting_df = rbindlist(lapply(gene_set_list, find_extra_effects))

#     return(plotting_df)
# }

# Plot additive and higher-order epistatic interactions
# The from-WT and pairwise points come from the M=2 model, which represents an average of all M=3 models
plot_multi_gene_effects = function(df, strict=TRUE, interactions_to_plot=NULL){
# TODO: add option to plot only interactions that have an extra effect

    alpha_palette = c("TRUE" = 1, "FALSE" = 0.1)
    fill_palette = c("WT (no epistasis)" = hue_pal()(3)[1], "Pairwise" = hue_pal()(3)[3], "Higher-order" = hue_pal()(3)[2])
    size_palette = c("WT (no epistasis)" = 5, "Pairwise" = 5, "Higher-order" = 3)

    if(!is.null(interactions_to_plot)){
        df = df %>% filter(pairwise_combo %in% interactions_to_plot) %>%
                    mutate(pairwise_combo = factor(pairwise_combo, interactions_to_plot))
    }

    df = df %>%
            mutate(group_label = ifelse(group == "WT", "WT (no epistasis)", ifelse(group == "single", "Pairwise", "Higher-order")),
                    #extra_effect = ifelse(is.na(extra_effect), FALSE, extra_effect),
                    extra_effect_strict = ifelse(is.na(extra_effect_strict), FALSE, extra_effect_strict)) %>%
            rowwise() %>% mutate(show_effect = extra_effect_strict) %>% ungroup()
    
    plot = 
        df %>%
        ggplot(aes(x = pairwise_combo, y=gamma_mle)) + 
            geom_errorbar(aes(ymin = gamma_ci_low, ymax = gamma_ci_high, alpha= (group != "double" | show_effect)),
                            width=0, linewidth = 0.2, position = position_nudge(df$nudge_dist)) + 
            geom_point(aes(fill=group_label, size=group_label, alpha= (group != "double" | show_effect)), 
                            shape = 21, position = position_nudge(df$nudge_dist)) +
            geom_text(data = df %>% filter(group=="double", show_effect),
                                        aes(label=genes_to_label), 
                                        position = position_nudge(df %>% filter(group=="double", show_effect) %>% pull(label_nudge_dist)),
                                        check_overlap = TRUE) +
            labs(x="Mutated gene [somatic genotype]", y="Scaled selection coefficient") +
            scale_alpha_manual(values = alpha_palette, name="Significant Difference in Selection") +
            scale_fill_manual(values = fill_palette, name="Genetic Context") +
            scale_size_manual(values = size_palette, name="Genetic Context") +
            guides(alpha="none") + 
            scale_y_continuous(labels=scientific_expr) +
            theme_classic() +
            theme(plot.title = element_text(size = 24, hjust=0.5),
                        axis.title = element_text(size = 20),
                        axis.text = element_text(size = 16),
                        legend.position = "right",#c(0.6,0.2),
                        # legend.direction="horizontal",
                        legend.key.size = unit(1.5, 'cm'),
                        legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18),
                        panel.grid.major.y = element_line(color="gray",linewidth=0.75, linetype=3)) +
            coord_flip()
    
    return(plot)
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

#' Helper Functions

# With modification from https://groups.google.com/g/ggplot2/c/a_xhMoQyxZ4
fancy_scientific = function(l) {
    l = format(l, scientific = TRUE)
    l = gsub("^0(.*)e(.*)", "0", l) # Just show 0 when is 0e...
    l = gsub("^(.*)e", "'\\1'e", l)
    l = gsub("\'","",l)
    l = gsub("\\+","",l) # modification to remove unnecessary pluses
    l = gsub("e", "%*%10^", l)
    return(l)
}

scientific_expr = function(l){parse(text=fancy_scientific(l))}

log_labels = function(l) {
    l = format(l, scientific = TRUE)
    l = gsub("\\+","",l)
    l = gsub("^1e", "10^", l)
    return(parse(text=l))
}

add_default_theme = function(p){
    p = p + 
        theme(axis.title = element_text(size = 18), 
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 20, hjust = 0.5))
}

sym_diff = function(s1,s2){unique(c(setdiff(s1,s2), setdiff(s2,s1)))}

# Welch's t-test
t_test = function(mu1, mu2, sd1, sd2, n1, n2, mu0 = 0, detailed=F){
    t_stat = (mu1 - mu2 - mu0) / sqrt(((sd1^2)/n1) + ((sd2^2)/n2))
    df = (((sd1^2)/n1 + (sd2^2)/n2)^2) / (((sd1^4)/(n1^2*(n1 - 1))) + ((sd2^4)/(n2^2*(n2 - 1))))
    p_val = 2 * pt(-abs(t_stat), df)

    if(!detailed){
        return(p_val)
    } else{
        return(list(t_stat, p_val))
    }
}

head_and_tail = function(df, n=6, n_head=n, n_tail=n){rbind(head(df, n_head),tail(df, n_tail))}

# From https://stackoverflow.com/a/34096575 
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


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
