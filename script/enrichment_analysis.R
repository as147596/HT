enrichment_analysis<-function(list, background, pathway_data) {
  library(dplyr)
  pathway_data <- pathway_data %>% filter(metabolite %in% background)
  
  total_metabolite <- length(unique(pathway_data$metabolite))
  
  pathway_enrichment <- pathway_data %>%
    group_by(pathway) %>%
    summarize(
      Pathway_Size = n(),
      Overlap_Size = sum(metabolite %in% list),
      P_Value = phyper(
        q = sum(metabolite %in% list) - 1,
        m = length(metabolite),
        n = total_metabolite - length(metabolite),
        k = length(list),
        lower.tail = FALSE
      )
    ) %>%
    arrange(P_Value)
  
  return(pathway_enrichment)
}

