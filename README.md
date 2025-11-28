# Gut Microbiome and Metabolite Flux Analysis
## Introduction
This repository contains a complete set of scripts for gut microbiome and metabolic flux analysis.

## Main scripts include:

**microbiome_diversity_taxa_diff.R:**
Performs gut microbiome alpha and beta diversity analyses, along with differential taxa identification.

**flux_preFiltering.R:**
Filters samples prior to metabolic flux prediction to ensure data quality.

**key_metabolite_identification.py:**
Applies DoubleML (Double Machine Learning) to infer causal relationships between metabolic fluxes and blood pressure.

**key_metabolite_identification.R:**
Conducts key metabolic flux analysis and generates related visualizations.

**microbiome_differential_network.R:**
Constructs differential microbial networks and identifies key metabolite-associated guilds, followed by visualization.

**key_metabolite_associated_species_validation.py:**
Utilizes DoubleML to perform causal inference of key guilds’ metabolic fluxes on blood pressure within the cross-feeding network.

**key_metabolite_associated_species_validation.R:**
Visualizes and validates differential metabolic fluxes of key guilds in the cross-feeding network.

**cross_feeding_network_analysis.R:**
Analyzes causal relationships between cross-feeding network edge weights and blood pressure, and performs pathway enrichment analysis of differential metabolites.
