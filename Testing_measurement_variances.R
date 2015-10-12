## This is is one of several files containing scripts and functions used in processing and analysis of data for Matthew Dufort's Ph.D. dissertation at the University of Minnesota, titled "Coexistence, Ecomorphology, and Diversification in the Avian Family Picidae (Woodpeckers and Allies)."

# this file contains scripts to check the within-species or within-sex-within-species? and among-species variances for each morphological variables
# the purpose of this is to detect measurements that may have greater measurement error or greater variation within populations and species

## load morphological data generated in "Morphology_data_processing.R" scripts
load(file='Picidae_data_for_distribution_morphology_evolution.RData')


library(lme4)  # load lme4 package


## fit ANOVA models to the each measurement, with sex nested within species as the explanatory variables

aov.Maxilla_depth_ramp_absent <- lmer(Maxilla_depth_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Maxilla_depth_ramp_absent
summary(aov.Maxilla_depth_ramp_absent)  # output summary of ANOVA model
(summary(aov.Maxilla_depth_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Maxilla_depth_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Maxilla_length_ramp_absent <- lmer(Maxilla_length_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Maxilla_length_ramp_absent
summary(aov.Maxilla_length_ramp_absent)  # output summary of ANOVA model
(summary(aov.Maxilla_length_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Maxilla_length_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Cranium_GL_ramp_absent <- lmer(Cranium_GL_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Cranium_GL_ramp_absent
summary(aov.Cranium_GL_ramp_absent)  # output summary of ANOVA model
(summary(aov.Cranium_GL_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Cranium_GL_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Cranium_OC_ramp_absent <- lmer(Cranium_OC_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Cranium_OC_ramp_absent
summary(aov.Cranium_OC_ramp_absent)  # output summary of ANOVA model
(summary(aov.Cranium_OC_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Cranium_OC_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Mandible_mean_ramp_absent <- lmer(Mandible_mean_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Mandible_mean_ramp_absent
summary(aov.Mandible_mean_ramp_absent)  # output summary of ANOVA model
(summary(aov.Mandible_mean_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Mandible_mean_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Humerus_mean <- lmer(Humerus_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Humerus_mean
summary(aov.Humerus_mean)  # output summary of ANOVA model
(summary(aov.Humerus_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Humerus_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Ulna_mean <- lmer(Ulna_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Ulna_mean
summary(aov.Ulna_mean)  # output summary of ANOVA model
(summary(aov.Ulna_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Ulna_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Femur_mean <- lmer(Femur_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Femur_mean
summary(aov.Femur_mean)  # output summary of ANOVA model
(summary(aov.Femur_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Femur_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Tibiotarsus_mean <- lmer(Tibiotarsus_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Tibiotarsus_mean
summary(aov.Tibiotarsus_mean)  # output summary of ANOVA model
(summary(aov.Tibiotarsus_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Tibiotarsus_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Tarsometatarsus_mean <- lmer(Tarsometatarsus_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Tarsometatarsus_mean
summary(aov.Tarsometatarsus_mean)  # output summary of ANOVA model
(summary(aov.Tarsometatarsus_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Tarsometatarsus_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Pygostyle_width <- lmer(Pygostyle_width ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Pygostyle_width
summary(aov.Pygostyle_width)  # output summary of ANOVA model
(summary(aov.Pygostyle_width)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Pygostyle_width"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.Pygostyle_length <- lmer(Pygostyle_length ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$all_inds, na.action=na.omit)  # fit ANOVA to Pygostyle_length
summary(aov.Pygostyle_length)  # output summary of ANOVA model
(summary(aov.Pygostyle_length)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$all_inds[["Pygostyle_length"]], na.rm=TRUE)  # calculate percentage of variance not explained by model


## rerun these ANOVA models only for length and depth of maxilla, using only cases where the sex is "M" or "F"

picidae.morph_combined.log.reduced_var_inds.mf_only <- subset(picidae.morph_combined.log.reduced_var_inds$all_inds, subset=(Sex %in% c("M","F")))  # subset the morphological data, reatining only cases where the Sex is "M" or "F"

aov.mf_only.Maxilla_depth_ramp_absent <- lmer(Maxilla_depth_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.mf_only, na.action=na.omit)  # fit ANOVA to Maxilla_depth_ramp_absent, using only cases where the sex is "M" or "F"
summary(aov.mf_only.Maxilla_depth_ramp_absent)  # output summary of ANOVA model
(summary(aov.mf_only.Maxilla_depth_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.mf_only[["Maxilla_depth_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.mf_only.Maxilla_length_ramp_absent <- lmer(Maxilla_length_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.mf_only, na.action=na.omit)  # fit ANOVA to Maxilla_length_ramp_absent, using only cases where the sex is "M" or "F"
summary(aov.mf_only.Maxilla_length_ramp_absent)  # output summary of ANOVA model
(summary(aov.mf_only.Maxilla_length_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.mf_only[["Maxilla_length_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model


## rerun these ANOVA models only for length and depth of maxilla, using data with ramp_absent imputed by stochastic imputation (see Morphology_data_processing.R)

aov.imputed.Maxilla_depth_ramp_absent <- lmer(Maxilla_depth_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, na.action=na.omit)  # fit ANOVA to Maxilla_depth_ramp_absent, using data with ramp_absent imputed by stochastic imputation
summary(aov.imputed_resids_mf.Maxilla_depth_ramp_absent)  # output summary of ANOVA model
(summary(aov.imputed.Maxilla_depth_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.imputed$all_inds[["Maxilla_depth_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.imputed.Maxilla_length_ramp_absent <- lmer(Maxilla_length_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.imputed$all_inds, na.action=na.omit)  # fit ANOVA to Maxilla_length_ramp_absent, using data with ramp_absent imputed by stochastic imputation
summary(aov.imputed_resids_mf.Maxilla_length_ramp_absent)  # output summary of ANOVA model
(summary(aov.imputed.Maxilla_length_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.imputed$all_inds[["Maxilla_length_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model


## rerun these ANOVA models only for length and depth of maxilla and length and width of pygostyle, using complete cases only

aov.complete_ind_only.Maxilla_depth_ramp_absent <- lmer(Maxilla_depth_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$complete_ind_only, na.action=na.omit)  # fit ANOVA to Maxilla_depth_ramp_absent, using complete cases only
summary(aov.complete_ind_only.Maxilla_depth_ramp_absent)  # output summary of ANOVA model
(summary(aov.complete_ind_only.Maxilla_depth_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$complete_ind_only[["Maxilla_depth_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.complete_ind_only.Maxilla_length_ramp_absent <- lmer(Maxilla_length_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$complete_ind_only, na.action=na.omit)  # fit ANOVA to Maxilla_length_ramp_absent, using complete cases only
summary(aov.complete_ind_only.Maxilla_length_ramp_absent)  # output summary of ANOVA model
(summary(aov.complete_ind_only.Maxilla_length_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$complete_ind_only[["Maxilla_length_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.complete_ind_only.Pygostyle_length <- lmer(Pygostyle_length ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$complete_ind_only, na.action=na.omit)  # fit ANOVA to Pygostyle_length, using complete cases only
summary(aov.complete_ind_only.Pygostyle_length)  # output summary of ANOVA model
(summary(aov.complete_ind_only.Pygostyle_length)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$complete_ind_only[["Pygostyle_length"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.complete_ind_only.Pygostyle_width <- lmer(Pygostyle_width ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds$complete_ind_only, na.action=na.omit)  # fit ANOVA to Pygostyle_width, using complete cases only
summary(aov.complete_ind_only.Pygostyle_width)  # output summary of ANOVA model
(summary(aov.complete_ind_only.Pygostyle_width)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds$complete_ind_only[["Pygostyle_width"]], na.rm=TRUE)  # calculate percentage of variance not explained by model



## rerun ANOVA analyses only for length and depth of the maxilla, with complete cases only, with ramp_absent imputed using stochastic imputation

aov.complete_ind_only.imputed.Maxilla_depth_ramp_absent <- lmer(Maxilla_depth_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.imputed$complete_ind_only, na.action=na.omit)  # fit ANOVA to Maxilla_depth_ramp_absent, using complete cases only
summary(aov.complete_ind_only.imputed.Maxilla_depth_ramp_absent)  # output summary of ANOVA model
(summary(aov.complete_ind_only.imputed.Maxilla_depth_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.imputed$complete_ind_only[["Maxilla_depth_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.complete_ind_only.imputed.Maxilla_length_ramp_absent <- lmer(Maxilla_length_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.imputed$complete_ind_only, na.action=na.omit)  # fit ANOVA to Maxilla_length_ramp_absent, using complete cases only
summary(aov.complete_ind_only.imputed.Maxilla_length_ramp_absent)  # output summary of ANOVA model
(summary(aov.complete_ind_only.imputed.Maxilla_length_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.imputed$complete_ind_only[["Maxilla_length_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model


### rerun these ANOVA models, excluding the species in my community-level analyses, because those species have larger sample sizes spread across a number of populations

## generate a data frame without the species used in community-level analyses
comm_species <- c("Colaptes_auratus", "Colaptes_cafer", "Melanerpes_carolinus", "Melanerpes_erythrocephalus", "Melanerpes_formicivorus", "Dryobates_pubescens", "Leuconotopicus_villosus", "Leuconotopicus_borealis", "Dryocopus_pileatus", "Sphyrapicus_ruber", "Sphyrapicus_thyroideus")  # create a vector of taxon names included in community analyses
picidae.morph_combined.log.reduced_var_inds.no_comm_species <- picidae.morph_combined.log.reduced_var_inds$all_inds[!(picidae.morph_combined.log.reduced_var_inds$My_genus_species %in% comm_species),]  # create a new data frame with those community species excluded

aov.no_comm_species.Maxilla_depth_ramp_absent <- lmer(Maxilla_depth_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Maxilla_depth_ramp_absent, using complete cases only
summary(aov.no_comm_species.Maxilla_depth_ramp_absent)  # output summary of ANOVA model
(summary(aov.no_comm_species.Maxilla_depth_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Maxilla_depth_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Maxilla_length_ramp_absent <- lmer(Maxilla_length_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Maxilla_length_ramp_absent, using complete cases only
summary(aov.no_comm_species.Maxilla_length_ramp_absent)  # output summary of ANOVA model
(summary(aov.no_comm_species.Maxilla_length_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Maxilla_length_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Cranium_GL_ramp_absent <- lmer(Cranium_GL_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Cranium_GL_ramp_absent, using complete cases only
summary(aov.no_comm_species.Cranium_GL_ramp_absent)  # output summary of ANOVA model
(summary(aov.no_comm_species.Cranium_GL_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Cranium_GL_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Cranium_OC_ramp_absent <- lmer(Cranium_OC_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Cranium_OC_ramp_absent, using complete cases only
summary(aov.no_comm_species.Cranium_OC_ramp_absent)  # output summary of ANOVA model
(summary(aov.no_comm_species.Cranium_OC_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Cranium_OC_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Mandible_mean_ramp_absent <- lmer(Mandible_mean_ramp_absent ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Mandible_mean_ramp_absent, using complete cases only
summary(aov.no_comm_species.Mandible_mean_ramp_absent)  # output summary of ANOVA model
(summary(aov.no_comm_species.Mandible_mean_ramp_absent)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Mandible_mean_ramp_absent"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Humerus_mean <- lmer(Humerus_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Humerus_mean, using complete cases only
summary(aov.no_comm_species.Humerus_mean)  # output summary of ANOVA model
(summary(aov.no_comm_species.Humerus_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Humerus_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Ulna_mean <- lmer(Ulna_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Ulna_mean, using complete cases only
summary(aov.no_comm_species.Ulna_mean)  # output summary of ANOVA model
(summary(aov.no_comm_species.Ulna_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Ulna_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Femur_mean <- lmer(Femur_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Femur_mean, using complete cases only
summary(aov.no_comm_species.Femur_mean)  # output summary of ANOVA model
(summary(aov.no_comm_species.Femur_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Femur_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Tibiotarsus_mean <- lmer(Tibiotarsus_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Tibiotarsus_mean, using complete cases only
summary(aov.no_comm_species.Tibiotarsus_mean)  # output summary of ANOVA model
(summary(aov.no_comm_species.Tibiotarsus_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Tibiotarsus_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Tarsometatarsus_mean <- lmer(Tarsometatarsus_mean ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Tarsometatarsus_mean, using complete cases only
summary(aov.no_comm_species.Tarsometatarsus_mean)  # output summary of ANOVA model
(summary(aov.no_comm_species.Tarsometatarsus_mean)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Tarsometatarsus_mean"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Pygostyle_width <- lmer(Pygostyle_width ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Pygostyle_width, using complete cases only
summary(aov.no_comm_species.Pygostyle_width)  # output summary of ANOVA model
(summary(aov.no_comm_species.Pygostyle_width)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Pygostyle_width"]], na.rm=TRUE)  # calculate percentage of variance not explained by model

aov.no_comm_species.Pygostyle_length <- lmer(Pygostyle_length ~ 1 + (1 | My_genus_species/Sex), data = picidae.morph_combined.log.reduced_var_inds.no_comm_species, na.action=na.omit)  # fit ANOVA to Pygostyle_length, using complete cases only
summary(aov.no_comm_species.Pygostyle_length)  # output summary of ANOVA model
(summary(aov.no_comm_species.Pygostyle_length)$sigma^2) / var(picidae.morph_combined.log.reduced_var_inds.no_comm_species[["Pygostyle_length"]], na.rm=TRUE)  # calculate percentage of variance not explained by model
