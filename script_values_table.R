source('library.R')
source('alphas_distribution.R')
source('sigmas_distribution.R')
source('correlation_distribution.R')
source('Marginal_Likelihoods_BF.R')

parameters <- cbind(alphas, sigmas, correlations, likelihoods)
distributions <- parameters
distributions$BF <- 2 * (distributions$Lh_dep - distributions$Lh_indep)
distributions <- round(distributions, digits = 2)

row_names <- as.vector(results.indep.world)
row.names(distributions) <- row_names

#providing a column corresponding to the file name (First trait used in the analysis, second trait used in the analysis, and Phylogeny)
distributions$Traits <- c('X: Nominal words Y: Verbs (Indo-European 1)', 'X: Nominal words Y: Verbs (Indo-European 2)', 'X: Nominal words Y: Verbs (Austronesian)', 
  'X: articles Y: arguments', 'X: articles Y: arguments_non_core', 'X: articles Y: aspect', 'X: articles Y: clause_v', 'X: articles Y: mood', 'X: articles Y: negation', 'X: articles Y: tense', 'X: articles Y: transitivity',
                              
                              'X: case Y: arguments', 'X: case Y: arguments_non_core', 'X: case Y: aspect', 'X: case Y: clause_v', 'X: case Y: mood', 'X: case Y: negation', 'X: case Y: tense', 'X: case Y: transitivity',
                              
                              'X: gender Y: arguments', 'X: gender Y: arguments_non_core', 'X: gender Y: aspect', 'X: gender Y: clause_v', 'X: gender Y: mood', 'X: gender Y: negation', 'X: gender Y: tense', 'X: gender Y: transitivity',
                              
                              'X: Nominal words Y: Verbs (Global tree)',
                              
                              'X: number Y: arguments', 'X: number Y: arguments_non_core', 'X: number Y: aspect', 'X: number Y: clause_v', 'X: number Y: mood', 'X: number Y: negation', 'X: number Y: tense', 'X: number Y: transitivity',
                              
                              'X: possession Y: arguments', 'X: possession Y: arguments_non_core', 'X: possession Y: aspect', 'X: possession Y: clause_v', 'X: possession Y: mood', 'X: possession Y: negation', 'X: possession Y: tense', 'X: possession Y: transitivity', 'X: Nominal words Y: Verbs (Sino-Tibetan)')

distributions <- distributions %>% dplyr::relocate(Traits)

opening_bracket <- "("
closing_bracket <- ")"

#putting lower and upper intervals into one cell
distributions$correlation_intervals <- paste (distributions$lower_correlation_dep, distributions$upper_correlation_dep, sep = " ", collapse = NULL) 
distributions$correlation_int <- paste (opening_bracket, distributions$correlation_interval, closing_bracket, sep = "", collapse = NULL) 
distributions$correlation <- paste (distributions$median_correlation_dep, distributions$correlation_int, sep = " ", collapse = NULL) 

distributions$alpha1_dep_intervals <- paste (distributions$lower_alpha1_dep, distributions$upper_alpha1_dep, sep = " ", collapse = NULL) 
distributions$alpha1_dep_int <- paste (opening_bracket, distributions$alpha1_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha1_dep <- paste (distributions$median_alpha1_dep, distributions$alpha1_dep_int, sep = " ", collapse = NULL) 

distributions$alpha1_indep_intervals <- paste (distributions$lower_alpha1_indep, distributions$upper_alpha1_indep, sep = " ", collapse = NULL) 
distributions$alpha1_indep_int <- paste (opening_bracket, distributions$alpha1_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha1_indep <- paste (distributions$median_alpha1_indep, distributions$alpha1_indep_int, sep = " ", collapse = NULL) 

distributions$alpha2_dep_intervals <- paste (distributions$lower_alpha2_dep, distributions$upper_alpha2_dep, sep = " ", collapse = NULL) 
distributions$alpha2_dep_int <- paste (opening_bracket, distributions$alpha2_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha2_dep <- paste (distributions$median_alpha2_dep, distributions$alpha2_dep_int, sep = " ", collapse = NULL) 

distributions$alpha2_indep_intervals <- paste (distributions$lower_alpha2_indep, distributions$upper_alpha2_indep, sep = " ", collapse = NULL) 
distributions$alpha2_indep_int <- paste (opening_bracket, distributions$alpha2_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$alpha2_indep <- paste (distributions$median_alpha2_indep, distributions$alpha2_indep_int, sep = " ", collapse = NULL) 

distributions$sigma1_dep_intervals <- paste (distributions$lower_sigma1_dep, distributions$upper_sigma1_dep, sep = " ", collapse = NULL) 
distributions$sigma1_dep_int <- paste (opening_bracket, distributions$sigma1_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma1_dep <- paste (distributions$median_sigma1_dep, distributions$sigma1_dep_int, sep = " ", collapse = NULL) 

distributions$sigma1_indep_intervals <- paste (distributions$lower_sigma1_indep, distributions$upper_sigma1_indep, sep = " ", collapse = NULL) 
distributions$sigma1_indep_int <- paste (opening_bracket, distributions$sigma1_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma1_indep <- paste (distributions$median_sigma1_indep, distributions$sigma1_indep_int, sep = " ", collapse = NULL) 

distributions$sigma2_dep_intervals <- paste (distributions$lower_sigma2_dep, distributions$upper_sigma2_dep, sep = " ", collapse = NULL) 
distributions$sigma2_dep_int <- paste (opening_bracket, distributions$sigma2_dep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma2_dep <- paste (distributions$median_sigma2_dep, distributions$sigma2_dep_int, sep = " ", collapse = NULL) 

distributions$sigma2_indep_intervals <- paste (distributions$lower_sigma2_indep, distributions$upper_sigma2_indep, sep = " ", collapse = NULL) 
distributions$sigma2_indep_int <- paste (opening_bracket, distributions$sigma2_indep_intervals, closing_bracket, sep = "", collapse = NULL) 
distributions$sigma2_indep <- paste (distributions$median_sigma2_indep, distributions$sigma2_indep_int, sep = " ", collapse = NULL)


distributions_subset <- distributions[c(1, 29:31, 34, 37, 40, 43, 46, 49, 52, 55, 58)]

distributions_subset_long <- distributions_subset %>% 
  gather(Parameter, Value, c(2:3, 6:13)) %>% 
  separate(col = Parameter, into = c("Par", "Model"), sep = "_") %>%
  spread(key = "Par",
         value = "Value") %>%
  rename(Marginal_Likelihood = Lh) %>%
  dplyr::select(Traits, Model, BF, alpha1, sigma1, alpha2, sigma2, correlation)

distributions_subset_long$correlation <- ifelse(distributions_subset_long$Model == "dep", distributions_subset_long$correlation, " ")

#reorganizing the rows: placing the values for Nominal words and Verbs at the top of the dataframe with other traits combinations following them
NwV_subset <- distributions_subset_long[grep("X: Nominal words Y: Verb", distributions_subset_long$Traits), ]

distributions_subset_long_without <- distributions_subset_long[-grep("X: Nominal words Y: Verb", distributions_subset_long$Traits), ]

distributions_subset_final <- rbind(NwV_subset, distributions_subset_long_without)

write.csv(distributions_subset_final, file = here("output","Table_D.csv"), row.names = FALSE)


#Creating a heatmap out of the correlation coefficients for interdomain pairs, while taking into account Bayes Factor values: when lower than 2, correlation is assigned to be zero as the relationship is not supported; otherwise the correlation coefficient is taken from another column.)
heatmap_subset <- distributions[-grep("X: Nominal words Y: Verbs", distributions$Traits), ]
heatmap_subset$cor_heatmap <- ifelse(heatmap_subset$BF < 2, 0, heatmap_subset$median_correlation_dep)
vector <- heatmap_subset$cor_heatmap
vector <- vector[1:40]
m <- matrix(vector, nrow=5, byrow=TRUE, dimnames = list(c("articles","case","gender", "number", "possession"), c("core arg", "non-core arg", "aspect", "clause", "mood", "negation", "tense", "transitivity"))) 
heatmap <- pheatmap(m, fontsize = 16, treeheight_row = 0, treeheight_col = 0)
ggsave(file = here("output","heatmap.png"), plot = heatmap, dpi=300, dev="png")
