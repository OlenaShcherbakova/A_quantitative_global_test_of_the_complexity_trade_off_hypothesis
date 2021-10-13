# Code accompanying the paper *A quantitative global test of the complexity trade-off hypothesis* by Olena Shcherbakova, Volker Gast, Damián E. Blasi‬, Hedvig Skirgård, Russell D. Gray, and Simon J. Greenhill

**Overview of the files:**

Input files can be found in the “data” folder:
 - GB.tsv: the list of 302 Grambank languages with computed metric scores and with values for each Grambank feature
 - languages_glottolog.csv stems from Glottolog 4.4 and offers information on languages, such as language family (Family_ID), Latitude, and Longitude 
 - folder “phylogenies” contains 5 folders with respective phylogenies with taxa and posterior.trees files in each. These files were obtained from DPlace (Kirby et al. 2016) repository: https://github.com/D-PLACE/dplace-data.

Files within “results_bayestraits” folder contain the input and output files for analyses carried out within BayesTraits .

**Overview of the scripts:**

Before replicating the study or parts of it, run library.R first. This script loads R packages and custom functions, in particular for loading trees, reading in files, and calculating metric scores. 
The rest of the scripts can be run in no particular order.

 - setup.R provides code for creating the input files that were used in the analyses carried out with the help of BayesTraits

 - world_tree_transformation.R: transforms the global tree (Jäger 2018) from newick (world.tre) into nexus format (posterior.trees), so that the name of the global tree is the same as those of individual trees, making navigating between phylogenies in the scripts easier. However, both tree files are already available in “data/phylogenies/world”.

 - preparing_world_subset.R is used 1) to subset the initial input file (GB.tsv; 302 languages) to the sample of languages (244 languages) that are matched with the global tree (this script is used within some of the figure-generating scripts) and 2) to create a list (Languages_Appendix.csv) of 244 languages with corresponding metric scores that used in the analyses and plots based on the global tree

Figure- and table-generating scripts: 
 - script_scatter_plot.R (Figure 1)
 - script_histogram_plot.R (Figure 2)
 - script_map_plot.R (Figure 3)
 - script_tree_plot.R (Figure 4)
 - script_histograms_set_plot.R (Figure A in Supplementary Materials)
 - script_values_table.R runs the number of scripts that read and summarize the files from BayesTraits analyses output (alphas_distribution.R, sigmas_distribution.R, correlation_distribution.R, Marginal_Likelihoods_BF.R) to generate Table D and a heatmap with interdomain correlations on the global tree (Figure 5)
 - measuring_phylo_signal.R generates Table C with two measures of phylogenetic signal for nominal and verbal scores across phylogenies

References

Kirby, Kathryn R., Russell D. Gray, Simon J. Greenhill, Fiona M. Jordan, Stephanie Gomes-Ng, Hans-Jörg Bibiko, Damián E. Blasi, Carlos A. Botero, Claire Bowern, Carol R. Ember, Dan Leehr, Bobbi S. Low, Joe McCarter, William Divale & Michael C. Gavinet. D-PLACE: A global database of cultural, linguistic and environmental diversity. *PloS one* 11.7 (2016): e0158391.
