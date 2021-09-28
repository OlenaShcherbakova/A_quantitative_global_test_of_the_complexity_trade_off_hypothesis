#preparing world subset to be used for analyses and plots based on global tree + a spreadsheet of relevant languages

source("library.R")

gb <- load_data("data/GB.tsv")
taxa <- read.csv("data/phylogenies/world/taxa.csv") #loading ASJP file (v. 17)
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("languages_glottolog.csv", header=TRUE)[c("ID", "Family_ID", "Latitude", "Longitude", "Name")]

# merge datasets
gb_world_subset <- merge(gb, language_glottolog, by.x="Glottocode", by.y="ID")
gb_world_subset <- gb_world_subset %>%
dplyr::filter(!is.na(taxon))

#load tree
TREEFILE <- "data/phylogenies/world/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]
#merge data with tree and set up tree
gb_world_subset <- gb_world_subset[gb_world_subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_world_subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_world_subset <- gb_world_subset %>%
mutate(across(2:18, round, 2))

gb_world_subset_spreadsheet <- gb_world_subset %>%
  rename(Nominal_words_score= Nominal_words_complexity, Verbal_score = Verbal_complexity) %>%
  select(Name, Glottocode, Family_ID, Nominal_words_score, Verbal_score)

write.csv(gb_world_subset_spreadsheet, file = here("output","Languages_Appendix.csv"), row.names = FALSE)
