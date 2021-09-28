# Setup for Bayestraits analyses

#script written by Simon Greenhill and modified by Olena Shcherbakova

#PHYLOGENY: Bouckaert et al 2012

# load data
library(here)
source(here("library.R"))
gb <- load_data(here('data/GB.tsv'))

# what are our analyses. This is a list of pairs that we want to analyse together
ANALYSES <- list(
  "Nom_words_Verbs"= c("Nominal_words_complexity_10", "Verbal_complexity_10")
)

COMMAND_DEPENDENT <- "4
2
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.dependent.log
Run
"

COMMAND_INDEPENDENT <- "4
2
TestCorrel
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.independent.log
Run
"

# the template for the cluster scheduler / job runner
SLURM_TEMPLATE <- "#!/bin/sh
#SBATCH --job-name bayestraits
#SBATCH --get-user-env
#SBATCH --output run.%%j.log
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1

BayesTraitsV3 %s %s < %s
"

# and the trees we care about...
PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

basedir <- basename(PHYLOGENIES[1])
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[1])


taxa <- read.csv(TAXAS[1])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Verbal_complexity_10", "Nominal_words_complexity_10", "case", "number", "gender", "possession", "articles", "arguments", "transitivity", "markers_arguments_non_core", "clause_v", "negation", "tense", "aspect", "mood"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))




#PHYLOGENY: Chang et al 2015
gb <- load_data(here('data/GB.tsv'))

# what are our analyses. This is a list of pairs that we want to analyse together
ANALYSES <- list(
  "Nom_words_Verbs"= c("Nominal_words_complexity_10", "Verbal_complexity_10")
)

COMMAND_DEPENDENT <- "4
2
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.dependent.log
Run
"

COMMAND_INDEPENDENT <- "4
2
TestCorrel
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.independent.log
Run
"

# the template for the cluster scheduler / job runner
SLURM_TEMPLATE <- "#!/bin/sh
#SBATCH --job-name bayestraits
#SBATCH --get-user-env
#SBATCH --output run.%%j.log
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1

BayesTraitsV3 %s %s < %s
"

# and the trees we care about...
PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

basedir <- basename(PHYLOGENIES[2])
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[2])


taxa <- read.csv(TAXAS[2])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Verbal_complexity_10", "Nominal_words_complexity_10", "case", "number", "gender", "possession", "articles", "arguments", "transitivity", "markers_arguments_non_core", "clause_v", "negation", "tense", "aspect", "mood"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))



#PHYLOGENY: Gray et al 2009
gb <- load_data(here('data/GB.tsv'))

# what are our analyses. This is a list of pairs that we want to analyse together
ANALYSES <- list(
  "Nom_words_Verbs"= c("Nominal_words_complexity_10", "Verbal_complexity_10")
)

COMMAND_DEPENDENT <- "4
2
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.dependent.log
Run
"

COMMAND_INDEPENDENT <- "4
2
TestCorrel
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.independent.log
Run
"

# the template for the cluster scheduler / job runner
SLURM_TEMPLATE <- "#!/bin/sh
#SBATCH --job-name bayestraits
#SBATCH --get-user-env
#SBATCH --output run.%%j.log
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1

BayesTraitsV3 %s %s < %s
"

# and the trees we care about...
PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

basedir <- basename(PHYLOGENIES[3])
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[3])


taxa <- read.csv(TAXAS[3])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Verbal_complexity_10", "Nominal_words_complexity_10", "case", "number", "gender", "possession", "articles", "arguments", "transitivity", "markers_arguments_non_core", "clause_v", "negation", "tense", "aspect", "mood"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))



#PHYLOGENY: Jäger 2018 (world/global tree)
gb <- load_data(here('data/GB.tsv'))

# what are our analyses. This is a list of pairs that we want to analyse together
ANALYSES <- list(
  "Nom_words_Verbs"= c("Nominal_words_complexity_10", "Verbal_complexity_10"),
  "case_arguments" = c("case", "arguments"),
  "case_transitivity" = c("case", "transitivity"),
  "case_arguments_non_core" = c("case", "markers_arguments_non_core"),
  "case_clause_v" = c("case", "clause_v"),
  "case_negation" = c("case", "negation"),
  "case_tense" = c("case", "tense"),
  "case_aspect" = c("case", "aspect"),
  "case_mood" = c("case", "mood"),
  
  "number_arguments" = c("number", "arguments"),
  "number_transitivity" = c("number", "transitivity"),
  "number_arguments_non_core" = c("number", "markers_arguments_non_core"),
  "number_clause_v" = c("number", "clause_v"),
  "number_negation" = c("number", "negation"),
  "number_tense" = c("number", "tense"),
  "number_aspect" = c("number", "aspect"),
  "number_mood" = c("number", "mood"),
  
  "gender_arguments" = c("gender", "arguments"),
  "gender_transitivity" = c("gender", "transitivity"),
  "gender_arguments_non_core" = c("gender", "markers_arguments_non_core"),
  "gender_clause_v" = c("gender", "clause_v"),
  "gender_negation" = c("gender", "negation"),
  "gender_tense" = c("gender", "tense"),
  "gender_aspect" = c("gender", "aspect"),
  "gender_mood" = c("gender", "mood"),
  
  "possession_arguments" = c("possession", "arguments"),
  "possession_transitivity" = c("possession", "transitivity"),
  "possession_arguments_non_core" = c("possession", "markers_arguments_non_core"),
  "possession_clause_v" = c("possession", "clause_v"),
  "possession_negation" = c("possession", "negation"),
  "possession_tense" = c("possession", "tense"),
  "possession_aspect" = c("possession", "aspect"),
  "possession_mood" = c("possession", "mood"),
  
  "articles_arguments" = c("articles", "arguments"),
  "articles_transitivity" = c("articles", "transitivity"),
  "articles_arguments_non_core" = c("articles", "markers_arguments_non_core"),
  "articles_clause_v" = c("articles", "clause_v"),
  "articles_negation" = c("articles", "negation"),
  "articles_tense" = c("articles", "tense"),
  "articles_aspect" = c("articles", "aspect"),
  "articles_mood" = c("articles", "mood")
)

COMMAND_DEPENDENT <- "4
2
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.dependent.log
Run
"

COMMAND_INDEPENDENT <- "4
2
TestCorrel
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.independent.log
Run
"

# the template for the cluster scheduler / job runner
SLURM_TEMPLATE <- "#!/bin/sh
#SBATCH --job-name bayestraits
#SBATCH --get-user-env
#SBATCH --output run.%%j.log
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1

BayesTraitsV3 %s %s < %s
"

# and the trees we care about...
PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

basedir <- basename(PHYLOGENIES[4])
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[4])


taxa <- read.csv(TAXAS[4])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Verbal_complexity_10", "Nominal_words_complexity_10", "case", "number", "gender", "possession", "articles", "arguments", "transitivity", "markers_arguments_non_core", "clause_v", "negation", "tense", "aspect", "mood"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))




#PHYLOGENY: Zhang et al 2020
gb <- load_data(here('data/GB.tsv'))

# what are our analyses. This is a list of pairs that we want to analyse together
ANALYSES <- list(
  "Nom_words_Verbs"= c("Nominal_words_complexity_10", "Verbal_complexity_10")
)

COMMAND_DEPENDENT <- "4
2
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.dependent.log
Run
"

COMMAND_INDEPENDENT <- "4
2
TestCorrel
ScaleTrees
Burnin 1000
Iterations 10000000
Stones 1000 100000
LogFile %s.independent.log
Run
"

# the template for the cluster scheduler / job runner
SLURM_TEMPLATE <- "#!/bin/sh
#SBATCH --job-name bayestraits
#SBATCH --get-user-env
#SBATCH --output run.%%j.log
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1

BayesTraitsV3 %s %s < %s
"

# and the trees we care about...
PHYLOGENIES <- list.files(here('data/phylogenies'), full.names=TRUE)
TAXAS <- list.files(here('data/phylogenies'), pattern="taxa.csv", full.names = TRUE, recursive = TRUE)

basedir <- basename(PHYLOGENIES[5])
dir.create(basedir)

# 2. load the trees
trees <- load_trees(PHYLOGENIES[5])


taxa <- read.csv(TAXAS[5])
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Verbal_complexity_10", "Nominal_words_complexity_10", "case", "number", "gender", "possession", "articles", "arguments", "transitivity", "markers_arguments_non_core", "clause_v", "negation", "tense", "aspect", "mood"))
gb.subset <- gb.subset[complete.cases(gb.subset),]
gb.subset <- gb[gb$taxon %in% trees[[1]]$tip.label, ]
rownames(gb.subset) <- gb.subset$taxon


for (ax in names(ANALYSES)) {
  var1 <- ANALYSES[[ax]][[1]]
  var2 <- ANALYSES[[ax]][[2]]
  
  filename <- sprintf("%s/features_%s.txt", basedir, ax)
  
  write.table(
    gb.subset[ANALYSES[[ax]]],
    file = filename,
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
  )
  
  # 5. write the command files.
  cmdi <- sprintf("%s/features_%s.independent.cmd", basedir, ax)
  cmdd <- sprintf("%s/features_%s.dependent.cmd", basedir, ax)
  writeLines(
    sprintf(COMMAND_INDEPENDENT, sprintf("features_%s", ax)), con=cmdi
  )
  writeLines(
    sprintf(COMMAND_DEPENDENT, sprintf("features_%s", ax)), con=cmdd
  )
  
  # 6. write the slurm command files for the cluster.
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdi)),
    con=sprintf("%s/features_%s.independent.job", basedir, ax)
  )
  writeLines(
    sprintf(SLURM_TEMPLATE, 'posterior.trees', basename(filename), basename(cmdd)),
    con=sprintf("%s/features_%s.dependent.job", basedir, ax)
  )
  
}
# 6. prune trees to match subset
to_remove <- setdiff(trees[[1]]$tip.label, gb.subset$taxon)
trees <- lapply(trees, drop.tip, to_remove)

# 7. write the trees
write.nexus(trees, file = sprintf("%s/posterior.trees", basedir))