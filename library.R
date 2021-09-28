#library script

library(ape)
library(tidyverse)
library(geiger)
library(ggplot2)
library(dplyr)
library(phytools)
library(here)
library(ggpubr)
library(viridis)
library(mapdata)
library(maps)
library(maptools)
library(grDevices)
library(rethinking)
library(bayestraitr)
library(ggtree)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
source('world_tree_transformation.R')

#script written by Simon Greenhill and Olena Shcherbakova

# helper function to convert and load trait data into a named vector for
# plotting.
get_trait_vector <- function(tree, data, variable) {
  x <- data[tree$tip.label, variable]
  x[is.na(x)] <- 0  # set NA's to zero to enable plotting.
  names(x) <- tree$tip.label
  x
}

Nominal_words_complexity <- list(
  #cases in nouns and pronouns:
  c("GB070", "GB072", "GB071", "GB073"),
  #number in nouns
  c("singular", "dual", "plural", "trial", "paucal", "GB031", "GB185", "GB186"),
  #gender in nouns
  c("gender_nouns", "gender_pronouns", "GB170", "GB198", "GB172", "GB171", "GB057"),
  #possession marked on nouns:
  c("possession_on_possessor", "possession_on_possessed"),
  #articles
  c("GB020", "GB021")
)


Verbal_complexity <- list(
  #marking of arguments:
  c("S_arg", "A_arg", "P_arg"),
  #transitivity-related features:
  c("GB113", "GB124", "GB149", "GB155", "passive", "antipassive"),
  #negation:
  c("negation"),
  #tenses:
  c("tense"),
  #aspect:
  c("aspect"),
  #mood:
  c("mood"),
  #markers_arguments_non_core:
  c("GB103", "GB104", "GB108", "reflexivity", "reciprocity"),
  #syntactic features:
  c("GB151", "GB152")
)



# A function to create a metric from a dataset and a given set of recodings.
create_metric <- function(df, recodings) {
  # apply rowMeans to all the sets of variables in the recodings list
  scores <- sapply(
    recodings,
    function(df, varlist) { rowMeans(df[varlist], na.rm = TRUE) },
    df=df, simplify=TRUE
  )
  rowMeans(scores, na.rm = TRUE)
}

# Function to load trees from DPLACE-data repository
# can be used either with or without renameto = 'glottocode'
load_trees <- function(dirname, type='posterior', mappingfile='taxa.csv', renameto=NA) {
  # check file type
  if (type == 'summary') {
    treefile <- file.path(dirname, 'summary.trees')
  }
  else if (type == 'posterior') {
    treefile <- file.path(dirname, 'posterior.trees')
  } else {
    stop(paste("Unknown Tree Type:", type))
  }
  
  # check file exists
  if (file.exists(treefile) == FALSE) {
    stop(paste("Invalid file:", treefile))
  }
  
  trees <- ape::read.nexus(treefile)
  if (class(trees) == 'phylo') { trees <- c(trees) ; class(trees) <- 'multiPhylo' }
  
  # make full path if just given taxa.csv
  if (mappingfile == 'taxa.csv') { mappingfile <- file.path(dirname, mappingfile) }
  
  if (file.exists(mappingfile) & is.na(renameto) == FALSE) {
    mapping <- read.csv(mappingfile, header = TRUE, stringsAsFactors = FALSE, na.string="")
    
    # check the required columns exist
    if ('taxon' %in% colnames(mapping) == FALSE) stop(paste('column `taxon` not in', mappingfile))
    
    if (renameto %in% colnames(mapping) == FALSE) stop(paste('colname', renameto, 'not in', mappingfile))
    
    trees <- ape::.uncompressTipLabel(trees)
    
    for (i in 1:length(trees)){
      # remove tips not in `taxon` mapping
      missing <- trees[[i]]$tip.label[trees[[i]]$tip.label %in% mapping[['taxon']] == FALSE]
      if (length(missing) > 0) {
        trees[[i]] <- ape::drop.tip(trees[[i]], missing)
      }
      
      # remove tips not in `renameto` mapping
      missing <- mapping[is.na(mapping[[renameto]]), 'taxon']
      if (length(missing) > 0) {
        trees[[i]] <- ape::drop.tip(trees[[i]], missing)
      }
      
      # handle duplicate rename tips
      dupes <- mapping[duplicated(mapping[[renameto]], incomparables=NA), ]
      if (nrow(dupes)) {
        warning(paste("Removing ", nrow(dupes), "tips that will be duplicated after rename:", paste(dupes[['taxon']], collapse=", ")))
        trees[[i]] <- ape::drop.tip(trees[[i]], dupes[['taxon']])
      }
      
      # rename tips
      matches <- match(trees[[i]]$tip.label, mapping[['taxon']])
      trees[[i]]$tip.label <- mapping[matches, renameto]
    }
    trees <- ape::.compressTipLabel(trees, ref=mapping[matches, renameto])
  }
  trees
}

load_data <- function(filename="data/GB.tsv") {
  grambank <- read.csv(filename, header = TRUE, sep = '\t', stringsAsFactors=FALSE)
  colnames(grambank)[colnames(grambank)=="Language_ID"] <- "Glottocode"  # rename column
  
  #removing languages that have NAs for Grambank features included in the metric
  grambank <- subset(x=grambank, select=c("Glottocode", 
                                          "GB070", "GB072", "GB071", "GB073", 
                                          "GB042", "GB316", 
                                          "GB043", "GB317", 
                                          "GB044", "GB318", 
                                          "GB165", "GB319", 
                                          "GB166", "GB320", 
                                          "GB057",
                                          "GB184", "GB031", "GB185", "GB186", 
                                          "GB051", "GB052", "GB053", "GB054", "GB192", 
                                          "GB170", "GB196", "GB197", "GB030", 
                                          "GB198", "GB172", "GB171", 
                                          "GB430", "GB432", "GB431", "GB433", 
                                          "GB020", "GB021", 
                                          "GB089", "GB090", "GB091", "GB092", "GB093", "GB094",
                                          "GB113", "GB124", "GB149", "GB155", "GB147", "GB302", "GB148", "GB303",
                                          "GB107", "GB298", "GB299",
                                          "GB082", "GB083", "GB084", "GB121", "GB521",
                                          "GB086", "GB120", "GB520",
                                          "GB312", "GB119", "GB519",
                                          "GB103", "GB104", "GB108", "GB114", "GB305", "GB115", "GB306",
                                          "GB151", "GB152"))
  grambank <- na.omit(grambank)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB042[i]), as.numeric(grambank$GB316[i])), na.rm = T)
    if(summ > 0 ){grambank$singular[i] <- 1}
    else(grambank$singular[i] <- 0)
  }
  grambank$singular <- as.numeric(grambank$singular)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB043[i]), as.numeric(grambank$GB317[i])), na.rm = T)
    if(summ > 0 ){grambank$dual[i] <- 1}
    else(grambank$dual[i] <- 0)
  }
  grambank$dual <- as.numeric(grambank$dual)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB044[i]), as.numeric(grambank$GB318[i])), na.rm = T)
    if(summ > 0 ){grambank$plural[i] <- 1}
    else(grambank$plural[i] <- 0)
  }
  grambank$plural <- as.numeric(grambank$plural)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB165[i]), as.numeric(grambank$GB319[i])), na.rm = T)
    if(summ > 0 ){grambank$trial[i] <- 1}
    else(grambank$trial[i] <- 0)
  }
  grambank$trial <- as.numeric(grambank$trial)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB166[i]), as.numeric(grambank$GB320[i])), na.rm = T)
    if(summ > 0 ){grambank$paucal[i] <- 1}
    else(grambank$paucal[i] <- 0)
  }
  grambank$paucal <- as.numeric(grambank$paucal)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB431[i]), as.numeric(grambank$GB433[i])), na.rm = T)
    if(summ > 0 ){grambank$possession_on_possessed[i] <- 1}
    else(grambank$possession_on_possessed[i] <- 0)
  }
  grambank$possession_on_possessed <- as.numeric(grambank$possession_on_possessed)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB430[i]), as.numeric(grambank$GB432[i])), na.rm = T)
    if(summ > 0 ){grambank$possession_on_possessor[i] <- 1}
    else(grambank$possession_on_possessor[i] <- 0)
  }
  grambank$possession_on_possessor <- as.numeric(grambank$possession_on_possessor)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB051[i]), as.numeric(grambank$GB052[i]), as.numeric(grambank$GB053[i]), as.numeric(grambank$GB054[i]), as.numeric(grambank$GB192[i])), na.rm = T)
    if(summ > 0 ){grambank$gender_nouns[i] <- 1}
    else(grambank$gender_nouns[i] <- 0)
  }
  grambank$gender_nouns <- as.numeric(grambank$gender_nouns)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB196[i]), as.numeric(grambank$GB197[i]), as.numeric(grambank$GB030[i])), na.rm = T)
    if(summ > 0 ){grambank$gender_pronouns[i] <- 1}
    else(grambank$gender_pronouns[i] <- 0)
  }
  grambank$gender_pronouns <- as.numeric(grambank$gender_pronouns)
  
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB089[i]), as.numeric(grambank$GB090[i])), na.rm = T)
    if(summ > 0 ){grambank$S_arg[i] <- 1}
    else(grambank$S_arg[i] <- 0)
  }
  grambank$S_arg <- as.numeric(grambank$S_arg)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB091[i]), as.numeric(grambank$GB092[i])), na.rm = T)
    if(summ > 0 ){grambank$A_arg[i] <- 1}
    else(grambank$A_arg[i] <- 0)
  }
  grambank$A_arg <- as.numeric(grambank$A_arg)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB093[i]), as.numeric(grambank$GB094[i])), na.rm = T)
    if(summ > 0 ){grambank$P_arg[i] <- 1}
    else(grambank$P_arg[i] <- 0)
  }
  grambank$P_arg <- as.numeric(grambank$P_arg)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB147[i]), as.numeric(grambank$GB302[i])), na.rm = T)
    if(summ > 0 ){grambank$passive[i] <- 1}
    else(grambank$passive[i] <- 0)
  }
  grambank$passive <- as.numeric(grambank$passive)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB148[i]), as.numeric(grambank$GB303[i])), na.rm = T)
    if(summ > 0 ){grambank$antipassive[i] <- 1}
    else(grambank$antipassive[i] <- 0)
  }
  grambank$antipassive <- as.numeric(grambank$antipassive)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB107[i]), as.numeric(grambank$GB298[i]), as.numeric(grambank$GB299[i])), na.rm = T)
    if(summ > 0 ){grambank$negation[i] <- 1}
    else(grambank$negation[i] <- 0)
  }
  grambank$negation <- as.numeric(grambank$negation)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB082[i]), as.numeric(grambank$GB083[i]), as.numeric(grambank$GB084[i]), as.numeric(grambank$GB121[i]), as.numeric(grambank$GB521[i])), na.rm = T)
    if(summ > 0 ){grambank$tense[i] <- 1}
    else(grambank$tense[i] <- 0)
  }
  grambank$tense <- as.numeric(grambank$tense)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB114[i]), as.numeric(grambank$GB305[i])), na.rm = T)
    if(summ > 0 ){grambank$reflexivity[i] <- 1}
    else(grambank$reflexivity[i] <- 0)
  }
  grambank$reflexivity <- as.numeric(grambank$reflexivity)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB115[i]), as.numeric(grambank$GB306[i])), na.rm = T)
    if(summ > 0 ){grambank$reciprocity[i] <- 1}
    else(grambank$reciprocity[i] <- 0)
  }
  grambank$reciprocity <- as.numeric(grambank$reciprocity)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB086[i]), as.numeric(grambank$GB120[i]), as.numeric(grambank$GB520[i])), na.rm = T)
    if(summ > 0 ){grambank$aspect[i] <- 1}
    else(grambank$aspect[i] <- 0)
  }
  grambank$aspect <- as.numeric(grambank$aspect)
  
  for(i in 1:nrow(grambank)){
    summ <- sum(c(as.numeric(grambank$GB312[i]), as.numeric(grambank$GB119[i]), as.numeric(grambank$GB519[i])), na.rm = T)
    if(summ > 0 ){grambank$mood[i] <- 1}
    else(grambank$mood[i] <- 0)
  }
  grambank$mood <- as.numeric(grambank$mood)
  
  # setup metrics..
  
  #Metric 1. Verbal domain
  grambank$Verbal_complexity <- create_metric(grambank, Verbal_complexity)
  
  #Metric 2. Nominal words domain
  grambank$Nominal_words_complexity <- create_metric(grambank, Nominal_words_complexity)
  
  #multiplying scores by 10 (otherwise BayesTraits complains about the values)
  grambank$Verbal_complexity_10 <- grambank$Verbal_complexity * 10
  grambank$Nominal_words_complexity_10 <- grambank$Nominal_words_complexity * 10
  
  grambank$case <- rowMeans(grambank[,c("GB070", "GB071", "GB072", "GB073")])
  grambank$number <- rowMeans(grambank[,c("singular", "dual", "plural", "trial", "paucal", "GB031", "GB185", "GB186", "GB057")])
  grambank$gender <- rowMeans(grambank[,c("gender_nouns", "gender_pronouns", "GB170", "GB198", "GB172", "GB171")])
  grambank$possession <- rowMeans(grambank[,c("possession_on_possessor", "possession_on_possessed")])
  grambank$articles <- rowMeans(grambank[,c("GB020", "GB021")])
  
  grambank$arguments <- rowMeans(grambank[,c("S_arg", "A_arg", "P_arg")])
  grambank$transitivity <- rowMeans(grambank[,c("GB113", "GB124", "GB149", "GB155", "passive", "antipassive")])
  grambank$markers_arguments_non_core <- rowMeans(grambank[,c("GB103", "GB104", "GB108", "reflexivity", "reciprocity")])
  grambank$clause_v <- rowMeans(grambank[,c("GB151", "GB152")])
  
  sample <- subset(x = grambank, select = c("Glottocode", "Verbal_complexity_10", "Nominal_words_complexity_10", "Verbal_complexity", "Nominal_words_complexity", "case", "number", "gender", "possession", "articles", "arguments", "transitivity", "markers_arguments_non_core", "clause_v", "negation", "tense", "aspect", "mood"))
  
  rownames(sample) <- sample$Glottocode
  
  sample
  
}
