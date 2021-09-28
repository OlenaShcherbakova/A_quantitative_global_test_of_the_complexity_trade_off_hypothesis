#measuring phylogenetic signal

#global tree
source("library.R")
gb <- load_data()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/world/taxa.csv") 

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/world/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)


#measuring phylogenetic signal (global):
tree.Nominal_words_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Nominal_words_complexity), 'Glottocode']))
Nominal_words_complexity <- get_trait_vector(tree.Nominal_words_complexity, gb_geo.subset, 'Nominal_words_complexity')
physig_Nw_world_l <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="lambda", test=TRUE)

lambda_Nw_world_l <- physig_Nw_world_l[1][["lambda"]]
LR_Nw_world_l <- 2*(physig_Nw_world_l$logL-physig_Nw_world_l$logL0) #performing likelihood ratio test
P_lambda_Nw_world_l <- physig_Nw_world_l$P

physig_Nw_world_K <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="K", test=TRUE)
K_Nw_world_K <- physig_Nw_world_K[1][["K"]]
P_Nw_world_K <- physig_Nw_world_K$P

tree.Verbal_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Verbal_complexity), 'Glottocode']))
Verbal_complexity <- get_trait_vector(tree.Verbal_complexity, gb_geo.subset, 'Verbal_complexity')
physig_V_world_l <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="lambda", test=TRUE)
lambda_V_world_l <- physig_V_world_l[1][["lambda"]]
LR_V_world_l <- 2*(physig_V_world_l$logL-physig_V_world_l$logL0) #performing likelihood-ratio test
P_lambda_V_world_l <- physig_V_world_l$P

physig_V_world_K <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="K", test=TRUE)
K_V_world_K <- physig_V_world_K[1][["K"]]
P_V_world_K <- physig_V_world_K$P

Nw_world <- c(physig_Nw_world_l$logL, physig_Nw_world_l$logL0, LR_Nw_world_l, P_lambda_Nw_world_l, K_Nw_world_K, P_Nw_world_K)
V_world <- c(physig_V_world_l$logL, physig_V_world_l$logL0, LR_V_world_l, P_lambda_V_world_l, K_V_world_K, P_V_world_K)



#Indo-European tree (1): Bouckaert et al 2012
source("library.R")
gb <- load_data()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/bouckaert_et_al2012/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/bouckaert_et_al2012/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)

#measuring phylogenetic signal (ie (1:bouckaert et al 2012)):
tree.Nominal_words_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Nominal_words_complexity), 'Glottocode']))
Nominal_words_complexity <- get_trait_vector(tree.Nominal_words_complexity, gb_geo.subset, 'Nominal_words_complexity')
physig_Nw_ie1_l <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="lambda", test=TRUE)
lambda_Nw_ie1_l <- physig_Nw_ie1_l[1][["lambda"]]
LR_Nw_ie1_l <- 2*(physig_Nw_ie1_l$logL-physig_Nw_ie1_l$logL0) #performing likelihood ratio test
P_lambda_Nw_ie1_l <- physig_Nw_ie1_l$P

physig_Nw_ie1_K <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="K", test=TRUE)
K_Nw_ie1_K <- physig_Nw_ie1_K[1][["K"]]
P_Nw_ie1_K <- physig_Nw_ie1_K$P

tree.Verbal_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Verbal_complexity), 'Glottocode']))
Verbal_complexity <- get_trait_vector(tree.Verbal_complexity, gb_geo.subset, 'Verbal_complexity')
physig_V_ie1_l <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="lambda", test=TRUE)
lambda_V_ie1_l <- physig_V_ie1_l[1][["lambda"]]
LR_V_ie1_l <- 2*(physig_V_ie1_l$logL-physig_V_ie1_l$logL0) #performing likelihood ratio test
P_lambda_V_ie1_l <- physig_V_ie1_l$P

physig_V_ie1_K <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="K", test=TRUE)
K_V_ie1_K <- physig_V_ie1_K[1][["K"]] 
P_V_ie1_K <- physig_V_ie1_K$P

Nw_ie1 <- c(physig_Nw_ie1_l$logL, physig_Nw_ie1_l$logL0, LR_Nw_ie1_l, P_lambda_Nw_ie1_l, K_Nw_ie1_K, P_Nw_ie1_K)
V_ie1 <- c(physig_V_ie1_l$logL, physig_V_ie1_l$logL0, LR_V_ie1_l, P_lambda_V_ie1_l, K_V_ie1_K, P_V_ie1_K)




#Indo-European tree (2): Chang et al 2015
source("library.R")
gb <- load_data()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/chang_et_al2015/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/chang_et_al2015/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)
#measuring phylogenetic signal (ie (2:chang et al 2015)):
tree.Nominal_words_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Nominal_words_complexity), 'Glottocode']))
Nominal_words_complexity <- get_trait_vector(tree.Nominal_words_complexity, gb_geo.subset, 'Nominal_words_complexity')
physig_Nw_ie2_l <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="lambda", test=TRUE)
lambda_Nw_ie2_l <- physig_Nw_ie2_l[1][["lambda"]]
LR_Nw_ie2_l <- 2*(physig_Nw_ie2_l$logL-physig_Nw_ie2_l$logL0) #performing likelihood ratio test
P_lambda_Nw_ie2_l <- physig_Nw_ie2_l$P

physig_Nw_ie2_K <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="K", test=TRUE)
K_Nw_ie2_K <- physig_Nw_ie2_K[1][["K"]] 
P_Nw_ie2_K <- physig_Nw_ie2_K$P

tree.Verbal_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Verbal_complexity), 'Glottocode']))
Verbal_complexity <- get_trait_vector(tree.Verbal_complexity, gb_geo.subset, 'Verbal_complexity')
physig_V_ie2_l <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="lambda", test=TRUE)
lambda_V_ie2_l <- physig_V_ie2_l[1][["lambda"]]
LR_V_ie2_l <- 2*(physig_V_ie2_l$logL-physig_V_ie2_l$logL0) #performing likelihood ratio test
P_lambda_V_ie2_l <- physig_V_ie2_l$P

physig_V_ie2_K <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="K", test=TRUE)
K_V_ie2_K <- physig_V_ie2_K[1][["K"]] 
P_V_ie2_K <- physig_V_ie2_K$P

Nw_ie2 <- c(physig_Nw_ie2_l$logL, physig_Nw_ie2_l$logL0, LR_Nw_ie2_l, P_lambda_Nw_ie2_l, K_Nw_ie2_K, P_Nw_ie2_K)
V_ie2 <- c(physig_V_ie2_l$logL, physig_V_ie2_l$logL0, LR_V_ie2_l, P_lambda_V_ie2_l, K_V_ie2_K, P_V_ie2_K)



#Austronesian tree
source("library.R")
gb <- load_data()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/gray_et_al2009/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/gray_et_al2009/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)

#measuring phylogenetic signal (Austronesian):
tree.Nominal_words_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Nominal_words_complexity), 'Glottocode']))
Nominal_words_complexity <- get_trait_vector(tree.Nominal_words_complexity, gb_geo.subset, 'Nominal_words_complexity')
physig_Nw_a_l <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="lambda", test=TRUE)
lambda_Nw_a_l <- physig_Nw_a_l[1][["lambda"]]
LR_Nw_a_l <- 2*(physig_Nw_a_l$logL-physig_Nw_a_l$logL0) #performing likelihood ratio test
P_lambda_Nw_a_l <- physig_Nw_a_l$P

physig_Nw_a_K <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="K", test=TRUE)
K_Nw_a_K <- physig_Nw_a_K[1][["K"]] 
P_Nw_a_K <- physig_Nw_a_K$P

tree.Verbal_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Verbal_complexity), 'Glottocode']))
Verbal_complexity <- get_trait_vector(tree.Verbal_complexity, gb_geo.subset, 'Verbal_complexity')
physig_V_a_l <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="lambda", test=TRUE)
lambda_V_a_l <- physig_V_a_l[1][["lambda"]]
LR_V_a_l <- 2*(physig_V_a_l$logL-physig_V_a_l$logL0) #performing likelihood ratio test
P_lambda_V_a_l <- physig_V_a_l$P

physig_V_a_K <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="K", test=TRUE)
K_V_a_K <- physig_V_a_K[1][["K"]] 
P_V_a_K <- physig_V_a_K$P

Nw_a <- c(physig_Nw_a_l$logL, physig_Nw_a_l$logL0, LR_Nw_a_l, P_lambda_Nw_a_l, K_Nw_a_K, P_Nw_a_K)
V_a <- c(physig_V_a_l$logL, physig_V_a_l$logL0, LR_V_a_l, P_lambda_V_a_l, K_V_a_K, P_V_a_K)




#Sino-Tibetan tree
source("library.R")
gb <- load_data()

#loading ASJP file (v. 17) for world tree or a taxa file for a corresponding phylogeny
taxa <- read.csv("data/phylogenies/zhang_et_al2020/taxa.csv")

gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/zhang_et_al2020/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon

#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)

#measuring phylogenetic signal (Sino-Tibetan):
tree.Nominal_words_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Nominal_words_complexity), 'Glottocode']))
Nominal_words_complexity <- get_trait_vector(tree.Nominal_words_complexity, gb_geo.subset, 'Nominal_words_complexity')
physig_Nw_st_l <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="lambda", test=TRUE)
lambda_Nw_st_l <- physig_Nw_st_l[1][["lambda"]]
LR_Nw_st_l <- 2*(physig_Nw_st_l$logL-physig_Nw_st_l$logL0) #performing likelihood ratio test
P_lambda_Nw_st_l <- physig_Nw_st_l$P

physig_Nw_st_K <- phytools::phylosig(tree.Nominal_words_complexity, Nominal_words_complexity, method="K", test=TRUE)
K_Nw_st_K <- physig_Nw_st_K[1][["K"]] 
P_Nw_st_K <- physig_Nw_st_K$P

tree.Verbal_complexity <- drop.tip(tree, as.vector(gb_geo.subset[is.na(gb_geo.subset$Verbal_complexity), 'Glottocode']))
Verbal_complexity <- get_trait_vector(tree.Verbal_complexity, gb_geo.subset, 'Verbal_complexity')
physig_V_st_l <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="lambda", test=TRUE)
lambda_V_st_l <- physig_V_st_l[1][["lambda"]]
LR_V_st_l <- 2*(physig_V_st_l$logL-physig_V_st_l$logL0) #performing likelihood ratio test
P_lambda_V_st_l <- physig_V_st_l$P

physig_V_st_K <- phytools::phylosig(tree.Verbal_complexity, Verbal_complexity, method="K", test=TRUE)
K_V_st_K <- physig_V_st_K[1][["K"]] 
P_V_st_K <- physig_V_st_K$P

Nw_st <- c(physig_Nw_st_l$logL, physig_Nw_st_l$logL0, LR_Nw_st_l, P_lambda_Nw_st_l, K_Nw_st_K, P_Nw_st_K)
V_st <- c(physig_V_st_l$logL, physig_V_st_l$logL0, LR_V_st_l, P_lambda_V_st_l, K_V_st_K, P_V_st_K)



#Making a table out of two measures of phylogenetic signal
physig <- as.data.frame(rbind(Nw_world, V_world, Nw_a, V_a, Nw_st, V_st, Nw_ie1, V_ie1, Nw_ie2, V_ie2))
colnames(physig) <- c("logL", "logL0", "LR (lambda)", "p-value (lambda)", "K", "p-value (K)")
physig <- round(physig, digits=2)
rownames(physig) <- c("Nominal words (World)", "Verbs (World)", "Nominal words (Austronesian)", "Verbs (Austronesian)", "Nominal words (Sino-Tibetan)", "Verbs (Sino-Tibetan)", "Nominal words (Indo-European 1)", "Verbs (Indo-European 1)", "Nominal words (Indo-European 2)", "Verbs (Indo-European 2)")
phylogeny <- as.data.frame(c(rep(c("World"), times=2), rep(c("Austronesian"), times=2), rep(c("Sino-Tibetan"), times=2), rep(c("Indo-European 1"), times=2), rep(c("Indo-European 2"), times=2)))
colnames(phylogeny) <- "Phylogeny"
gram_coding <- as.data.frame(c(rep(c("Nominal words", "Verbs"), times=5)))
colnames(gram_coding) <- "Grammatical coding"
physig <- cbind(phylogeny, gram_coding, physig)
write.csv(physig, file=here("output", "Table_C.csv"), row.names = FALSE)


