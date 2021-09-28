#Visualizing the scores on the global tree
source('library.R')

# load data
gb <- load_data()
taxa <- read.csv("data/phylogenies/world/taxa.csv") #loading ASJP file (v. 17)
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]
language_glottolog <- read.csv("languages_glottolog.csv")

# merge datasets
gb.subset <- merge(gb, language_glottolog[c("Glottocode", "Family_ID", "Latitude", "Longitude")])
gb.subset <- subset(x = gb.subset, select = c("taxon", "Glottocode", "Family_ID", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#for later double checking if languages of the same family are assigned correctly, i.e. to the same clade
gb.subset %>% 
#  mutate(fam = as.numeric(if_else(Family_ID == "indo1319", "1", "0"))) %>%
  mutate(fam = as.numeric(if_else(Family_ID == "aust1307", "1", "0"))) -> gb.subset


gb.subset$fam <- ifelse(gb.subset$Family_ID == 'aust1307', 'Austronesian',
                  ifelse(gb.subset$Family_ID == 'indo1319', 'Indo-European',
                         ifelse(gb.subset$Family_ID == 'sino1245', 'Sino-Tibetan', 'other')))

#load tree
TREEFILE <- "data/phylogenies/world/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$Glottocode
tree$tip.label <- gb_geo.subset[gb_geo.subset$taxon %in% tree$tip.label, 'Glottocode']
gb_geo.subset %>% dplyr::relocate(Glottocode) -> gb_geo.subset
#checking if duplicates are present
dupes <- gb_geo.subset[duplicated(gb_geo.subset$Glottocode), ]
head(dupes)

# ggtree
# if you get warnings about this:
# Error in DataMask$new(.data, caller_env) : 
#  argument "caller_env" is missing, with no default
# 
# then run:
# > remotes::install_github("YuLab-SMU/ggtree")
# > remotes::install_github("YuLab-SMU/tidytree")
midpointNw <- (max(gb_geo.subset$Nominal_words_complexity) - min(gb_geo.subset$Nominal_words_complexity))/2 + min(gb_geo.subset$Nominal_words_complexity)
p.nc <- ggtree(tree, aes(color=Nominal_words_complexity), layout = 'circular', branch.length='none') %<+% gb_geo.subset +
  geom_tippoint(aes(color=Nominal_words_complexity)) + 
  ggtitle("a. Nominal Words") +
  scale_colour_gradient2(low="#2166AC", mid="#FFFFBF", high="#B2182B", midpoint= midpointNw, "Score") +
  theme_tree() + theme(legend.position = "bottom", plot.title = element_text(size = 14), legend.title = element_text(size = 14),  legend.text = element_text(size = 12), legend.key.size = unit(0.8, 'cm'))

midpointV <- (max(gb_geo.subset$Verbal_complexity)-min(gb_geo.subset$Verbal_complexity))/2 + min(gb_geo.subset$Verbal_complexity)
p.vc <- ggtree(tree, aes(color=Verbal_complexity), layout = 'circular', branch.length='none') %<+% gb_geo.subset +
  geom_tippoint(aes(color=Verbal_complexity)) + 
  ggtitle("b. Verbs") +
  scale_colour_gradient2(low="#2166AC", mid="#FFFFBF", high="#B2182B", midpoint= midpointV, "Score") +
  theme_tree() + theme(legend.position = "bottom", plot.title = element_text(size = 14), legend.title = element_text(size = 14),  legend.text = element_text(size = 12), legend.key.size = unit(0.8, 'cm'))

p.nc + p.vc
ggsave(filename = here("output","plot_tree_global_3_coloured.png"), p.nc + p.vc, dpi=300, dev="png")


#plotting scores on individual family trees
#PHYLOGENY Bouckaert et al 2012
gb <- load_data()
taxa <- read.csv("data/phylogenies/bouckaert_et_al2012/taxa.csv") #loading a taxa file for a corresponding phylogeny
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]

#subsetting
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/bouckaert_et_al2012/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]  # get one tree only

#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon
dupes <- gb_geo.subset[duplicated(gb_geo.subset$taxon), ]
head(dupes)


# make a palette for plotting
mycolors <- colorRampPalette(brewer.pal(6, "RdYlBu")[2: 6])

get_tree_colors <- function(vec, pal=mycolors) {
  v <- findInterval(vec, sort(vec))
  pal(length(vec))[v]
}

# plot VC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Verbal_complexity']))

# plot NWC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Nominal_words_complexity']))

# getting correct tip labels after tree was ladderized
get_tips_after_ladderize <- function(tree) {
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}

png(filename = here("output","mirror_tree_Nominal_words_Verbs_ie1.png"), width=1024, height=1024)
layout(matrix(1:3, 1, 3), widths=c(0.4, 0.2, 0.4))
# subplot 1
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Nominal_words_complexity'),
  leg.txt="Nominal words",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

# subplot 2 -- labels
plot.new(); 
plot.window(xlim=c(-0.1,0.1), ylim=get("last_plot.phylo", envir=.PlotPhyloEnv)$y.lim)
text(x=rep(0, Ntip(tree)), y=1:Ntip(tree), get_tips_after_ladderize(tree), cex=2.5)

# subplot 3
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Verbal_complexity'),
  direction="leftwards",
  leg.txt="Verbs",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

x <- dev.off()
graphics.off()



#PHYLOGENY Chang et al 2015
gb <- load_data()
taxa <- read.csv("data/phylogenies/chang_et_al2015/taxa.csv") #loading a taxa file for a corresponding phylogeny
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]

#subsetting
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/chang_et_al2015/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]  # get one tree only



#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon
dupes <- gb_geo.subset[duplicated(gb_geo.subset$taxon), ]
head(dupes)


# make a palette for plotting
mycolors <- colorRampPalette(brewer.pal(6, "RdYlBu")[2: 6])

get_tree_colors <- function(vec, pal=mycolors) {
  v <- findInterval(vec, sort(vec))
  pal(length(vec))[v]
}

# plot VC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Verbal_complexity']))

# plot NWC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Nominal_words_complexity']))

# bad ape, no banana.
get_tips_after_ladderize <- function(tree) {
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}


png(filename = here("output","mirror_tree_Nominal_words_Verbs_ie2.png"), width=1024, height=1024)

layout(matrix(1:3, 1, 3), widths=c(0.4, 0.2, 0.4))
# subplot 1
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Nominal_words_complexity'),
  leg.txt="Nominal words",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

# subplot 2 -- labels
plot.new(); 
plot.window(xlim=c(-0.1,0.1), ylim=get("last_plot.phylo", envir=.PlotPhyloEnv)$y.lim)
text(x=rep(0, Ntip(tree)), y=1:Ntip(tree), get_tips_after_ladderize(tree), cex=2.5)

# subplot 3
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Verbal_complexity'),
  direction="leftwards",
  leg.txt="Verbs",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

x <- dev.off()
graphics.off()




#PHYLOGENY Gray et al 2009
gb <- load_data()
taxa <- read.csv("data/phylogenies/gray_et_al2009/taxa.csv") #loading a taxa file for a corresponding phylogeny
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]

#subsetting
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/gray_et_al2009/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]  # get one tree only



#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon
dupes <- gb_geo.subset[duplicated(gb_geo.subset$taxon), ]
head(dupes)


# make a palette for plotting
mycolors <- colorRampPalette(brewer.pal(6, "RdYlBu")[2: 6])

get_tree_colors <- function(vec, pal=mycolors) {
  v <- findInterval(vec, sort(vec))
  pal(length(vec))[v]
}

# plot VC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Verbal_complexity']))

# plot NWC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Nominal_words_complexity']))

get_tips_after_ladderize <- function(tree) {
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}


png(filename = here("output","mirror_tree_Nominal_words_Verbs_a.png"), width=1024, height=1024)

layout(matrix(1:3, 1, 3), widths=c(0.4, 0.2, 0.4))
# subplot 1
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Nominal_words_complexity'),
  leg.txt="Nominal words",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

# subplot 2 -- labels
plot.new(); 
plot.window(xlim=c(-0.1,0.1), ylim=get("last_plot.phylo", envir=.PlotPhyloEnv)$y.lim)
text(x=rep(0, Ntip(tree)), y=1:Ntip(tree), get_tips_after_ladderize(tree), cex=2.3)

# subplot 3
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Verbal_complexity'),
  direction="leftwards",
  leg.txt="Verbs",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

x <- dev.off()
graphics.off()



#Phylogeny: Zhang et al 2020
gb <- load_data()
taxa <- read.csv("data/phylogenies/zhang_et_al2020/taxa.csv") #loading a taxa file for a corresponding phylogeny
gb$taxon <- taxa$taxon[match(gb$Glottocode, taxa$glottocode)]

#subsetting
gb.subset <- subset(x = gb, select = c("taxon", "Glottocode", "Nominal_words_complexity", "Verbal_complexity"))
gb.subset <- gb.subset[complete.cases(gb.subset),]

#load tree
TREEFILE <- "data/phylogenies/zhang_et_al2020/"
tree <- sample(load_trees(TREEFILE), 1)[[1]]  # get one tree only


#merge data with tree and set up tree
gb_geo.subset <- gb.subset[gb.subset$taxon %in% tree$tip.label, ]
to_remove <- setdiff(tree$tip.label, gb_geo.subset$taxon)
tree <- drop.tip(tree, to_remove)

gb_geo.subset<-gb_geo.subset[match(tree$tip.label, gb_geo.subset$taxon),]
rownames(gb_geo.subset) <- gb_geo.subset$taxon
dupes <- gb_geo.subset[duplicated(gb_geo.subset$taxon), ]
head(dupes)


# make a palette for plotting
mycolors <- colorRampPalette(brewer.pal(6, "RdYlBu")[2: 6])

get_tree_colors <- function(vec, pal=mycolors) {
  v <- findInterval(vec, sort(vec))
  pal(length(vec))[v]
}

# plot VC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Verbal_complexity']))

# plot NWC
base::plot(tree, cex = 0.5, align.tip.label = TRUE, label.offset=0.1, font=1)
tiplabels(pch=22, bg=get_tree_colors(gb_geo.subset[tree$tip.label, 'Nominal_words_complexity']))

# bad ape, no banana.
get_tips_after_ladderize <- function(tree) {
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}


png(filename = here("output","mirror_tree_Nominal_words_Verbs_st.png"), width=1024, height=1024)

layout(matrix(1:3, 1, 3), widths=c(0.4, 0.2, 0.4))
# subplot 1
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Nominal_words_complexity'),
  leg.txt="Nominal words",
  lwd=2, outline=FALSE,
  ftype="off", sig=2, fsize=c(1.9), res=300
)

# subplot 2 -- labels
plot.new(); 
plot.window(xlim=c(-0.1,0.1), ylim=get("last_plot.phylo", envir=.PlotPhyloEnv)$y.lim)
text(x=rep(0, Ntip(tree)), y=1:Ntip(tree), get_tips_after_ladderize(tree), cex=2.5)

# subplot 3
contMap(
  tree, get_trait_vector(tree, gb_geo.subset, 'Verbal_complexity'),
  direction="leftwards",
  leg.txt="Verbs",
  lwd=2, outline=FALSE,
  ftype="off", fsize=c(1.9), res=300
)

x <- dev.off()
graphics.off()
