source('library.R')

languages_glottolog <- read.csv("data/languages_glottolog.csv", header=TRUE)
gb <- load_data()
data <- merge(gb, languages_glottolog)

#Only taking languages from the separate phylogenies which were used for the analyses
taxons <- list(
  'Indo-European' = "data/phylogenies/bouckaert_et_al2012/taxa.csv",
  'Indo-European' = "data/phylogenies/chang_et_al2015/taxa.csv",
  'Austronesian' = "data/phylogenies/gray_et_al2009/taxa.csv",
  'Sino-Tibetan' = "data/phylogenies/zhang_et_al2020/taxa.csv"
)

taxons.df <- data.frame(Family=NULL, Language=NULL)
for (t in names(taxons)) {
  taxons.df <- rbind(
    taxons.df, data.frame(
      Family=t,
      Glottocode=read.csv(taxons[[t]], header=TRUE)$glottocode
    )
  )
}

taxons.df <- unique( taxons.df )


taxons <- taxons.df$Family
names(taxons) <- taxons.df$Glottocode

languages <- c(
  data[data$Verbal_complexity <= 0.25, 'Glottocode'],
  data[data$Verbal_complexity >= 0.93, 'Glottocode'],
  data[data$Nominal_words_complexity >= 0.68, 'Glottocode']
)


data$Subset <- as.vector(taxons[data$Glottocode])
data[is.na(data$Subset), 'Subset'] <- 'World'

palette <- c(RColorBrewer::brewer.pal(3, 'Dark2'), "gray40")

P <- ggplot(data, aes(Verbal_complexity, Nominal_words_complexity, color=Subset)) +
  geom_point() +
  geom_text_repel(
    data=data[data$Glottocode %in% languages,], aes(label=Name)
  ) +
  scale_color_manual(values=palette) +
  xlab("Verbs") + xlim(0, 1) +
  ylab("Nominal words") + ylim(0, 1) +
  theme_bw() +  guides(color="none")
P

# subplot
PA <- ggplot(data , aes(Verbal_complexity, Nominal_words_complexity, color=Subset, group=Subset)) +
  geom_point(color="gray90") + 
  geom_point(data=data[data$Subset =='Austronesian',], color=palette[[1]]) +
  xlab("") + 
  xlim(0, 1) +
  ylab("") + 
  ylim(0, 1) +
  theme_bw(base_size=10) +
  ggtitle("Austronesian") +
  guides(color="none")
PA

PI <- ggplot(data , aes(Verbal_complexity, Nominal_words_complexity, color=Subset, group=Subset)) +
  geom_point(color="gray90") + 
  geom_point(data=data[data$Subset == 'Indo-European',], color=palette[[2]]) +
  xlab("") + 
  xlim(0, 1) +
  ylab("") + 
  ylim(0, 1) +
  theme_bw(base_size=10) +
  ggtitle("Indo-European") +
  guides(color="none")
PI


PS <- ggplot(data, aes(Verbal_complexity, Nominal_words_complexity, color=Subset, group=Subset)) +
  geom_point(color="gray90") + 
  geom_point(data=data[data$Subset == 'Sino-Tibetan',], color=palette[[3]]) +
  xlab("") + 
  xlim(0, 1) +
  ylab("") + 
  ylim(0, 1) +
  theme_bw(base_size=10) +
  ggtitle("Sino-Tibetan") +
  guides(color="none")
PS

(P / (PA | PS  | PI)) + plot_layout(heights = c(3, 1))
ggsave(filename=here("output", "scatter_plot_302_L.png"), dev="png", width=7, height=5, dpi=300)
