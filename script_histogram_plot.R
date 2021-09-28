#plotting a histogram
source('library.R')
source('preparing_world_subset.R')


compare <- gb_world_subset %>%
  pivot_longer(c(Nominal_words_complexity, Verbal_complexity), names_to="complexity_type", values_to="Complexity_score")

h1 <- ggplot(compare, aes(x=Complexity_score, fill=complexity_type, color=complexity_type)) +
  geom_histogram(position="identity", alpha=0.5, binwidth=0.04) + 
  scale_fill_manual(name=" ", values = c("#00AFBB", "#E7B800"), labels = c("Nominal words", "Verbs")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"), guide="none") + ylim(0, 60) + xlim(0, 1.01) + theme_classic() + 
  theme(legend.position = c(0.3, 0.95), legend.direction = "vertical") +
  xlab("Score") +  ylab("languages") + ylim(0, 40) +
  theme(axis.text = element_text(size = 17)) + theme(axis.title = element_text(size = 18)) +
  theme(legend.text = element_text(size = 17)) + theme(legend.title = element_text(size = 13))
h1

ggsave(h1, filename = here("output","two_overlapping_histograms.png"), dpi=300, dev="png")
