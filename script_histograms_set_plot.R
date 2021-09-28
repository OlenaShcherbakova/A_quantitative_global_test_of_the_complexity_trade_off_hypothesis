#plotting variation in each feature group of nominal words and verbs metrics

source('library.R')
source('preparing_world_subset.R')

gb_world_subset <- gb_world_subset %>%
  mutate(across(5:18, type.convert))  %>%
  mutate(across(5:18, round, 2))


h1 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$case)),], aes(case)) + geom_histogram(stat="count", color = "#00AFBB", fill = "#00AFBB", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Case") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
h1

h2 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$number)),], aes(number)) + geom_histogram(stat="count", color = "#00AFBB", fill = "#00AFBB", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Number") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.12, 0.25, 0.38, 0.5, 0.62, 0.75, 0.87, 1))
h2

h3 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$gender)),], aes(gender)) + geom_histogram(stat="count", color = "#00AFBB", fill = "#00AFBB", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Gender") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.17, 0.33, 0.5, 0.67, 0.83, 1))
h3

h4 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$possession)),], aes(possession)) + geom_histogram(stat="count", color = "#00AFBB", fill = "#00AFBB", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Possession") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.5, 1))
h4

h5 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$articles)),], aes(articles)) + geom_histogram(stat="count", color = "#00AFBB", fill = "#00AFBB", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Articles") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.5, 1))
h5

h6 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$arguments)),], aes(arguments)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Arguments") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20))+ scale_x_continuous(breaks = c(0, 0.33, 0.67, 1))
h6

h7 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$transitivity)),], aes(transitivity)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Transitivity") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.17, 0.33, 0.5, 0.67, 0.83, 1))
h7

h8 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$markers_arguments_non_core)),], aes(markers_arguments_non_core)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Non-core arguments") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
h8

h9 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$clause_v)),], aes(clause_v)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Clause-related") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 0.5, 1))
h9

h10 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$negation)),], aes(negation)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Negation") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 1))
h10

h11 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$tense)),], aes(tense)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Tense") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 1))
h11

h12 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$aspect)),], aes(aspect)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Aspect") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 1))
h12

h13 <- ggplot(gb_world_subset[!(is.na(gb_world_subset$mood)),], aes(mood)) + geom_histogram(stat="count", color = "#E7B800", fill = "#E7B800", alpha=0.4)+ ylim(0, 360) + xlim(0, 1) + xlab("Mood") + ylab("languages") + theme_classic() + theme(text = element_text(size = 20)) + scale_x_continuous(breaks = c(0, 1))
h13

ggsave(ggarrange(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, ncol=4, nrow=4, font.label = list(size = 5, face = "plain", color ="black")), width=24, height=15, filename = here("output","domains_histograms_set.png"), dpi=300, dev="png")