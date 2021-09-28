#visualizing complexity on the world map
source('library.R')
source('preparing_world_subset.R')

world <- map_data("world", wrap=c(-25,335), ylim=c(-56,80), margin=T)

lakes <- map_data("lakes", wrap=c(-25,335), col="white", border="gray", ylim=c(-55,65), margin=T)

#shifting the longlat of the dataframe to match the pacific centered map
gb_world_subset <- gb_world_subset %>%
  mutate(Longitude = if_else(Longitude <= -25, Longitude + 360, Longitude))

#Basemap
basemap <- ggplot(gb_world_subset) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3)  +
  theme(legend.position="none",
        panel.grid.major = element_blank(), #all of these lines are just removing default things like grid lines, axises etc
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())   +
  coord_map(projection = "vandergrinten", ylim=c(-56,67))

basemap <- ggplot(gb_world_subset) +
  geom_polygon(data=world, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="gray87", size = 0.5) +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),
               colour="gray87",
               fill="white", size = 0.3)  +
  theme(
    panel.grid.major = element_blank(), #all of these lines are just removing default things like grid lines, axises etc
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.line = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())   +
  coord_map(projection = "vandergrinten", ylim=c(-56,67))

midpointNw <- (max(gb_world_subset$Nominal_words_complexity) - min(gb_world_subset$Nominal_words_complexity))/2 + min(gb_world_subset$Nominal_words_complexity)

m1 <- basemap + geom_point(
  aes(x=Longitude, y=Latitude, color=Nominal_words_complexity),
  stat = "identity", 
)

m1 <- m1 + ggtitle("Nominal words") + scale_colour_gradient2(low="#2166AC", mid="#FFFFBF", high="#B2182B", midpoint= midpointNw, "Score") + theme(text = element_text(size = 10), legend.title=element_blank()) + theme(plot.margin = unit(c(0,0,0,0), "lines")) #scale_color_viridis(option="A")
m1

midpointV <- (max(gb_world_subset$Verbal_complexity)-min(gb_world_subset$Verbal_complexity))/2 + min(gb_world_subset$Verbal_complexity)

m2 <- basemap + geom_point(
  aes(x=Longitude, y=Latitude, color=Verbal_complexity),
  stat = "identity", 
)
m2 <- m2 + ggtitle("Verbs") + scale_colour_gradient2(low="#2166AC", mid="#FFFFBF", high="#B2182B", midpoint= midpointV, "Score") + theme(text = element_text(size = 10), legend.title=element_blank()) + theme(plot.margin = unit(c(0,0,0,0), "lines"))
m2

ggsave(ggarrange(m1, m2, ncol=1, nrow=2), width=5, height=4.5, filename = here("output", "two_maps_3_col.png"), dpi=300, dev="png")

