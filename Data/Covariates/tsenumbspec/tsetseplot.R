#plot tsetse fly data
library(ggplot2)
library(afrilearndata)

africa <- africontinent
#extract the border as longitude and latitude values
border <- st_coordinates(st_geometry(africa))

#remove Madagascar
border <- border[!(border[,1] > 40 & border[,1] < 60 & border[,2] > -30 & border[,2] < -10),]
border <- border[,1:2]
#change to data frame
border_df <- data.frame(long = border[,1], lat = border[,2])

# Load the data
tsetse<-raster("Data/Covariates/tsenumbspec")
tsetse[tsetse > 1] <- 1
# Convert raster to data frame for ggplot
tsetse_df <- as.data.frame(rasterToPoints(tsetse), xy = TRUE)


#plot tsetse data with border
ggplot() +
  geom_raster(data = tsetse_df, aes(x = x, y = y, fill = tsenumbspec)) +
  scale_fill_gradient(low = "white", high = "blue", na.value = "transparent") +
  geom_polygon(data = border_df, aes(x = long, y = lat),fill=NA, color = "black") +
  coord_fixed() +
  labs(title = "Tsetse Fly Distribution in Africa",
       x = "Longitude",
       y = "Latitude",
       fill = "Tsetse Fly Presence") +
  theme_minimal()
# Save the plot
ggsave("Data/Covariates/tsenumbspec/tsetse_plot.pdf", width = 10, height = 6)
