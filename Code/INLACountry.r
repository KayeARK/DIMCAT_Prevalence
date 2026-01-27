library(readxl)
library(ggplot2)
library(INLA)
library(sp)
library(sf)
library(afrilearndata)
library(concaveman)

data <- read_excel("Data/ContAtlas_v2/250318_AT_Presence.xls")

###### TIDY DATA, REMOVE ROWS WHERE PREVALENCE CANNOT BE INFERRED

data <- data[!is.na(data$Number_of_infections) & !is.na(data$Number_of_animal_tested) & !is.na(data$Diagnostic_CATEGORY),]

#remove any rows where Number_of_infections/Number_of_animal_tested is greater than 1
data <- data[data$Number_of_infections/data$Number_of_animal_tested <= 1,]
######

positive=ceiling(data$Number_of_infections)
positive[positive>0]=1
latitudes=data$Latitude
longitudes=data$Longitude
sample_size=data$Number_of_animal_tested

countries_to_infer=c("Nigeria")
countries_to_infer <- sort(countries_to_infer)

#if Africa is in this list then perform the calculations for the whole continent
if("Africa" %in% countries_to_infer){
  border <- st_coordinates(st_geometry(africontinent))

#remove Madagascar
  border <- border[!(border[,1] > 40 & border[,1] < 60 & border[,2] > -30 & border[,2] < -10),]
  border <- border[,1:2]
  #remove Madagascar from africountries
  africountries <- africountries[!(africountries$name=="Madagascar"),]
  country <- africountries
} else{
v=africountries$name %in% countries_to_infer
country <- africountries[which(v==TRUE),]

}

#remove all countries except the ones in countries_to_infer
africountries_toplot <- africountries[which(africountries$name %in% countries_to_infer),]

country <- st_coordinates(st_geometry(country))

poly.df <- data.frame("long"=country[,1],"lat"=country[,2])
poly.sf <- st_as_sf(poly.df, coords = c("long","lat"))
poly <- concaveman(poly.sf, length_threshold = 0,concavity=1.1)
point.df <- data.frame("long"=longitudes,"lat"=latitudes)
point.sf <- st_as_sf(point.df, coords = c("long","lat"))
long_lat_concave <- st_coordinates(st_geometry(poly))
border <- long_lat_concave

inside <- st_intersects(point.sf,poly)
inside <- inside[]==1

latitudes <- latitudes[inside]
longitudes <- longitudes[inside]
sample_size <- sample_size[inside]
positive <- positive[inside]

#remove NA values
ind <- !is.na(latitudes)
latitudes <- latitudes[ind]
longitudes <- longitudes[ind]
sample_size <- sample_size[ind]
positive <- positive[ind]

prevalence <- positive/sample_size

# put latitude and longitude next to each other
coo <- cbind(longitudes, latitudes)

mesh <- inla.mesh.2d(
  loc = coo, offset = c(50, 100),
  cutoff = 1, max.edge = c(30, 60)
)

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

indexs <- inla.spde.make.index("s", spde$n.spde)

A <- inla.spde.make.A(mesh = mesh, loc = coo)

bb <- bbox(border)
x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 200)
y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 200)
coop <- as.matrix(expand.grid(x, y))

ind <- point.in.polygon(
  coop[, 1], coop[, 2],
  border[, 1], border[, 2]
)
coop <- coop[which(ind == 1), ]


Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = positive, numtrials=sample_size),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA, numtrials=NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ 0 + b0 + f(s, model = spde)

res <- inla(formula,
  data = inla.stack.data(stk.full),
  family = "binomial", Ntrials = numtrials,
  control.compute=list(return.marginals.predictor=TRUE),
  control.predictor = list(link=1,
    compute = TRUE,
    A = inla.stack.A(stk.full)
  )
)

index <- inla.stack.index(stk.full, tag = "pred")$data

pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

# plot

dpm <- rbind(
  data.frame(
    Latitude = coop[, 1], Longitude = coop[, 2],
    value = pred_mean, variable = "Mean"
  ),
  data.frame(
    Latitude = coop[, 1], Longitude = coop[, 2],
    value = pred_ll, variable = "2.5th percentile"
  ),
  data.frame(
    Latitude = coop[, 1], Longitude = coop[, 2],
    value = pred_ul, variable = "97.5th percentile"
  )
)
dpm$variable <- as.factor(dpm$variable)


ggplot(dpm) + geom_tile(aes(Latitude, Longitude, fill = value)) + geom_polygon(data =  border, aes(x = X, y = Y), color = "black",fill=NA)+
  facet_grid(~factor(variable,levels=c("2.5th percentile","Mean","97.5th percentile"))) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Prevalence",
    low = "blue", high = "orange",
    limits = c(0, 1)
  ) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) 
ggsave("Images/Predictions_SpecificCountry.png",width = 12, height = 8)


#for each element in coo find the index of the closest point in coop using l2 norm
closest_point_index <- numeric(nrow(coo))
for(i in 1:nrow(coo)){
  closest_point_index[i] <- which.min(rowSums((coop-t(matrix(coo[i,],2,nrow(coop))))^2))
}

marginals <-  array(NA,c(43,2,nrow(coo)))
for(i in 1:nrow(coo)){
  marg <- res$marginals.fitted.values[index][[closest_point_index[i]]]
  marginals[,,i] <- marg
}

#generate 4 random numbers between 1 and length(coo)/2
random_numbers <- sample(1:(length(coo)/2),4)

#plot a line graph of the first marginal
for(i in 1:(length(coo)/2)){
marginal <- marginals[,,i]
v[i]<-(prevalence[i]-marginal[which.max(marginal[,2])])
}
png("Images/Posteriors_Africa/Hist.png")
hist(v,main="",xlab="Difference between true prevalence and mode of the prevalence posterior",ylab="Frequency",breaks=30,xlim=c(-1,1))
dev.off()
#plot a line graph of the first marginal
k<-1
for(i in random_numbers){
marginal <- marginals[,,i]
png(paste("Images/Posteriors_Africa/Marginal_",k,".png",sep=""))
plot(marginal, type="l", xlab="Prevalence", ylab="Density", main="",xlim=c(min(min(marginal[,1],prevalence[i])),max(max(marginal[,1],prevalence[i]))))
abline(v=prevalence[i], col="red")
dev.off()
k<-k+1
print(k)
}
