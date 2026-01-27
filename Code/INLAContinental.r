library(readxl)
#library(geoR)
#library(ggplot2)
library(INLA)
library(sp)
library(afrilearndata)

africa <- africontinent
#extract the border as longitude and latitude values
border <- st_coordinates(st_geometry(africa))

#remove Madagascar
border <- border[!(border[,1] > 40 & border[,1] < 60 & border[,2] > -30 & border[,2] < -10),]
border <- border[,1:2]

data <- read_excel("Data/ContAtlas_v2/250318_AT_Presence.xls")

###### TIDY DATA, REMOVE ROWS WHERE PREVALENCE CANNOT BE INFERRED

data <- data[!is.na(data$Number_of_infections) & !is.na(data$Number_of_animal_tested) & !is.na(data$Diagnostic_CATEGORY),]

#remove any rows where Number_of_infections/Number_of_animal_tested is greater than 1
data <- data[data$Number_of_infections/data$Number_of_animal_tested <= 1,]
######

positive=ceiling(data$Number_of_infections)
#positive[positive>0]=1
latitudes=data$Latitude
longitudes=data$Longitude
sample_size=data$Number_of_animal_tested

# put latitude and longitude next to each other
coo <- cbind(longitudes, latitudes)

mesh <- inla.mesh.2d(
  loc = coo, offset = c(50, 100),
  cutoff = 1, max.edge = c(30, 60)
)

pdf("Images/AfricaMesh.pdf")
plot(mesh)
dev.off()
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
    east = coop[, 1], north = coop[, 2],
    value = pred_mean, variable = "Mean"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ll, variable = "2.5th percentile"
  ),
  data.frame(
    east = coop[, 1], north = coop[, 2],
    value = pred_ul, variable = "97.5th percentile"
  )
)
dpm$variable <- as.factor(dpm$variable)


ggplot(dpm) + geom_tile(aes(east, north, fill = value)) + geom_polygon(data =  border, aes(x = X, y = Y), color = "black",fill=NA)+
facet_grid(~factor(variable,levels=c("2.5th percentile","Mean","97.5th percentile"))) +
  facet_wrap(~variable, nrow = 1) +
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
ggsave("Images/Predictions.png",width = 12, height = 8)
