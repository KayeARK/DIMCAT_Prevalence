library(readxl)
library(geoR)
library(ggplot2)
library(INLA)
library(sp)
library(afrilearndata)

#####

data <- read_excel("Data/ContAtlas_v2/250318_AT_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data <- data[!is.na(data$Trypanosome_detection) | !is.na(data$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data$Trypanosome_detection <- as.factor(data$Trypanosome_detection)
data$Trypanosome_detection <- factor(data$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data$Trypanosome_detection <- as.numeric(data$Trypanosome_detection) - 1
presence <- data$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_AT.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data$Longitude, data$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Tsi_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Tsi.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Tb_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Tb.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Tc_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Tc.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Te_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Te.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Tg_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Tg.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Tsu_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Tsu.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

#####

data_s <- read_excel("Data/ContAtlas_v2/250318_Tv_Presence.xls")
#remove rows where Trypanosome_detection or Number_of_infections are empty
data_s <- data_s[!is.na(data_s$Trypanosome_detection) | !is.na(data_s$Number_of_infections),]

#colour the data points by Trypanosome_detection. No is black, yes is red
data_s$Trypanosome_detection <- as.factor(data_s$Trypanosome_detection)
data_s$Trypanosome_detection <- factor(data_s$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data_s$Trypanosome_detection <- as.numeric(data_s$Trypanosome_detection) - 1

#data_c <- data
#remove rows from data_c that are in data_s
#data_c <- data_c[!(data_c$Longitude %in% data_s$Longitude & data_c$Latitude %in% data_s$Latitude),]
#set Trypanosome_detection to 0 for these rows
#data_c$Trypanosome_detection <- 0
#combine data_c and data_s
#data_s <- rbind(data_c, data_s)
presence <- data_s$Trypanosome_detection

pdf("Code/AbsencePresence/Results/Images/DataPlot_Tv.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(data_s$Longitude, data_s$Latitude, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

