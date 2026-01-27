library(ggplot2)

BEN<-read.csv("Code/Prevalence/Bovine BCT/Results/Benin/BeninCovariates.csv")
BFA<-read.csv("Code/Prevalence/Bovine BCT/Results/Burkina Faso/Burkina FasoCovariates.csv")
ETH<-read.csv("Code/Prevalence/Bovine BCT/Results/Ethiopia/EthiopiaCovariates.csv")
KEN<-read.csv("Code/Prevalence/Bovine BCT/Results/Kenya/KenyaCovariates.csv")
MWI<-read.csv("Code/Prevalence/Bovine BCT/Results/Malawi/MalawiCovariates.csv")
MOZ<-read.csv("Code/Prevalence/Bovine BCT/Results/Mozambique/MozambiqueCovariates.csv")
NGA<-read.csv("Code/Prevalence/Bovine BCT/Results/Nigeria/NigeriaCovariates.csv")
TGO<-read.csv("Code/Prevalence/Bovine BCT/Results/Togo/TogoCovariates.csv")
UGA<-read.csv("Code/Prevalence/Bovine BCT/Results/Uganda/UgandaCovariates.csv")

covariates<-c("elevation","precipitation","human_fp","pop_den","tree","grassland",
"shrub","cropland","built","bare","water","wetland","mangrove","moss","cattle","tavg","tmin","tmax","tsetse")


top_3_cov<-c(BEN[1:3,1],BFA[1:3,1],ETH[1:3,1],KEN[1:3,1],MWI[1:3,1],MOZ[1:3,1],NGA[1:3,1],TGO[1:3,1],UGA[1:3,1])
count_top_3<-c()

all_cov<-c(BEN[,1],BFA[,1],ETH[,1],KEN[,1],MWI[,1],MOZ[,1],NGA[,1],TGO[,1],UGA[,1])
count_all<-c()

for (i in 1:length(covariates)){

    count_top_3[i]<-length(grep(covariates[i],top_3_cov))
    count_all[i]<-length(grep(covariates[i],all_cov))

}

BEN_df <- as.data.frame(BEN, stringsAsFactors = FALSE)
colnames(BEN_df) <- c("covariate", "value")
BEN_df$value <- as.numeric(BEN_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(BEN_df$covariate, full_df$covariate)] <- BEN_df$value
BEN_final <- as.matrix(full_df)
BEN_cov <-as.numeric(BEN_final[,2])/sum(as.numeric(BEN_final[,2]))

BFA_df <- as.data.frame(BFA, stringsAsFactors = FALSE)
colnames(BFA_df) <- c("covariate", "value")
BFA_df$value <- as.numeric(BFA_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(BFA_df$covariate, full_df$covariate)] <- BFA_df$value
BFA_final <- as.matrix(full_df)
BFA_cov <-as.numeric(BFA_final[,2])/sum(as.numeric(BFA_final[,2]))

ETH_df <- as.data.frame(ETH, stringsAsFactors = FALSE)
colnames(ETH_df) <- c("covariate", "value")
ETH_df$value <- as.numeric(ETH_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(ETH_df$covariate, full_df$covariate)] <- ETH_df$value
ETH_final <- as.matrix(full_df)
ETH_cov <-as.numeric(ETH_final[,2])/sum(as.numeric(ETH_final[,2]))

KEN_df <- as.data.frame(KEN, stringsAsFactors = FALSE)
colnames(KEN_df) <- c("covariate", "value")
KEN_df$value <- as.numeric(KEN_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(KEN_df$covariate, full_df$covariate)] <- KEN_df$value
KEN_final <- as.matrix(full_df)
KEN_cov <-as.numeric(KEN_final[,2])/sum(as.numeric(KEN_final[,2]))

MWI_df <- as.data.frame(MWI, stringsAsFactors = FALSE)
colnames(MWI_df) <- c("covariate", "value")
MWI_df$value <- as.numeric(MWI_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(MWI_df$covariate, full_df$covariate)] <- MWI_df$value
MWI_final <- as.matrix(full_df)
MWI_cov <-as.numeric(MWI_final[,2])/sum(as.numeric(MWI_final[,2]))

MOZ_df <- as.data.frame(MOZ, stringsAsFactors = FALSE)
colnames(MOZ_df) <- c("covariate", "value")
MOZ_df$value <- as.numeric(MOZ_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(MOZ_df$covariate, full_df$covariate)] <- MOZ_df$value
MOZ_final <- as.matrix(full_df)
MOZ_cov <-as.numeric(MOZ_final[,2])/sum(as.numeric(MOZ_final[,2]))

NGA_df <- as.data.frame(NGA, stringsAsFactors = FALSE)
colnames(NGA_df) <- c("covariate", "value")
NGA_df$value <- as.numeric(NGA_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(NGA_df$covariate, full_df$covariate)] <- NGA_df$value
NGA_final <- as.matrix(full_df)
NGA_cov <-as.numeric(NGA_final[,2])/sum(as.numeric(NGA_final[,2]))

TGO_df <- as.data.frame(TGO, stringsAsFactors = FALSE)
colnames(TGO_df) <- c("covariate", "value")
TGO_df$value <- as.numeric(TGO_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(TGO_df$covariate, full_df$covariate)] <- TGO_df$value
TGO_final <- as.matrix(full_df)
TGO_cov <-as.numeric(TGO_final[,2])/sum(as.numeric(TGO_final[,2]))

UGA_df <- as.data.frame(UGA, stringsAsFactors = FALSE)
colnames(UGA_df) <- c("covariate", "value")
UGA_df$value <- as.numeric(UGA_df$value)
full_df <- data.frame(covariate = covariates, stringsAsFactors = FALSE)
full_df$value <- 0
full_df$value[match(UGA_df$covariate, full_df$covariate)] <- UGA_df$value
UGA_final <- as.matrix(full_df)
UGA_cov <-as.numeric(UGA_final[,2])/sum(as.numeric(UGA_final[,2]))

cov_all<-BEN_cov+BFA_cov+ETH_cov+KEN_cov+MWI_cov+MOZ_cov+NGA_cov+TGO_cov+UGA_cov
cov_all<-cov_all/sum(cov_all)



#### TOP 3

df <- data.frame(covariate = covariates, count = count_top_3, stringsAsFactors = FALSE)
df <- df[order(-df$count), ]
df$covariate <- factor(df$covariate, levels = df$covariate)

ggplot(df, aes(x = covariate, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "",
       x = "Covariate",
       y = "Count") +
  theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Code/Prevalence/Bovine BCT/Results/Top_3_covariates.png",width = 12, height = 8)


#### ALL COVS

df <- data.frame(covariate = covariates, count = count_all, stringsAsFactors = FALSE)
df <- df[order(-df$count), ]
df$covariate <- factor(df$covariate, levels = df$covariate)

ggplot(df, aes(x = covariate, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "",
       x = "Covariate",
       y = "Count") +
  theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Code/Prevalence/Bovine BCT/Results/All_covariates.png",width = 12, height = 8)


##### COVARIATE IMPORTANCE

df <- data.frame(covariate = covariates, count = cov_all, stringsAsFactors = FALSE)
df <- df[order(-df$count), ]
df$covariate <- factor(df$covariate, levels = df$covariate)

ggplot(df, aes(x = covariate, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "",
       x = "Covariate",
       y = "Importance") +
  theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Code/Prevalence/Bovine BCT/Results/Covariate_importance.png",width = 12, height = 8)
