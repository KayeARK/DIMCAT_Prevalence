library(readxl)
library(writexl)

# Load the data
data_geo <- read_excel("Data/240708_FAO_Continental_Atlas_AT_1990_2020.xlsx", sheet="Geo_data")
data_epi <- read_excel("Data/240708_FAO_Continental_Atlas_AT_1990_2020.xlsx", sheet="Epi_data")

longitudes <- data_geo$LONG
latitudes <- data_geo$LAT
area <- data_geo$AREA_TYPE_AT

epi_ID <- data_epi$LOCATION_ID
geo_ID <- data_geo$LOCATION_ID


for (i in 1:length(epi_ID)){
  for (j in 1:length(geo_ID)){
    if (epi_ID[i] == geo_ID[j]){
      data_epi$LONG[i] <- data_geo$LONG[j]
      data_epi$LAT[i] <- data_geo$LAT[j]
      data_epi$AREA_TYPE_AT[i] <- data_geo$AREA_TYPE_AT[j]
    }
  }
}

#remove rows where the TPR is empty
data_epi <- data_epi[!is.na(data_epi$TPR),]

#remove rows where the T is empty
data_epi <- data_epi[!is.na(data_epi$T),]

#remove rows where T is greater than the sample size
data_epi <- data_epi[data_epi$T <= data_epi$SAMPLE_SIZE,]

#remove rows where the sample size is empty
data_epi <- data_epi[!is.na(data_epi$SAMPLE_SIZE),]

#remove rows where the latitude or longitude is empty
data_epi <- data_epi[!is.na(data_epi$LAT),]
data_epi <- data_epi[!is.na(data_epi$LONG),]

#sum the number of positive cases and the sample size for each location keeping latitude and longitude
#data_epi <- aggregate(cbind(T, SAMPLE_SIZE) ~ LAT + LONG, data_epi, sum)

#save the data as an excel file

write_xlsx(data_epi, "Data/RefinedData.xlsx")