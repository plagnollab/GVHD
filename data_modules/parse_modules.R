library(gdata)


data <- read.xls("data_modules/Kasela_et_al_Plos.xls", skip = 1, header = TRUE)

library(dplyr)
data <- dplyr::select(data, SNPName, HGNCName) %>% distinct()
