library(readxl)

sharma <- read_excel("../data/Additional_data/potato_breeds_in_sharma_2018.xlsx", skip = 2)
load("Samples")
sum(toupper(sharma$Name) %in% toupper(Samples$VARIETY))
