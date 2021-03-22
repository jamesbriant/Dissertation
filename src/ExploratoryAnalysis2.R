####################
# Required Packages
####################

# install.packages("lubridate")
# install.packages("ggplot2")
# install.packages("data.table")
# install.packages("ggrepel")
# install.packages("dplyr")
# install.packages("data.table")
# install.packages("tidyverse")
# 
# #install.packages("digest")
# #install.packages("glue")
# 
# if(!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("dkahle/ggmap")

#####################
# Load Libraries
#####################

library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(tidyverse)
library(ggmap)

#####################
# Code
#####################

# read in the raw data
#raw.data <- read.csv("data/20210320135655.csv") # 1st walk
#raw.data <- read.csv("data/20210322154634.csv") # circle walk
raw.data <- read.csv("data/20210322162522.csv") # straight walk
df <- data.frame(x = raw.data$lon,
                 y = raw.data$lat)

# get google cloud key
source("functions/GoogleKeyInfo.R")
MY_KEY <- GetKey()
ggmap::register_google(key = MY_KEY)

ggmap(get_googlemap(center = c(mean(raw.data$lon), 
                               mean(raw.data$lat)),
                    zoom = 16, 
                    size = c(640, 320),
                    scale = 2,
                    maptype ='terrain',
                    color = 'color',
                    markers = df,
                    path = df))



summary(raw.data)



