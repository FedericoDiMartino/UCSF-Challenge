## Load packages
library(plyr)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

#### Load in data
## Set working directory
setwd("./data")

## Create file names
filenames <- list.files(pattern="*.csv")

## Read in files

for (i in 1:length(filenames)) assign(filenames[i], read.csv(filenames[i]))
## length(filenames) = 129 if run first time, 130 after unsupp.csv is created

## Create list without ViralLoads.csv
community_period_object_list <- mget(filenames[1:128])

## Create namelist without ViralLoads.csv
community_period_name_list <- filenames[1:128]


## Function to generate lists of community names and prop_unsupp_t for each period
names_periods <- function(z){
    
    
    ## Get name of community and period    
    
    name_time <-   z   
    
    name_pattern <- "[^_]*" ## Everything up to "_"
    "community name" <- str_extract(name_time, name_pattern)
    
    period_pattern <- "\\d" ## Only numbers
    period_id <- str_extract(name_time, period_pattern)
    
    
    if(period_id == "0"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_0")    
        return(output)
        
        
        
        
    }else if (period_id == "1"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_1")    
        return(output)
        
        
        
    }else if (period_id == "2"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_2")    
        return(output)
        
        
        
    }
    else if (period_id == "3"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_3")    
        return(output)
        
        
        
    }else
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = period_id) 
    return(output)
    
}
##

## Run function 
left_columns <- lapply(community_period_name_list, names_periods)

## Collate
left_columns <- ldply(left_columns, data.frame)

## Function to calculate proportion of patients who have an unsuppressed viral load
## for each community_period combination
calculate_proportions <- function(x, ViralLoads) {   
    
    ## ViralLoads = ViralLoads.csv always
    
    ## join ViralLoads onto community_period table 
    y <- left_join(x, ViralLoads, by = "braceletid")
    
    ## Convert dates into POSIXct date-time objects
    y[,"chcdate"] <- ymd(y[,"chcdate"])
    y[,"trdate"] <- ymd(y[,"trdate"])    
    y[,"date"] <- ymd(y[,"date"])
    
    ## Calculating chcstart and trkend
    chcstart_t <- min(y[,"chcdate"] , na.rm = TRUE)
    trkend_t <- max( y[,"trdate"] , na.rm = TRUE)
    
    ## Some chcstart values are an obviously incorrect 1899-01-01
    ## decided to interpret 1899 as what should have been an NA
    ## filter them out unless they have a trdate
    y <- filter(y, ! (chcdate == ymd("1899-01-01") & is.na(trdate)) )
    
    
    ## If they have a trdate, set 1899 chcdate to NA,
    y[,"chcdate"][which(y[,"chcdate"] == ymd("1899-01-01"))] <- NA
    
    
    ## Since the thing we want measure is defined by the period it occurs in,
    ## filter out all observations where:
    ## chcdate AND trdate == NA 
    ## VL != NA AND date == NA (i.e. VL test has been done but we don't know when)
    
    y <- filter(y, !( is.na(chcdate) & is.na(trdate) ) )
    
    y <- filter(y, !( !is.na(VL) & is.na(date) ) )
    
    ## Re-calculate chcstart
    chcstart_t <- min(y[,"chcdate"] , na.rm = TRUE)
    
    
    ## In case of multiple measurements, only viral load measurement performed 
    ## closest to chcstart (i.e. the earliest) is desired, remove others
    ## y <- group_by(y, date)
    
    y <- arrange(y, date)
    y <- distinct(y, searchid , .keep_all = TRUE)
    
    
    
    ## Calculating unsupp VL
    
    ## Filter out cases where VL test has been performed (date != NA)
    ## and VL = NA (i.e test success unknown, we can't know whether VL is unsuppresed or not)
    y <- filter(y, !( is.na(VL) & !is.na(date) ) )
    
    ## Add unsupp_test column, determining value of unsupp_t
    y <- mutate(y, unsupp_t = if_else(HIV == 1 & VL > 500 & !is.na(VL), 1, 0) )
    
    ##
    
    prop_unsupp_calculated <- mean(y[,"unsupp_t"])
    
    
    
    return(prop_unsupp_calculated)
    
}

##

## Run function for every combination
right_column <- lapply(community_period_object_list, calculate_proportions, ViralLoads.csv )

## Collate
right_column <- ldply(right_column, data.frame)
## Rename column, then remove useless column
names(right_column) <- c(".id" ,  "prop_unsupp_calculated")
right_column <- select(right_column, "prop_unsupp_calculated")

## Join columns
joined <- bind_cols(left_columns, right_column)

## reshape data into required format, rename column
unsupp <- spread(joined, prop_unsupp_t, prop_unsupp_calculated )
names(unsupp) <- c("community name", "prop_unsupp_0",  "prop_unsupp_1",
                   "prop_unsupp_2",  "prop_unsupp_3")

## write to file
write.csv(unsupp, "unsupp.csv")