## Load packages
library(plyr)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

## set seed
set.seed(12345)

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


## Function to generate lists of community names, prop_unsupp_t and
## HIV_pos_given_ART_t for each period
names_periods <- function(z){
    
    
    ## Get name of community and period    
    
    name_time <-   z   
    is
    name_pattern <- "[^_]*" ## Everything up to "_"
    "community name" <- str_extract(name_time, name_pattern)
    
    period_pattern <- "\\d" ## Only numbers
    period_id <- str_extract(name_time, period_pattern)
    
    
    if(period_id == "0"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_0",
                             "HIV_pos_given_ART_t" = "HIV_pos_given_ART_0")    
        return(output)
        
        
        
        
    }else if (period_id == "1"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_1",
                             "HIV_pos_given_ART_t" = "HIV_pos_given_ART_1")    
        return(output)
        
        
        
    }else if (period_id == "2"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_2",
                             "HIV_pos_given_ART_t" = "HIV_pos_given_ART_2")    
        return(output)
        
        
        
    }
    else if (period_id == "3"){
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = "prop_unsupp_3",
                             "HIV_pos_given_ART_t" = "HIV_pos_given_ART_3")    
        return(output)
        
        
        
    }else
        
        output <- data_frame("community name" = `community name`,
                             "prop_unsupp_t" = period_id,
                             "HIV_pos_given_ART_t" = "HIV_pos_given_ART_t") 
    return(output)
    
}
##

## Run function 
left_columns <- lapply(community_period_name_list, names_periods)

## Collate
left_columns <- ldply(left_columns, data.frame)

## Function to calculate proportion of patients who have an unsuppressed viral load
## and ART given/HIV positive ratio
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
    
    ## Calculating ART given/HIV positive ratio
    y_only_HIV_POS <- filter(y, HIV == 1 & !is.na(ART))
    HIV_pos_given_ART_calculated <- (sum(y_only_HIV_POS[,"ART"], na.rm = TRUE) / 
                                     sum(y_only_HIV_POS[,"HIV"], na.rm = TRUE)
                                    )
    
    ## return calculated values
    output <-  data.frame(prop_unsupp_calculated, HIV_pos_given_ART_calculated)
    return(output)
}

##

## Run function for every combination
right_columns <- lapply(community_period_object_list, calculate_proportions, ViralLoads.csv )

## Collate
right_columns <- ldply(right_columns, data.frame)
## Remove useless column
right_columns <- select(right_columns, "prop_unsupp_calculated",
                        "HIV_pos_given_ART_calculated")

## Join columns
joined <- bind_cols(left_columns, right_columns)



#### reshape data into required format, rename column
joined_unsupp <- select(joined, -HIV_pos_given_ART_t, -HIV_pos_given_ART_calculated)
joined_HIV <- select(joined, -prop_unsupp_t, -prop_unsupp_calculated)


unsupp <- spread(joined_unsupp, prop_unsupp_t, prop_unsupp_calculated )
HIV <- spread(joined_HIV, HIV_pos_given_ART_t, HIV_pos_given_ART_calculated)

names(unsupp)[1] <- "community name"
names(HIV)[1] <- "community name"

rejoined <- left_join(unsupp, HIV, by = "community name")

## Building linear model to answer:
## Suppose we changed our data simulation so that all patients who are HIV positive
## at time 1, 2, or 3 are immediately treated with antiretroviral therapy (ART).
## The data generating process is otherwise unchanged (including treatment at time 0).
## In the resulting data, what would be the total population proportion of patients with an unsuppressed viral load at time 3? Provide a
## single estimate and a 95% confidence interval. Very briefly describe your methodology.

## Variables selected for model are three we will have altered in the prediction
## and two (prop_unsupp_0 + HIV_pos_given_ART_0) which will be calculated for whole population
## in existing data, which will be used in prediction

predictmodel <- lm(prop_unsupp_3 ~ prop_unsupp_0 + HIV_pos_given_ART_0 
                 + HIV_pos_given_ART_1 + HIV_pos_given_ART_2
                 + HIV_pos_given_ART_3 , data = rejoined)



####
## Calculate prop_unsupp_0, HIV_pos_given_ART_0 for combination of every community
## Select and combine data for all communities at t = 0
t_0_pattern <- "[0]"
communities_t_0 <- grep(t_0_pattern, community_period_name_list, value = TRUE)
community_t_0_object_list <-mget(communities_t_0)
communities_t_0 <- ldply(community_t_0_object_list, data.frame)
communities_t_0 <- left_join(communities_t_0 ,ViralLoads.csv, by = "braceletid")

### Now Calculate prop_unsupp_0, HIV_pos_given_ART_0
## Taken and modified from body of earlier function

## Convert dates into POSIXct date-time objects
communities_t_0[,"chcdate"] <- ymd(communities_t_0[,"chcdate"])
communities_t_0[,"trdate"] <- ymd(communities_t_0[,"trdate"])    
communities_t_0[,"date"] <- ymd(communities_t_0[,"date"])



## Some chcstart values are an obviously incorrect 1899-01-01
## decided to interpret 1899 as what should have been an NA
## filter them out unless they have a trdate
communities_t_0 <- filter(communities_t_0, ! (chcdate == ymd("1899-01-01") & is.na(trdate)) )


## If they have a trdate, set 1899 chcdate to NA,
communities_t_0[,"chcdate"][which(communities_t_0[,"chcdate"] == ymd("1899-01-01"))] <- NA


## Since the thing we want measure is defined by the period it occurs in,
## filter out all observations where:
## chcdate AND trdate == NA 
## VL != NA AND date == NA (i.e. VL test has been done but we don't know when)

communities_t_0 <- filter(communities_t_0, !( is.na(chcdate) & is.na(trdate) ) )

communities_t_0 <- filter(communities_t_0, !( !is.na(VL) & is.na(date) ) )

## Calculating chcstart and trkend
chcstart_0 <- min(communities_t_0[,"chcdate"] , na.rm = TRUE)
trkend_0 <- max( communities_t_0[,"trdate"] , na.rm = TRUE)


## In case of multiple measurements, only viral load measurement performed 
## closest to chcstart (i.e. the earliest) is desired, remove others


communities_t_0 <- arrange(communities_t_0, date)
communities_t_0 <- distinct(communities_t_0, searchid , .keep_all = TRUE)



## Calculating unsupp VL

## Filter out cases where VL test has been performed (date != NA)
## and VL = NA (i.e test success unknown, we can't know whether VL is unsuppresed or not)
communities_t_0 <- filter(communities_t_0, !( is.na(VL) & !is.na(date) ) )

## Add unsupp_test column, determining value of unsupp_t
communities_t_0 <- mutate(communities_t_0, unsupp_t = if_else(HIV == 1 & VL > 500 & !is.na(VL), 1, 0) )

##

communities_t_0_prop_unsupp <- mean(communities_t_0[,"unsupp_t"])

## Calculating ART given/HIV positive ratio
communities_t_0_only_HIV_POS <- filter(communities_t_0, HIV == 1 & !is.na(ART))
communities_t_0_HIV_pos_given_ART <- (sum(communities_t_0_only_HIV_POS[,"ART"], na.rm = TRUE) / 
                                     sum(communities_t_0_only_HIV_POS[,"HIV"], na.rm = TRUE)
                                )

##############################


predictdata <- data.frame(prop_unsupp_0 = communities_t_0_prop_unsupp, 
                              HIV_pos_given_ART_0 = communities_t_0_HIV_pos_given_ART,
                              HIV_pos_given_ART_1 = 1, HIV_pos_given_ART_2 = 1,
                              HIV_pos_given_ART_3 = 1)

prediction <- predict(predictmodel, newdata = predictdata, interval = "confidence", level = 0.95)
prediction
prediction[3] - prediction[2]


### Estimate is 0.05230103 ± 0.07384524