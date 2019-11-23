
# Functions for returning means pluss minus SD or SE
sem <- function(x) sqrt(var(x, na.rm=T)/nanaLength(x))
nanaLength = function(x) sum(!is.na(x))

#Function for getting the mean temperature of a time period
#Dates must be supplied as a lubridate time object
calc_meanTemp = function(df_temp, startDate, endDate, column) {
  if(is.na(startDate) | is.na(endDate) | is.na(column)) return(NA)
  return(
    unlist(df_temp[column])[df_temp$date < endDate & df_temp$date > startDate] %>% mean(na.rm=T)
  )
}  

extend_tempData = function(df_temp, date, columns){
  date=dmy(date)
  # get the last date of the temperature dataset
  date_last <- df_temp$date[nrow(df_temp)]
  message("Last recorded date is: ", date_last)
  
  # get the days between that date and the extension date and make this into a dataframe (rows being dates)
  betweenDates <- datesBetween(date, date_last)
  df_new <- data.frame(date=betweenDates)

  # for each column we want to extend:
  for (i in columns){
    # apply over the list of days, return the temperature of that date one year ago
    df_new[[i]] = apply(df_new,MARGIN=1,FUN=function(j){
      day = ymd(j[["date"]]) 
      day_prev = df_temp %>% filter(date == (day-years(1)))
      temp_prev = day_prev[[i]]
    })
  }

    
  message("Data extended to ",date)
  # rbind the new dataset to the old one
  return(bind_rows(df_temp,df_new))
  
}


# calculates thermal growth coefficient (from Jobling: The thermal growth coefficient (TGC) model of fish growth: a cautionary note)
calc_TGC = function(W0,W1,temp,time) {
  time = as.numeric(time)
  TGC = ((W1^(1/3)-W0^(1/3))/(temp*time)*1000)
}

#' improved TGC calculator. Returns a vector
adv_calc_TGC <- function(df_fish,df_temp,period)
{
  time <- dmy(df_fish[[glue("date.{period}")]]) - dmy(df_fish[[glue("date.{period-1}")]])
  W0   <- df_fish[[glue("weight.{period-1}")]]
  W1   <- df_fish[[glue("weight.{period}")]]
  
  temp <- mapply(calc_meanTemp, 
                startDate= dmy(df_fish[[glue("date.{period-1}")]]), 
                endDate  = dmy(df_fish[[glue("date.{period}")]]), 
                column   = df_fish$temp,
                MoreArgs = list(df_temp=df_temp))
  
  TGC <- calc_TGC(W0, W1, temp, time)
}

# predicts growth based on TGC, same model as above
predict_growth = function(W0, TGC, temp, time){
  time = as.numeric(time)
  return(
    as.numeric((W0^(1/3)+((TGC/1000)*(temp*time)) )^3)
  )
}


# returns the feed size for a list of fish weights
calc_feedSizeDist = function(weights){
  ret = case_when(weights <= 15 ~ 1.2,
                  weights > 15  & weights <= 30  ~ 1.7,
                  weights > 30  & weights <= 70  ~ 2.5,
                  weights > 70  & weights <= 125 ~ 3.5,
                  weights > 125 & weights <= 2000 ~ 5,
                  weights > 500 ~ 100,
                  TRUE ~ NaN)
  return(ret)
}

# uses the same function as above, but returns results as percentage for each feed size
calc_feedSizePercentage = function(weights) {
  weightClasses = calc_feedSizeDist(weights)
  dat = data.frame(weightClasses)
  totalLength = length(dat$weightClasses[!is.na(dat$weightClasses)])
  
  return(
    dat %>% group_by(weightClasses) %>% summarise( perc = length(weightClasses)/totalLength*100 )
  )
}

# calculates the total amount of feed given between two dates
calc_totalFeed = function(dates, feed, date1, date2) {
  return(
    sum( feed[ dates > date1 & dates < date2], na.rm=T )
  )
}

# calculates the amount of feed given a specific tank between two dates
calc_feedPrTank = function(df_feed,tankID, date1, date2) {
  data_feed_tank = df_feed[which(df_feed$tank == tankID),]
  amount = calc_totalFeed(data_feed_tank$date, data_feed_tank$feed, date1, date2)
  return(amount)
}

# alternative, "smoothed" version of funcion above. Uses feed data that gives feed pr every single day (averaged from original data) 
calc_feedPrTank2 = function(tankID, date1, date2) {
  data_feed_tank = data_feedDays[which(data_feedDays$tank == tankID),]
  amount = calc_totalFeed(data_feed_tank$date, data_feed_tank$feed, date1, date2)
  return(amount)
}

#' Function for gaining a list of all dates between two dates
#' sollution to this problem by user "yifyan" at stackoverflow.com
#' https://stackoverflow.com/questions/14450384/create-a-vector-of-all-days-between-two-dates
#' date_a and date_b must be lubridate-date-objects
datesBetween = function(date_a,date_b) { 
  require(lubridate)
  
  n_days <- interval(date_a,date_b)/days(1)
  dates = date_a + days(0:n_days)
  return(dates)
}

# function for converting the feed datasat to a "pr-day" type dataset
feedPrDate = function(df_feed,tank,date) {
  #BE AWARE: this function does not work well for the very last period of feeding
  
  #find the last date where this tank was "fed"
  prevDates = df_feed$date[df_feed$date <= date & df_feed[tank] != 0]
  prevDate  = max(prevDates,na.rm=T)
  #get the next date where this tank is being "fed"
  comingDates = df_feed$date[df_feed$date > date]
  nextDate = min(comingDates,na.rm=T)
  
  #get the amount of feed at the intial date
  amount = as.numeric(df_feed[tank][df_feed$date == prevDate,1])
  
  #get the number of days between last and next date
  numDays = as.numeric(nextDate-prevDate)
  
  #divide amount fed by number of dates
  dailyAmount = amount/numDays
  return(dailyAmount[1])
  
}

#' predict_weights
#' 
#' 
predict_weights <- function(df_fish, df_temp, date, startPeriod=999){
  message("Predicting weights...")
  # obtain all the periods so far as a vector (of numbers)
  periods <- colnames(df_fish) %>%
    str_extract("weight[.][0-9]*") %>%
    sub("weight.","",.) %>%
    numextract() %>%
    na.omit()
  
  # remove periods that are after the specified startPeriod
  periods = periods[periods <= startPeriod+1]
            
  
  # calculate TGC for all used periods
  tgc <- as.data.frame(sapply(periods[2:length(periods)], FUN=function(i){
    adv_calc_TGC(df_fish,df_temp,i)
  }))
  message("TGC calculated from periods: ",periods[2:length(periods)])

  # create one new tgc that is the mean of these TGC values, but weighting the last one five times
  tgc_weighted <- apply(tgc,1,FUN=function(x){
    x=na.omit(x)
    if(length(x)==0) return(NA)
    weighted.mean(x,c(rep(1,length(x)-1),5),)
  }) 
  
  date_cur   <- glue("date.{max(periods)}")
  temp = mapply(
    calc_meanTemp,
    startDate = dmy(df_fish[[date_cur]]),
    endDate   = date,
    column    = df_fish$temp,
    MoreArgs  = list(df_temp=df_temp))
  
  
  W0 <- df_fish[[glue("weight.{max(periods)}")]]
  weight_new <- glue("weight.{max(periods)+1}")
  
  df_fish[[weight_new]] <- predict_growth(
    W0,
    TGC  = tgc_weighted,
    time = date-dmy(df_fish[[date_cur]]),
    temp = temp)
  
  message(glue("Weights predicted for {weight_new}"))
  
  return(df_fish)
}

lastPeriod = function(df_fish){
    period <- colnames(df_fish) %>%
    str_extract("weight[.][0-9]*") %>%
    sub("weight.","",.) %>%
    numextract() %>%
    na.omit() %>%
    max()
    
    return(period)
}

#' funkyTranspose
#' https://stackoverflow.com/questions/6645524/what-is-the-best-way-to-transpose-a-data-frame-in-r-and-to-set-one-of-the-column
#' Credits: mortonjt and nzcoops
funkyTranspose = function(df){
  # Transpose table YOU WANT
  df_t <- t(df[,2:ncol(df)])
  # Set the column headings from the first column in the original table
  colnames(df_t) <- t(df[,1])
  return(df_t)
}

# write_size and feed table

write_size_and_feed_table = function(df_fish,df_temp,date_start,n_weeks)
{
  
  df = data.frame()
  date_now = 0
  
  for (i in 1:n_weeks)
  {
    date_now = date_start + weeks(i)
    message(date_now)
    df_fish_now = df_fish %>% predict_weights(df_temp,date_now)
    period = lastPeriod(df_fish_now)
    realPeriod = lastPeriod(df_fish) 
    df_fish_now$weight = df_fish_now[[glue("weight.{period}")]]
    
    
    proportions_feed = df_fish_now$weight %>% na.omit() %>% calc_feedSizePercentage() %>% funkyTranspose()
    print(proportions_feed)
    tank_biomass = df_fish_now %>% group_by(!!sym(glue("tank.{realPeriod}"))) %>% summarize(biomass=sum(weight,na.rm=T)) %>% na.omit() %>% funkyTranspose()
    print(tank_biomass)
    week_now = data.frame(week=week(date_now))
    
    newrow = cbind(week_now,tank_biomass,proportions_feed)
    df = df %>% rbind(newrow)
    
  }
  
  filename = glue("table - predicted biomass and feed size - {date_start} to {date_now}.csv")
  write.table(df, filename, sep=",", row.names=F)
  message("Table saved to ",filename)
  return(df)

  
}

# Create dates period 3 (predicted period)
# The "Time" variable is the time passed between period 2 and period 3, that is, the length of period 3
#data$date_p3 = today()
#data$time_p3 = data$date_p3 - data$date_p2

#Calculate mean temperatures for period one and three (p2, p3)
#Function "calc_meanTemp" is appplied for all fish (rows)
#data$calc_meanTemp_p2 = mapply(calc_meanTemp, startDate=data$date_p1, endDate=data$date_p2, treatment=data$temp)
#data$calc_meanTemp_p3 = mapply(calc_meanTemp, startDate=data$date_p2, endDate=data$date_p3, treatment=data$temp)

#Calculate Thermal growth coefficient TGC for period p2
#Function "calc_TGC" is applied separately for all fish, afterwards a mean TGC is given for all fish
#data$TGC = calc_TGC(data$weight_p1, data$weight_p2, data$calc_meanTemp_p2, data$time_p2)
#data$TGCmean = mean(data$TGC, na.rm=T)

#Predict the weights at the end of period 3 by using TGC from period 2  
#data$weight_p3Pred = predict_growth(data$weight_p2, data$TGCmean, data$calc_meanTemp_p3, data$time_p3)
#data$weight_p2Pred = predict_growth(data$weight_p1, data$TGCmean, data$calc_meanTemp_p2, data$time_p2)

