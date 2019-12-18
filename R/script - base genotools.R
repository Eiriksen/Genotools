

#' Converts messy names and ID's to tidy clean ones.
#'
#' For sorting out a vector with long and complicated identifiers or row names, where the true ID of a row is hidden in a string.\cr
#' E.g: Make "dirty" ID's like "A0006_3911_BT-F1_GTCGTCTA_run20190930N" turn into "clean" ID's like 3991_BT
#' @param vector A vector of "dirty" IDs
#' @param identifier ID's need to be formated with a number and following identifier, e.g "34_individuals2019" where "_individuals2019" is the identifier. Any entries not matching this format will be removed.
#' @param identifier_left Wether the identifier is on the left hand (T) or right-hand (R) side of the number
#' @param numLength if you want leading zeroes, use this parameter to specify the length of the number, e.g "8" for 00000342
#' @param prefix if you want a prefix in the new cleaned ID. Ex: "individuals2019_" will give you "individuals2019_0034". If not specified, the old identifier will be used instead. Set to NA if you only want the number.
#' @param na_remove if you want to remove any entries that don't follow your pattern (otherwise, they'll turn to NA)
#' @export
clean_ID = function(vector,identifier="", identifier_left=F, numLength=4, prefix, na_remove=F,numeric=F) {
  require(tidyverse)
  require(stringr)

  # SET THE REGULAR EXPRESSION
  if (!identifier_left) regExpr = paste("[0-9]{1,50}",identifier,sep="")
  else                  regExpr = paste(identifier,"[0-9]{1,50}",sep="")

  # Extract the ID's from the dirty ID's
  ID_dirty = vector
  ID_clean = ID_dirty %>% str_extract(regExpr)

  # Remove the old identifier (for now)
  ID_clean = ID_clean %>% sub(identifier,"",.)

  # Remove NA values
  if (na_remove) ID_clean = ID_clean[!is.na(ID_clean)]

  # Add leading zeroes
  if (numLength!=0) ID_clean[!is.na(ID_clean)] = ID_clean[!is.na(ID_clean)] %>% as.numeric() %>% sprintf(paste("%0",numLength,"d",sep=""),.)

  # Make the ID completely numeric
  if (numeric) ID_clean = as.numeric(ID_clean)

  # Add the new prefix
  if (exists("prefix")){
    if (is.na(prefix))       return(ID_clean)
    else                     ID_clean[!is.na(ID_clean)] = paste(prefix, ID_clean[!is.na(ID_clean)], sep="")
  }
  else if (identifier_left)  ID_clean[!is.na(ID_clean)] = paste(ID_clean[!is.na(ID_clean)], identifier, sep="")
  else if (!identifier_left) ID_clean[!is.na(ID_clean)] = paste(identifier, ID_clean[!is.na(ID_clean)], sep="")

  return(ID_clean)
}


#' In a dataframe, converts messy names and ID's to tidy clean ones.
#'
#' For sorting out column with long and complicated identifiers or row names, where the true ID of a row is hidden in a string.\cr
#' E.g: Make "dirty" ID's like "A0006_3911_BT-F1_GTCGTCTA_run20190930N" turn into "clean" ID's like 3991_BT
#' @param df The data frame
#' @param column The name of a column containing dirty IDs
#' @param identifier ID's need to be formated with a number and following identifier, e.g "34_individuals2019" where "_individuals2019" is the identifier. Any entries not matching this format will be removed.
#' @param identifier_left Wether the identifier is on the left hand (T) or right-hand (R) side of the number
#' @param numLength if you want leading zeroes, use this parameter to specify the length of the number, e.g "8" for 00000342
#' @param prefix if you want a prefix in the new cleaned ID. Ex: "individuals2019_" will give you "individuals2019_0034"
#' @param na_remove if you want to remove any rows that don't follow your pattern (otherwise, they'll turn to NA). Default is True.
#' @export
clean_ID_df = function(df, column_name, identifier="", identifier_left=F, numLength=F, prefix, na_remove=T, keep_name=F, numeric=F){
  require(tidyverse)
  require(stringr)

  # Ectract the dirty ID's
  ID_dirty = unlist(df[column_name])

  # Clean the ID
  ID_clean = clean_ID(ID_dirty, identifier, identifier_left, numLength, prefix,numeric=numeric)

  # Insert the cleaned ID's into the column
  df[column_name] = ID_clean

  # Remove NA values
  if (na_remove) df = df %>% remoNA(column_name)

  # Rename the old ID column
  # Check what name to use
  if (keep_name == F) column_name_new = "ID"
  else if (keep_name == T) column_name_new = column_name
  else column_name_new = keep_name
  # Rename the column to "ID"
  df = df %>% rename(!! column_name_new := !! column_name)

  return(df)
}



#' Converting sdy to F or M
#'
#' Used on dataframes, for determining sex based on SDY in a given column
#' @export
#'
determineSex = function(dataframe, column, cutoff) {
  dataframe = dataframe %>% group_by(ID, SEQRUN) %>% mutate(
    sex = SDY_to_sex(dataframe %>% select(matches(column)) %>% filter(dataframe$ID==ID) , cutoff)
  )
  # %>% select(-c(column))
  return(dataframe )
}


#' Set sex to NA if many SNPs missing.
#'
#' In a dataframe, sets sex to "NA" when a certain amount of SNP's are missing as NA
#' @export
#'
unSexBad = function(dataframe, column, sensitivity=0.35) {
  sex = unlist(dataframe[column])
  colNum = length(names(dataframe))

  na_prop <- apply(dataframe, 1, function(x) sum(is.na(x))/length(x))

  sex[na_prop > sensitivity] = "?"

  dataframe$sex = sex
  return(dataframe)
}

#' Rename genotypes based on a lookup table
#'
#' In a dataframe, rename genotype columns
#' @export
renameGenotypes = function(dataframe, LUT, not_genotypes=c()) {
  for (i in names(dataframe %>% select(-c(not_genotypes)))) {
    dataframe <- dataframe %>% renameGenotype(i, LUT)
  }
  dataframe
}

#' determinesex2
#' @keywords internal
determineSex2 = function(dataframe, column, cutoff) {
  dataframe = dataframe %>% group_by(ID) %>% mutate(
    sex = SDY_to_sex(dataframe %>% select(matches(column)) %>% filter(dataframe$ID==ID) , cutoff)
  )
  # %>% select(-c(column))
  return(dataframe )
}

#' SDY_to_sex
#' @keywords internal
SDY_to_sex = function(vector, cutoff) {
  sdy = mean(unlist(vector[1]), na.rm=T)

  if (is.na(sdy)) return(NA)
  else if (sdy <= cutoff) return("F")
  else return("M")
}

#' safeMerge
#' @keywords internal
safeMerge = function(vector){
  # Get the datatype of the vector
  type = typeof(vector)

  #1 remove NA values
  vector = vector[!is.na(vector)]
  #check if the remaning entries are equal

  #if they are, return one of them
  #if they're not, return NA

  if (length(unique(vector)) == 1) return(unique(vector))
  else return(convertType(NA,type))
}

#' renameGenotype
#' @keywords internal
renameGenotype = function(dataframe, column, LUT=c("1"="1 1","2"="1 2","3"="2 2")){
  genotype = dataframe[column] %>% unlist()

  col = LUT[genotype]
  col[is.na(col)] = "* *"
  dataframe[column] = col

  return(dataframe)
}

#' Cecks if certain columns exist in a dataset and returns an error message if not
#' @keywords internal
check_columns = function(dataset,columns,preMessage="Missing columns:"){
  message = c(preMessage)

  for (i in columns) {
    if(!i %in% colnames(dataset))
    {
      message = c(message,paste("Column",i,"is missing."))
    }

    if(length(message)>1) error(message)
  }
}

#' Converts anything to a number
#' @export
numextract <- function(string){
  require(stringr)
  as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

#' Removes rows with NA in a given column
#'
#' Removes NA rows (in a given column) from a dataset
#' @export
remoNA = function(dataset,column){
  return(dataset[which(!is.na(dataset[column])),])
}

#' Replace NA values with unique new identifier
#'
#' Changes all NA values to an unique identifier
#' @export
uniNA = function(values){
  uniques = cumsum(is.na(values))
  for (i in 1:length(values)){
    if (is.na(values[i])) {
      values[i]=paste("NA-",uniques[i],sep="")
    }
  }
  return(values)
}


#' Changes NA values in a dataframe to a given value
#' @export
  changeNA = function(dataset,value){
  dataset[is.na(dataset)] = value
  return(dataset)
}


#' Changes certain values in a list/vector to NA
#' @export
makeNA = function(values, which){
  for (i in which){
  values[values==i] = NA
  }
  return(values)
}



#' Converts a variabe from one type to another
#' @export
convertType = function(var,type){ #https://stackoverflow.com/questions/47410679/change-type-of-object-based-on-typeof-another-object
  unlist(lapply(var,paste0('as.',type)))
}

unSelect = function(df,...){
  return(df %>% select(-c(...)))
}

#' For looking up variables from one dataset and then add them to another one.
#'
#' Use to add a column (value) to a dataset (samples), from another dataset (lookup), based on an identifier that exists in both (id_lookup)
#' @param df_samples samples to look up
#' @param df_lookup dataframe to look up against
#' @param id_column common column between the two sets containing unique identifiers for rows
#' @param value the value that is looked up and added to df_samples
#' @example fishies <- fishies %>% lookup(df_birthdays, "fish_ID", "date_birth")
#' @export
lookup = function(df_samples, df_lookup, id_column, value_column,default=NA,overwrite=T){
  message("Looking up ",value_column," using ",id_column,"...")
  #check if the df_samples already has a column with /value/
  #if not, create one and fill it with NA
  if(!value_column %in% colnames(df_samples)){
    df_samples[[value_column]] = default
  }

  # for each row,
  # 1. get all matching rows (r_match)
  # 2. remove  all rows containing NA values
  # 3. select the first of these rows, and use the value from that one
  # 4. if the value is NA, use the value already present
  # 5. if there already is a value in the row that's being looked up, only overwrite if overwrite==T

  values = apply(df_samples, MARGIN=1, FUN=function(x){
    s_id = x[[id_column]]
    r_match = df_lookup %>%
              filter(!!sym(id_column)==s_id) %>%
              remoNA(id_column)

    # check if this item was found (and is not NA) (and the original is empty or overwrite==T)
    if (nrow(r_match)!=0 & !is.na(r_match[[value_column]][1]) & (is.na(x[[value_column]][1]) | overwrite==T)){
      r_match[[value_column]][1] %>% unlist()
    }
    else{
      #if not, use the value already present
      x[[value_column]][1] %>% unlist()
    }
  })
  message(typeof(values))
  df_samples[[value_column]] = values
  message("Done!")
  df_samples
}



#' Applies a function on the column of a dataframe and then returns that dataframe
#'
#' @param df A dataframe
#' @param column The name of the column (string) that we want apply the function to
#' @param fun The function we use on the column
#' @example dataframe2 <- dataframe1 %>% manipulate("lengths",convertInches)
#' @export
manipulate = function(df, column, fun){
  df[[column]] = fun(df[[column]])
  return(df)
}

#' perform
#' @keywords internal
perform = function(df, column, fun){
  return(fun(df[[column]]))
}

#' Selects random rows from a dataframe
#'
#' Takes a dataframe and a number n
#' Returns return n randomly selected rows from the dataframe (as a dataframe)
#' @export
selRandom = function(df, n) {
  rows = round(runif(n,0,nrow(df)))
  selection = df[rows,]
  return(selection)
}

