# Stuff to add later to clean_ID
#   option to not remove NA values
#   alternative function that only works on vectors/lists

#' clean_ID
#'
#' For cleaning "dirty" ID's and removing rows that doesn't match our ID-pattern
#' @param dataset the dataset we're using. Each row has it's own ID
#' @param column the name of the column containing the ID's
#' @param identifier ID's need to be formated with a number and following expresio, e.g "34_individuals2019" where "_individuals2019" is the expression. Any entries not matching this format will be removed.
#' @param trailing_ident Wether the expression if before (F) or after (T)of the ID number
#' @param numLength if you want leading zeroes, use this parameter to specify the length of the number, e.g "8" for 00000342
#' @param prefix if you want a prefix in the new cleaned ID. Ex: "individuals2019_" will give you "individuals2019_0034"s
#' @param keepName T: keeps the old column name, F: renames it to "ID", any other string: Renames the column to this string
#' @export
clean_ID = function(dataset, column, identifier="", trailing_ident=F, numLength=4, prefix="", numeric=F, keepName=F)
{
  # Extract the dirty ID's
  dirtyID = unlist(dataset[column])
  # set the regular identifier to be used based on the "trailing" parameter
  if (trailing_ident) regExpr = paste("[0-9]{1,50}",identifier,sep="")
  else          regExpr = paste(identifier,"[0-9]{1,50}",sep="")
  # Creat the clean ID using str_extract on the dirtyID together with the regular identifier
  cleanID = dirtyID %>% str_extract(regExpr)

  # Set the old column to be the new ID
  dataset[column] = cleanID

  # Check what name to use
  if (keepName == F) nColName = "ID"
  else if (keepName == T) nColName = column
  else nColName = keepName

  # Rename the column to "ID", and remove NA values (those not fitting the format)
  dataset = dataset %>% rename(!! nColName := !! column) %>% remoNA(nColName)

  # Remove the old identifier
  dataset[[nColName]] = dataset[[nColName]] %>% sub(identifier,"", .)

  # Add leading zeroes
  if (numLength !=0) dataset[[nColName]] = dataset[[nColName]] %>% as.numeric() %>% sprintf( paste("%0",numLength,"d",sep=""), .)

  # Make numeric (or not)
  if (numeric) dataset[[nColName]] = as.numeric(dataset[[nColName]])

  # Add the new prefix
  if (prefix!="") dataset[[nColName]] =dataset[[nColName]] %>% paste(prefix, ., sep="")

  return(dataset)
}



#' determineSex
#' For determining sex based on SDY
#' @export
#'
determineSex = function(dataframe, column, cutoff) {
  dataframe = dataframe %>% group_by(ID, SEQRUN) %>% mutate(
    sex = SDY_to_sex(dataframe %>% select(matches(column)) %>% filter(dataframe$ID==ID) , cutoff)
  )
  # %>% select(-c(column))
  return(dataframe )
}


#' unSexBad
#'
#' Sets sex to "NA" when a certain amount of SNP's are missing as NA
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

#' renameGenotypes
#'
#' rename genotype columns
#' @export
renameGenotypes = function(dataframe, LUT, not_genotypes=c()) {
  for (i in names(dataframe %>% select(-c(not_genotypes)))) {
    dataframe <- dataframe %>% renameGenotype(i, LUT)
  }
  dataframe
}


determineSex2 = function(dataframe, column, cutoff) {
  dataframe = dataframe %>% group_by(ID) %>% mutate(
    sex = SDY_to_sex(dataframe %>% select(matches(column)) %>% filter(dataframe$ID==ID) , cutoff)
  )
  # %>% select(-c(column))
  return(dataframe )
}

SDY_to_sex = function(vector, cutoff) {
  sdy = mean(unlist(vector[1]), na.rm=T)

  if (is.na(sdy)) return(NA)
  else if (sdy <= cutoff) return("F")
  else return("M")
}


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



renameGenotype = function(dataframe, column, LUT=c("1"="1 1","2"="1 2","3"="2 2")){
  genotype = dataframe[column] %>% unlist()

  col = LUT[genotype]
  col[is.na(col)] = "* *"
  dataframe[column] = col

  return(dataframe)
}

# Cecks if certain columns exist in a dataset and returns an error message if not
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

#' numextract
#' @export
numextract <- function(string){
  require(stringr)
  as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

#' remoNA
#' Removes NA rows (in a given column) from a dataset
#' @export
remoNA = function(dataset,column){
  return(dataset[which(!is.na(dataset[column])),])
}

#' uniNA
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

#' changeNA
#' Changes NA values in a dataframe to a given value
#' @export
  changeNA = function(dataset,value){
  dataset[is.na(dataset)] = value
  return(dataset)
}

#' makeNA
#' Changes certain values in a list/vector to NA
#' @export
makeNA = function(values, which){
  for (i in which){
  values[values==i] = NA
  }
  return(values)
}



#' convertType
#' Converts a variabe from one type to another
#' @export
convertType = function(var,type){ #https://stackoverflow.com/questions/47410679/change-type-of-object-based-on-typeof-another-object
  unlist(lapply(var,paste0('as.',type)))
}

unSelect = function(df,...){
  return(df %>% select(-c(...)))
}

#' lookup
#'
#' For looking up variables from one dataset and then add them to another one.
#' Use to add a column (value) to a dataset (samples), from another dataset (lookup), based on an identifier that exists in both (id_lookup)
#' @param df_samples samples to look up
#' @param df_lookup dataframe to look up against
#' @param id_column common column between the two sets containing unique identifiers for rows
#' @param value the value that is looked up and added to df_samples
#' @example fishies <- fishies %>% lookup(df_birthdays, "fish_ID", "date_birth")
#' @export
lookup = function(df_samples, df_lookup, id_column, value_column,default=NA){
  message("Looking up ",value_column," using ",id_column,"...")
  #check if the df_samples already has a column with /value/
  #if not, create one and fill it with NA
  if(!value_column %in% colnames(df_samples)){
    df_samples[[value_column]] = default
  }

  values = apply(df_samples, MARGIN=1, FUN=function(x){
    s_id = x[[id_column]]
    r_match = df_lookup %>% filter(!!sym(id_column)==s_id)
    # check if this item was found (and is not NA)
    if (nrow(r_match)!=0 & !is.na(r_match[[value_column]][1])){
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


#' manipulate
#' Applies a function on the column of a dataframe and then returns that dataframe
#' @param df A dataframe
#' @param column The name of the column (string) that we want apply the function to
#' @param fun The function we use on the column
#' @example dataframe2 <- dataframe1 %>% manipulate("lengths",convertInches)
#' @export
manipulate = function(df, column, fun){
  df[[column]] = fun(df[[column]])
  return(df)
}


