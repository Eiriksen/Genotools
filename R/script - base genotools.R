
#' clean_ID()
#'
#' For cleaning "dirty" ID's and removing rows that doesn't match our ID-pattern
#' @param dataset the dataset we're using. Each row has it's own ID
#' @param column the name of the column containing the ID's
#' @param identifier ID's need to be formated with a number and following expresio, e.g "34_individuals2019" where "_individuals2019" is the expression. Any entries not matching this format will be removed.
#' @param trailing_ident Wether the expression if before (F) or after (T)of the ID number
#' @param numLength if you want leading zeroes, use this parameter to specify the length of the number, e.g "8" for 00000342
#' @param prefix if you want a prefix in the new cleaned ID. Ex: "individuals2019_" will give you "individuals2019_0034"s
#' @export
clean_ID = function(dataset, column, identifier, trailing_ident=F, numLength=0, prefix="", numeric=F)
{
  # Extract the dirty ID's
  dirtyID = unlist(dataset[column])
  # set the regular identifier to be used based on the "trailing" parameter
  if (trailing_ident) regExpr = paste("[0-9]*",identifier,sep="")
  else          regExpr = paste(identifier,"[0-9]*",sep="")
  # Creat the clean ID using str_extract on the dirtyID together with the regular identifier
  cleanID = dirtyID %>% str_extract(regExpr)
  
  # Set the old column to be the new ID
  dataset[column] = cleanID
  # Rename the column to "ID", and remove NA values (thoes not fitting the format)
  dataset = dataset %>% rename("ID"=column) %>% remoNA("ID")
  
  # Remove the old identifier
  dataset$ID = dataset$ID %>% sub(identifier,"", .)
  
  # Add leading zeroes
  if (numLength !=0) dataset$ID = dataset$ID %>% as.numeric() %>% sprintf( paste("%0",numLength,"d",sep=""), .)
  
  # Make numeric (or not)
  if (numeric) dataset$ID = as.numeric(dataset$ID)
  
  # Add the new prefix
  if (prefix!="") dataset$ID =dataset$ID %>% paste(prefix, ., sep="")
  
  return(dataset)
}


#' determineSex
#' For determining sex based on SDY
#' @export
#' 
determineSex = function(dataframe, column, cutoff)
{
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
unSexBad = function(dataframe, column, sensitivity=0.35)
{
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


determineSex2 = function(dataframe, column, cutoff)
{
  dataframe = dataframe %>% group_by(ID) %>% mutate(
    sex = SDY_to_sex(dataframe %>% select(matches(column)) %>% filter(dataframe$ID==ID) , cutoff)
  )
  # %>% select(-c(column))
  return(dataframe )
}

SDY_to_sex = function(vector, cutoff)
{
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

#' convertType
#' Converts a variabe from one type to another
#' @export
convertType = function(var,type){ #https://stackoverflow.com/questions/47410679/change-type-of-object-based-on-typeof-another-object
  unlist(lapply(var,paste0('as.',type)))
}

unSelect = function(df,...){
  return(df %>% select(-c(...)))
}
