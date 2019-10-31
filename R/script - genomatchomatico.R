

#' calc_genetic_similarity
#'
#' takes two vectors of SNPs,
#' Each element in the vector should be the genotypes of one SNP, written as bases, e.g "CT" "CC" "AT" etc. Order, e.g. "CT" or "TC" is not taken into account.
#' It's important that the order of SNPs in the two vectors are the same.
#' Returns an integer between 1 and 0 (0 = 0% match, 1=100% match)
#' @param SNPs1 the first vector of SNPs
#' @param SNPs2 the second vector of SNPs
#' 
calc_genetic_similarity = function(sample1, sample2,naCutoff=15){
  require(tidyverse)
  require(stringr)
  require(tidyverse)
  
  #Error checking: similar length
  
  # For easier processing using sapply
  df = na.omit(data.frame(sample1,sample2))
  
  #iterate over each element in the vector
  #check how many alleles they have in common, (CT TC = 4, two matches,  CT CC = two matches, CC TT = nomatches), then store that number
  if (nrow(df) < naCutoff) return(0)
  
  matches = apply(df,1,FUN=function(x){
    snp1 = x[["sample1"]]
    snp2 = x[["sample2"]]
    
    # Heart of this function
    m=0
    if( substr(snp1,1,1)==substr(snp2,1,1)) m=m+1
    if( substr(snp1,2,2)==substr(snp2,2,2)) m=m+1
    m
  })
  # Return the proportion of matches vs total potential matches
  sum(matches)/(length(df$sample1)*2)
}


#' find closest match
#'
#' Using the find_genetic_similarity function, takes a table of genetic samples (rows: samples, columns: snp genotypes) and calculates 
#' @param sample the sample to find closest match to
#' @param lookup a dataframe of samples with SNP data (Each row one sample, each column a snp loki). First row must be ID info.
#' 
calculate_distance_to_sample = function(table, sample,notSNP){ #for one vs one
  require(tidyverse)
  SNPtable = table %>% select(-notSNP)
  table$distance = apply(SNPtable,1,FUN=function(x){
    SNPs = unlist(as.character(x))
    calc_genetic_similarity(sample,SNPs)
  })
  return(table)
}

#a=calculate_distance_to_sample(data_genotypes_p1, data_genotypes_p2 %>% filter(ID=="Offsp6050") %>% select(-c("ID","sex")) %>% as.character(), notSNP=c("ID","sex"))
#plot(a$distance)


find_closest_match = function(sample,table,notSNP) { #for one vs many
  require(tidyverse)
  distTable = table %>% calculate_distance_to_sample(sample,notSNP)
  highestDistance = max(distTable$distance)
  if(highestDistance < 0.1) return(NA)
  closestIndividual = distTable %>% filter(distance == highestDistance)
  closestIndividual
}

find_genetic_matches = function(samples,lookup,notSNP){ # for many vs many
  require(beepr)
  require(tidyverse)
  SNPsamples = samples %>% select(-notSNP)
  
  matches = pbapply(SNPsamples,1,FUN=function(x){
    closestMatch = find_closest_match(as.character(x), lookup, notSNP)
    
    matchID = NA
    matchSimilarity= NA
    
    if(!is.na(closestMatch)) {
      matchID = closestMatch$ID
      matchSimilarity = closestMatch$distance
    }
    
    return(list(matchID=matchID,matchSimilarity=matchSimilarity))
  }) 
  
  
  beep(8)
  matches = matches %>%   
    map(as.data.frame) %>%
    bind_rows() %>%
    select(matchID, matchSimilarity)
  
  samples = cbind(samples,matches)
  samples %>% select(c("ID","matchID","matchSimilarity"))
  
}
