#' get_advanced_genotypes
#' @export
get_advanced_genotypes <- function(offspring, families) {
  # Function for this specific pipeline, gives the correct heterozygote phenotype
  # apply following function to all offspring:

  apply(offspring,1,FUN=function(x) {
    vgll3geno = as.numeric(x[["c25_1441_SAC"]])
    fam = as.numeric(x[["family"]])
    pop = x[["population"]]
    if(is.na(vgll3geno)) return("NA")

    if (vgll3geno == 1) {
      return("EE")
    }
    if (vgll3geno == 2) {
      # find the genotype of the mother
      if (is.na(fam)) return("EL/LE")
      row_fam <- families %>% filter(ID_family == fam & Population==pop)
      damVGLL3 = row_fam$Dam.VGLL3
      if (damVGLL3 == "LL") return("LE")
      else if (damVGLL3 == "EE") return("EL")
      else return("EL/LE")
    }
    if (vgll3geno == 3) {
      return("LL")
    }
    return("none of the above")
  })

}

#' findFamily
#' @export
findFamily = function(ID_mams, ID_paps)
{
  family = unlist(data_families %>% filter(ID_ma == ID_mams & ID_pa == ID_paps) %>% select(ID_family))
  if (family %>% length() == 0) return(NA) else return(family)
}

#' find_familyID
#' @export
find_familyID = function(ID_mam, ID_pap, df_families){
  # df_families must be formated as
  # rows pr family
  # columns: ID_ma, ID_pa, ID_family

  family = df_families %>%
    filter(ID_ma == ID_mam & ID_pa == ID_pap) %>%
    select(ID_family) %>%
    unlist()
  if ( length(family) == 0) return(NA) else return(family)
}

#' find_familyIDs
#' @export
find_familyIDs = function(df, df_families){

  df$ID_family = apply(df,MARGIN=1,FUN=function(x){
    ID_pa = x[["ID_pa"]]
    ID_ma = x[["ID_ma"]]
    find_familyID(ID_ma, ID_pa, df_families)
  })

  df
}

#' get_familyInfo
#' @export
get_familyInfo = function(df,df_families,columns)
{
  for (i in columns){
    df[[i]] = apply(df,MARGIN=1,FUN=function(x){
      ID_fam = x[["ID_family"]]
      df_families %>% filter(ID_family == ID_fam) %>% select(i) %>% as.character() %>% unlist()
    })
  }
  df
}
