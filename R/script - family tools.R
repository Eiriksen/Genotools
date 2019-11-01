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
