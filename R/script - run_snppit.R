#' run_snppt
#'
#' Ver 0.1
#'
#' Function that simplifies using SNPPT.exe to calculate parentage based on SNP data
#' Requires that snppt.exe is installed in the working directory!
#'
#' Takes two datasets, one for offspring and one for parents. \cr
#' Both datasets must be formated as folows:
#' \itemize{
#' \item Individuals as rows
#' \item SNP data on columns, one column pr SNP
#' \item SNP genotype written numerically: 1=homozygous for allele a, 2=heterozygous, 3=homozygous for allele b
#' \item One column named "ID" that contains some ID for each individual,
#' \item Only for parents: One column named "population" that's either "M" or "F"
#' \item Only for parents: One column named "population" that gives that individual's population
#' }
#'
#' Important!: Any column that is not ID, sex (parents) or population (parents) will be interpreted as an SNP!!!
#'
#' @param data_offspring Dataset contaning offspring ID's and SNP genotypes (see above)
#' @param data_parents Dataset containing parent ID's, sex, population, and SNP genotypes (see above)
#' @param projectName Optional. A name that will be used for files created during the proces
#' @export
#' @examples run_snppt(offspring, parents, "Project_oct2019")
run_snppit <- function(data_offspring, data_parents, projectName="project1"){

  # Check that the parent set has columns "ID", "Sex" and "population"
  check_columns(data_parents,c("ID","sex","population"))
  check_columns(data_offspring,c("ID"))

  # Convert the genotype data in the dataset from normal numeric to snppt numeric type
  data_parents = data_parents %>%
    renameGenotypes(LUT=c("1"="1 1","2"="1 2","3"="2 2"), not_genotypes=c("ID","sex","population"))
  data_offspring = data_offspring %>%
    renameGenotypes(LUT=c("1"="1 1","2"="1 2","3"="2 2"), not_genotypes=c("ID"))

  # remove any column that is not equal between parents and offspring dataset (except population and sex)
  data_parents   <- data_parents %>% select( "ID","population","sex", data_offspring %>% names() %>% one_of() )
  data_offspring <- data_offspring %>% select( "ID", data_parents %>% names() %>% one_of() )

  #Filename to use for snpptfile
  snpptfile_name = paste("snpptfile -", projectName)

  #Write SNPPIT file based on parent and offspring data
  message("Attempting to write SNPPT settings file...")
  SNPPITfile(snpptfile_name,data_parents, data_offspring)
  message("SNPPT settings file written!")

  # Run snppit
  message("Starting SNPPT analysis...")
  location = paste("'",getwd(),"'",sep="") %>% str_replace_all("/","\\\\")
  system2("powershell", args=c("Set-location",location))
  system2("powershell", args=c(paste("./snppit -f '",snpptfile_name,".txt' --max-par-miss 100",sep="")))

  # Obtain snppit results and clean them
  message("Trying to obtain SNPPT results...")
  data_snppit = read.table("snppit_output_ParentageAssignments.txt", head=T, comment.char = "") %>%
    rename(ID_offspring=Kid, ID_pa=Pa, ID_ma=Ma, population=PopName)

  # Return results
  message("results ready!")
  return(data_snppit)
}



SNPPITfile = function(name,parents, offspring,parentSex=T,parentPopulation=T)
{
  # ASSUMES the following names for columns:
  # IDs:        "ID"         (string, no whitespace)
  # Sex:        "sex"        (F/M/?)
  # Population: "population" (string, no whitespace)
  # And no other columns other than loci!

  printLines = c()

  loci = names(parents)
  loci = loci [! loci %in% c("ID")]

  if (parentPopulation)
  {
    loci = loci[! loci %in% c("population")]
  }

  if (parentSex)
  {
    loci = loci[! loci %in% c("sex")]
    printLines = c(printLines, "POPCOLUMN_SEX")
  }

  printLines = c(paste("NUMLOCI", length(loci)),"MISSING_ALLELE *",printLines)



  for( i in loci)
  {
    printLines = c(printLines, paste(i, "0.005"))
  }

  snppitFile = file(paste(name,".txt",sep=""), "w")
  writeLines(printLines, snppitFile)

  if (parentPopulation)
  {
    for (i in unique(parents$population))
    {
      writeLines(paste("POP",i),snppitFile)
      write.table(parents[which(parents$population==i),] %>% select(-c("population")), snppitFile, sep="\t", quote=F, row.names = F, col.names=F)
    }
  }
  else
  {
    writeLines("POP Parentpool_0",snppitFile)
    write.table(parents, snppitFile, sep="\t", quote=F, row.names = F, col.names=F)
  }


  writeLines("OFFSPRING Offspring_pool ?",snppitFile)
  write.table(offspring, snppitFile, sep="\t", quote=F, row.names = F, col.names=F)

  close(snppitFile)
}
