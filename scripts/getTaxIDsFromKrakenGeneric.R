options(stringAsFactors=FALSE)
library(tidyverse)
args <- commandArgs(TRUE)

#inputTaxid <- 11320
#inputSp <- "Influenza_A_virus"
#inputTaxid <- 208893
#inputSp <- "RSV-A"
#inputTaxid <- 208895
#inputSp <- "RSV-B"
#inputTaxid <- 3052464
#inputSp <- "Dengue"
inputSp <- args[1]
inputTaxid <- args[2]

#databaseNames <- c("viral", "refBact", "pathogen", "standard")
databaseNames <- c("viral-neighbors-20200605-wo", "refseq-bacteria-20200607", "pathogen_strains_20200807", "standard_db_20231214")
allSpHits <- tibble()
for (dataID in databaseNames) {
  print(dataID) 
  currentDB <- read_delim(paste0(dataID, "/taxDB"), col_names = F)
  colnames(currentDB) <- c('taxID', "parentID", "taxa", "rank")
  spRow <- currentDB |>
    filter(taxID == inputTaxid)
  childrenRows <- currentDB |>
    filter(parentID == inputTaxid)
  spRows <- spRow |>
    add_case(childrenRows)
  preRows <- childrenRows
  continue <- T
  while (continue==T) {
    print(dim(spRows))
    nextChildrenRows <- currentDB |>
      filter(parentID %in% preRows$taxID)
    if (nrow(nextChildrenRows)>0) {
      spRows <- spRows |>
        add_case(nextChildrenRows)
      preRows <- nextChildrenRows
    } else {
      continue <- F
    }
  }
  write_delim(spRows, paste0(inputSp, "_", dataID, "_taxIDs.tsv"))
  if (nrow(allSpHits)==0) {
    allSpHits <- spRows
  } else {
    allSpHits <- allSpHits |>
      add_case(spRows)
  }
}

allChildrenRows <- allSpHits |>
  filter(taxID != inputTaxid)
 
uniqueTaxIDs <- allChildrenRows |>
  select(taxID) |>
  distinct(taxID)

uniqueParentIDs <- allChildrenRows |>
  select(parentID) |>
  distinct(parentID)

allIDs <- uniqueTaxIDs |>
  add_case(uniqueParentIDs |> rename(taxID=parentID)) |>
  distinct(taxID) |>
  arrange(taxID)

write_delim(allIDs, paste0(inputSp, "_allDatabase_taxIDs.tsv"))


