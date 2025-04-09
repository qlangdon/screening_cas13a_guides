#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
suppressMessages(library(Biostrings, quietly=TRUE))
suppressMessages(library(tidyverse, quietly=TRUE))

args <- commandArgs(TRUE)
spName <- args[1]
inputDF <- args[2]
fastaInput <- args[3]
#fastaInput <- "cladeH1_targets.fa"
#inputDF <- "measles_genome_cladeH1_guides.tsv"
#spName <- "Measles"

targetsIn <- readDNAStringSet(fastaInput)
targetsSeqs <- list(names(targetsIn), unname(as.character(targetsIn)))[[2]]
targetsNames <- list(names(targetsIn), unname(as.character(targetsIn)))[[1]]
targetDF <- tibble(guideName=targetsNames, target=targetsSeqs)

nameSplit <- str_split(inputDF, "_guides")
#spPrefix <- paste(nameSplit[[1]][1], nameSplit[[1]][2], sep="_")
spPrefix <- nameSplit[[1]][1]
#spName <- str_split(spPrefix, "_all")[[1]][1]

guidesDF <- read_delim(inputDF, show_col_types = FALSE)
guidesDF <- guidesDF |> mutate(target=toupper(target))
if ("refStart" %in% colnames(guidesDF)) {
  guidesDF <- guidesDF |> rename(start = refStart)
}
guidesCombo <- full_join(targetDF, guidesDF)

tax_ids <- read_delim(paste0("references/", spName, "_allDatabase_taxIDs.tsv"), delim="\t", show_col_types = FALSE)

print(spPrefix)

hammingDF <- tibble()

for (i in 1:nrow(guidesCombo)) {
  tempRow <- guidesCombo |>
    select(guideName, accession, start, target) |>
    slice(i)
  header <- colnames(tempRow)
  guideFile <- paste0(tempRow$guideName, ".tsv")
  for (testID in c("viral", "bacteria", "pathogen", "standard")) {
    krakenSumIn <- read_delim(paste0("3_classified_", testID, "/", guideFile), show_col_types = FALSE)
    krakenSumIn <- krakenSumIn |>
      filter(rank!="sequence") |>
      mutate(Target = ifelse(taxid %in% tax_ids$taxID, "On", "Off"))
    for (targetID in c("On", "Off")) {
      krakenSumTar <- krakenSumIn  |>
        filter(Target==targetID)
      for (j in 0:4) {
        tempRow <- tempRow |>
          add_column(count = nrow(krakenSumTar |> 
                                    filter(hamming_dist==j)))
        header <- c(header, paste0(testID, targetID, j))
        colnames(tempRow) <- header
        if (j==0) {
          namesZero <- NA
          if (nrow(krakenSumTar |> filter(hamming_dist==j)) > 0 & targetID=="Off") {
            zeroRows <- krakenSumTar |> 
              filter(hamming_dist==0) |>
              unite(taxWrank, taxName:rank, sep="-") |>
              select(taxWrank) |>
              summarise(names = paste0(taxWrank,collapse = ", "))
            namesZero <- zeroRows$names
          } 
        }
      }
    }
    tempRow <- tempRow |>
      add_column(names = namesZero)
    header <- c(header, paste0(testID, "NamesOff0"))
    colnames(tempRow) <- header
  }
  hammingDF <- rbind(hammingDF, tempRow)
}

guidesOut <- full_join(guidesCombo, hammingDF)
write_tsv(guidesOut, paste0(spPrefix, "_guides_wHammingSummary.tsv"))

print("viral On 0 = 0")
print(paste(nrow(hammingDF |> filter(viralOn0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viralOn0==0))/nrow(guidesDF)))
print("viral On 0 = 1")
print(paste(nrow(hammingDF |> filter(viralOn0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viralOn0==1))/nrow(guidesDF)))

print("viral Off 0 = 0")
print(paste(nrow(hammingDF |> filter(viralOff0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viralOff0==0))/nrow(guidesDF)))
print("viral Off 0 = 1")
print(paste(nrow(hammingDF |> filter(viralOff0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viralOff0==1))/nrow(guidesDF)))

print("bacteria On 0 = 0")
print(paste(nrow(hammingDF |> filter(bacteriaOn0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteriaOn0==0))/nrow(guidesDF)))
print("bacteria On 0 = 1")
print(paste(nrow(hammingDF |> filter(bacteriaOn0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteriaOn0==1))/nrow(guidesDF)))

print("bacteria Off 0 = 0")
print(paste(nrow(hammingDF |> filter(bacteriaOff0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteriaOff0==0))/nrow(guidesDF)))
print("bacteria Off 0 = 1")
print(paste(nrow(hammingDF |> filter(bacteriaOff0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteriaOff0==1))/nrow(guidesDF)))

print("pathogen On 0 = 0")
print(paste(nrow(hammingDF |> filter(pathogenOn0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogenOn0==0))/nrow(guidesDF)))
print("pathogen On 0 = 1")
print(paste(nrow(hammingDF |> filter(pathogenOn0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogenOn0==1))/nrow(guidesDF)))

print("pathogen Off 0 = 0")
print(paste(nrow(hammingDF |> filter(pathogenOff0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogenOff0==0))/nrow(guidesDF)))
print("pathogen Off 0 = 1")
print(paste(nrow(hammingDF |> filter(pathogenOff0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogenOff0==1))/nrow(guidesDF)))

print("standard On 0 = 0")
print(paste(nrow(hammingDF |> filter(standardOn0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standardOn0==0))/nrow(guidesDF)))
print("standard On 0 = 1")
print(paste(nrow(hammingDF |> filter(standardOn0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standardOn0==1))/nrow(guidesDF)))

print("standard Off 0 = 0")
print(paste(nrow(hammingDF |> filter(standardOff0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standardOff0==0))/nrow(guidesDF)))
print("standard Off 0 = 1")
print(paste(nrow(hammingDF |> filter(standardOff0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standardOff0==1))/nrow(guidesDF)))

