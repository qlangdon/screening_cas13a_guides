#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
library(tidyverse)

args <- commandArgs(TRUE)
inputDF <- args[1]
#inputDF <- "IAV_all_guides_withHairpin_wHumanAlignZero_kept0.9.tsv"

nameSplit <- str_split(inputDF, "_guides")
#spPrefix <- paste(nameSplit[[1]][1], nameSplit[[1]][2], sep="_")
spPrefix <- nameSplit[[1]][1]
spName <- str_split(spPrefix, "_all")[[1]][1]

guidesDF <- read_delim(inputDF)
guidesDF <- guidesDF |> mutate(target=toupper(target))
if ("refStart" %in% colnames(guidesDF)) {
  guidesDF <- guidesDF |> rename(start = refStart)
}

tax_ids <- read_delim(paste0("references/", spName, "_allDatabase_taxIDs.tsv"), delim="\t")


hammingDF <- tibble()

for (i in 1:nrow(guidesDF)) {
  tempRow <- guidesDF |>
    select(accession, start, target) |>
    slice(i)
  header <- colnames(tempRow)
  guideFile <- paste0(spPrefix, "_", tempRow$accession, "_", tempRow$start, ".tsv")
  for (testID in c("viral", "bacteria", "pathogen", "standard")) {
    krakenSumIn <- read_delim(paste0("3_classified_", testID, "/", guideFile))
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

guidesOut <- full_join(guidesDF, hammingDF)
write_tsv(guidesOut, paste0(spPrefix, "_guides_wHammingSummary.tsv"))

print(spPrefix)
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

