#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
library(tidyverse)

args <- commandArgs(TRUE)
inputDF <- args[1]
#inputDF <- "IAV_all_guides_withHairpin_wHumanAlignZero_kept0.9.tsv"

nameSplit <- str_split(inputDF, "_")
spPrefix <- paste(nameSplit[[1]][1], nameSplit[[1]][2], sep="_")

guidesDF <- read_delim(inputDF)
guidesDF <- guidesDF |> mutate(target=toupper(target))

if ("refStart" %in% colnames(guidesDF)) {
  guidesDF <- guidesDF |> rename(start = refStart)
}

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
      filter(rank!="sequence")
    for (j in 0:4) {
      tempRow <- tempRow |>
        add_column(count = nrow(krakenSumIn |> 
                                  filter(hamming_dist==j)))
      header <- c(header, paste0(testID, j))
      colnames(tempRow) <- header
      if (j==0) {
        namesZero <- NA
        if (nrow(krakenSumIn |> filter(hamming_dist==j)) > 0) {
          zeroRows <- krakenSumIn |> 
            filter(hamming_dist==0) |>
            unite(taxWrank, taxName:rank, sep="-") |>
            select(taxWrank) |>
            summarise(names = paste0(taxWrank,collapse = ", "))
          namesZero <- zeroRows$names
        } 
      }
    }
    tempRow <- tempRow |>
      add_column(names = namesZero)
    header <- c(header, paste0(testID, "Names0"))
    colnames(tempRow) <- header
  }
  hammingDF <- rbind(hammingDF, tempRow)
}

guidesOut <- full_join(guidesDF, hammingDF)
write_tsv(guidesOut, paste0(spPrefix, "_guides_wHammingSummary.tsv"))

print(spPrefix)
print("Viral 0 = 0")
print(paste(nrow(hammingDF |> filter(viral0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viral0==0))/nrow(guidesDF)))
print("Viral 0 = 1")
print(paste(nrow(hammingDF |> filter(viral0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viral0==1))/nrow(guidesDF)))

print("Viral 1 = 0")
print(paste(nrow(hammingDF |> filter(viral1==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viral1==0))/nrow(guidesDF)))
print("Viral 1 = 1")
print(paste(nrow(hammingDF |> filter(viral1==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(viral1==1))/nrow(guidesDF)))

print("bacteria 0 = 0")
print(paste(nrow(hammingDF |> filter(bacteria0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteria0==0))/nrow(guidesDF)))
print("bacteria 0 = 1")
print(paste(nrow(hammingDF |> filter(bacteria0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteria0==1))/nrow(guidesDF)))
print("bacteria 1 = 0")
print(paste(nrow(hammingDF |> filter(bacteria1==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteria1==0))/nrow(guidesDF)))
print("bacteria 1 = 1")
print(paste(nrow(hammingDF |> filter(bacteria1==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(bacteria1==1))/nrow(guidesDF)))

print("pathogen 0 = 0")
print(paste(nrow(hammingDF |> filter(pathogen0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogen0==0))/nrow(guidesDF)))
print("pathogen 0 = 1")
print(paste(nrow(hammingDF |> filter(pathogen0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogen0==1))/nrow(guidesDF)))
print("pathogen 1 = 0")
print(paste(nrow(hammingDF |> filter(pathogen1==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogen1==0))/nrow(guidesDF)))
print("pathogen 1 = 1")
print(paste(nrow(hammingDF |> filter(pathogen1==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(pathogen1==1))/nrow(guidesDF)))


print("standard 0 = 0")
print(paste(nrow(hammingDF |> filter(standard0==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standard0==0))/nrow(guidesDF)))
print("standard 0 = 1")
print(paste(nrow(hammingDF |> filter(standard0==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standard0==1))/nrow(guidesDF)))
print("standard 1 = 0")
print(paste(nrow(hammingDF |> filter(standard1==0)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standard1==0))/nrow(guidesDF)))
print("standard 1 = 1")
print(paste(nrow(hammingDF |> filter(standard1==1)), "of", nrow(guidesDF), 
            "so", nrow(hammingDF |> filter(standard1==1))/nrow(guidesDF)))