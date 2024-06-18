#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
library(tidyverse)
library(readxl)

#args <- commandArgs(TRUE)
#inputDF <- args[1]
#inputDF <- "IAV_all_guides_withHairpin_wHumanAlignZero_kept0.9.tsv"

guide_key <- read_xlsx("~/Dropbox (Gladstone)/screening_cas13a_guides/data/guides/HIV_set7.xlsx") %>% 
  arrange(start) %>%
  filter(Guide!="HIV_cr3") %>%
  select(Guide, 'Guide sequence', gene, start, end) 

#nameSplit <- str_split(inputDF, "_guides")
#spPrefix <- paste(nameSplit[[1]][1], nameSplit[[1]][2], sep="_")
#spPrefix <- nameSplit[[1]][1]

#guidesDF <- read_delim(inputDF)
#guidesDF <- guidesDF |> mutate(target=toupper(target))

#if ("refStart" %in% colnames(guidesDF)) {
#  guidesDF <- guidesDF |> rename(start = refStart)
#}

hammingDF <- tibble()

for (i in 1:nrow(guide_key)) {
  tempRow <- guide_key |>
    select(Guide, `Guide sequence`, start) |>
    slice(i)
  header <- colnames(tempRow)
  guideFile <- paste0(tempRow$Guide, ".tsv")
  for (testID in c("neighbors_vn", "neighbors_rb", "pathogen", "standard")) {
    krakenSumIn <- read_delim(paste0("~/Dropbox (Gladstone)/screening_cas13a_guides/results/guide_screening_pipeline/3_classified_", testID, "/", guideFile))
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

guidesOut <- full_join(guide_key, hammingDF)
write_tsv(guidesOut, paste0("HIV-1_guides_wHammingSummary.tsv"))

