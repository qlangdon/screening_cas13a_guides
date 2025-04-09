#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
library(tidyverse)

args <- commandArgs(TRUE)
inputPrefix <- args[1]
inputCutoff <- args[2]

guidesIn <- read_delim(paste0(inputPrefix, "_guides_filtered",inputCutoff, ".tsv"))
if ("refStart" %in% colnames(guidesIn)) {
  guidesIn <- guidesIn |> rename(start = refStart)
}

foldedIn <- read_delim("score_RNAfold_crRNAs.txt")
foldedKeep <- foldedIn |>
  select(target, has_hairpin, spacer_basepairs, MFE)

guidesOut <- full_join(guidesIn, foldedKeep) |>
  relocate(has_hairpin, .after =target)

write_tsv(guidesOut, paste0(inputPrefix, "_guides_withFoldScore.tsv"))

guidesFiltered <- guidesOut |>
   filter(has_hairpin==T)

print(paste("Kept", nrow(guidesFiltered), "of input", nrow(guidesOut), nrow(guidesFiltered)/nrow(guidesOut)))

write_tsv(guidesFiltered, paste0(inputPrefix, "_guides_withHairpin.tsv"))

guidesOutWindows <- guidesFiltered |>
  select(accession, start, target, spacer, strand, GC_content, A_content) |>
  rename(segment=accession)

write_tsv(guidesOutWindows, paste0(inputPrefix, "_guides_windows_withHairpin.txt"))
write_tsv(guidesOutWindows, "windows.txt")

guidesOutSpacer <- guidesOutWindows |>
  select(spacer)

write_tsv(guidesOutSpacer, paste0(inputPrefix, "_guides_spacers_withHairpin.txt"))
write_tsv(guidesOutSpacer, "spacers.txt", col_names = F)

print(paste0("Getting fasta of ", nrow(guidesOutWindows)))
guidesOutTargets <- c()
for (i in 1:nrow(guidesOutWindows)) {
  result1 <- paste0(">", guidesOutWindows$segment[i], "_", guidesOutWindows$start[i])
  result2 <- guidesOutWindows$target[i]
  guidesOutTargets <- c(guidesOutTargets, result1, result2)
  #print(paste(i, "of", nrow(guidesOutWindows)))
}

print("Writing fasta")
writeLines(guidesOutTargets, paste0(inputPrefix, "_guides_targets_withHairpin.fa"))
writeLines(guidesOutTargets, "targets.fa")
print("Done")



