#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
library(tidyverse)

args <- commandArgs(TRUE)
inputPrefix <- args[1]
cutoff <- as.numeric(args[2])

guidesIn <- read_delim(paste0(inputPrefix, "_guides_withHairpin.tsv"))

alignIn <- read_delim("alignment_cts_GRCh38_latest_rna.txt")
alignKeep <- alignIn |>
  rename(target=seq) |>
  select(target, ct)

guidesOut <- full_join(guidesIn, alignKeep) |>
  relocate(ct, .before =total_hit) |>
  rename(humanCount = ct)

write_tsv(guidesOut, paste0(inputPrefix, "_guides_withFoldScore_wHumanAlignCount.tsv"))

guidesFiltered <- guidesOut |>
  filter(humanCount==0)

print(paste("Kept", nrow(guidesFiltered), "of input", nrow(guidesOut), nrow(guidesFiltered)/nrow(guidesOut)))

write_tsv(guidesFiltered, paste0(inputPrefix, "_guides_withHairpin_wHumanAlignZero.tsv"))

if ("prop_subtypes_hit" %in% colnames(guidesFiltered)) {
  guidesFiltered <- guidesFiltered |> rename(p_subtypes_hit =prop_subtypes_hit)
}
if ("refStart" %in% colnames(guidesFiltered)) {
  guidesFiltered <- guidesFiltered |> rename(start = refStart)
}
if ("prop_unknown" %in% colnames(guidesFiltered)) {
  guidesFiltered <- guidesFiltered |> rename(p_unknown = prop_unknown )
}

guidesKeep <- guidesFiltered |>
  filter(if_any(starts_with("prop"),  ~ .x > cutoff))

print(paste("Using", cutoff, "Kept", nrow(guidesKeep), "of input", nrow(guidesFiltered), nrow(guidesKeep)/nrow(guidesFiltered)))

write_tsv(guidesKeep, paste0(inputPrefix, "_guides_withHairpin_wHumanAlignZero_kept", cutoff, ".tsv"))

guidesOutWindows <- guidesKeep |>
  select(accession, start, target, spacer, strand, GC_content, A_content) |>
  rename(segment=accession)

write_tsv(guidesOutWindows, paste0(inputPrefix, "_guides_windows_wHumanAlignZero_kept", cutoff, ".tsv"))
write_tsv(guidesOutWindows, "windows.txt")

guidesOutSpacer <- guidesOutWindows |>
  select(spacer)

write_tsv(guidesOutSpacer, paste0(inputPrefix, "_guides_spacers_wHumanAlignZero_kept", cutoff, ".tsv"))
write_tsv(guidesOutSpacer, "spacers.txt", col_names = F)

print(paste0("Getting fasta of ", nrow(guidesOutWindows)))
guidesOutTargets <- c()
for (i in 1:nrow(guidesOutWindows)) {
  result1 <- paste0(">", inputPrefix, "_", guidesOutWindows$segment[i], "_", guidesOutWindows$start[i])
  result2 <- toupper(guidesOutWindows$target[i])
  guidesOutTargets <- c(guidesOutTargets, result1, result2)
  #print(paste(i, "of", nrow(guidesOutWindows)))
}

#print("Writing fasta")
writeLines(guidesOutTargets, paste0(inputPrefix, "_guides_targets_wHumanAlignZero_kept", cutoff, ".fa"))
writeLines(guidesOutTargets, "targets.fa")
print("Done")



