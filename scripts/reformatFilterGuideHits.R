#! /usr/bin/Rscript
options(stringAsFactors=FALSE)
library(tidyverse)

args <- commandArgs(TRUE)
inputPrefix <- args[1]
cutoff <- as.numeric(args[2])

guideFile <- paste0(inputPrefix, "_guides_hit_summary.tsv")
guideDF <- read_tsv(guideFile)
print(paste("Input DF", nrow(guideDF)))

if ("prop_subtypes_hit" %in% colnames(guideDF)) {
  guideDF <- guideDF |> rename(p_subtypes_hit =prop_subtypes_hit)
}
if ("refStart" %in% colnames(guideDF)) {
  guideDF <- guideDF |> rename(start = refStart)
}

#guideDF <- guideDF |> rename(GC_content = GC)
#guideDF <- guideDF |> rename(A_content = A)

guidesKeep <- guideDF |>
  filter(if_any(starts_with("prop"),  ~ .x > cutoff))

print(paste("Filtered DF", nrow(guidesKeep)))

print(paste("Kept", nrow(guidesKeep), "of input", nrow(guideDF), nrow(guidesKeep)/nrow(guideDF)))

write_tsv(guidesKeep, paste0(inputPrefix, "_guides_filtered", cutoff, ".tsv"))

guidesOutWindows <- guidesKeep |>
  select(accession, start, target, spacer, strand, GC_content, A_content) |>
  rename(segment=accession)

write_tsv(guidesOutWindows, paste0(inputPrefix, "_guides_windows_filtered", cutoff, ".txt"))
write_tsv(guidesOutWindows, "windows.txt")

guidesOutSpacer <- guidesOutWindows |>
  select(spacer)

write_tsv(guidesOutSpacer, paste0(inputPrefix, "_guides_spacers_filtered", cutoff, ".txt"))
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
writeLines(guidesOutTargets, paste0(inputPrefix, "_guides_targets_filtered", cutoff, ".fa"))
writeLines(guidesOutTargets, "targets.fa")
print("Done")

