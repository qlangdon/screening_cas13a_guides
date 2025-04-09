library(tidyverse)
library(viridis)

setwd("~/Dropbox (Gladstone)/screening_cas13a_guides/RSV_datasets/")

rsvArefStrain <- "LR699737"
rsvAmeta <- read_tsv("RSV-A_nextclade_metadata_20230206.tsv")

dim(rsvAmeta)
length(unique(rsvAmeta$strain))
length(unique(rsvAmeta$clade))

rsvAmeta <- rsvAmeta |>
  separate(genbank_accession_rev, c("genbankID", "genbankV"), remove=F)
unique(rsvAmeta$genbankV)

rsvAout <- rsvAmeta |>
  mutate(length = alignmentEnd - alignmentStart) |>
  mutate(spp = "RSV-A") |>
  mutate(segment = 1) |>
  select(accession, length, alignmentStart, alignmentEnd, genbank_accession_rev, 
         spp, segment, clade, strain, country, host, date) |>
  rename(genbank_version = genbank_accession_rev) |>
  rename(subtype = clade)

#write_tsv(rsvAout, "RSV-A_sequence_description_parsed.txt")

ggplot(rsvAout) + geom_bar(aes(x=subtype, fill=subtype)) + theme_bw()

cladeCountA <- rsvAout %>%
  group_by(subtype) %>%
  summarise(subCount=n()) %>%
  arrange(desc(subCount)) 

cladeUniCountA <- rsvAout |>
  group_by(subtype) |>
  summarise(uniStrainCount = n_distinct(strain))

cladeCountA <- left_join(cladeCountA, cladeUniCountA)

cladeCountA <- cladeCountA |>
  mutate(totalCount=nrow(rsvAmeta)) |>
  mutate(propTot=subCount/totalCount) |>
  mutate(propUni = uniStrainCount/subCount)

cutoffA <- nrow(rsvAmeta)/length(unique(rsvAmeta$clade))

cladeCountA <- cladeCountA |>
  mutate(topSubtype = ifelse(subCount>cutoffA, T, F))

head(cladeCountA)
cladeCountA |>
  filter(topSubtype == F) |>
  select(subCount) |>
  sum()

topSubsA <- cladeCountA$subtype[which(cladeCountA$topSubtype==T)]

rsvAout <- rsvAout |>
  mutate(topSubtype = ifelse(subtype %in% topSubsA, T, F))

rsvAout <- rsvAout |>
  mutate(subtypeFolder = ifelse(subtype %in% topSubsA, paste0("subtype_", subtype), "subtype_allOthers"))

write_tsv(rsvAout, "RSV-A_strains_metadata.txt")
write_tsv(cladeCountA, "RSV-A_subtype_strain_counts.txt")

refLine <- rsvAout |> filter(accession==rsvArefStrain)
refStart <- refLine$alignmentStart
refEnd <- refLine$alignmentEnd

genomeAdf <- tibble()
for (i in (refStart:refEnd)) {
  winEnd <- i + 19
  if (winEnd > refEnd) {winEnd <- refEnd} 
  strainMatch <- rsvAout |>
    filter(alignmentStart<=i & alignmentEnd>=winEnd)
  tmpLine <- tibble(winStart=i, winEnd=winEnd, 
                    strainCount=nrow(strainMatch), strainUni=length(unique(strainMatch$strain)))
  tmpHeader <- c("winStart", "winEnd", "strainCount", "strainUni")
  for (clade in topSubsA) {
    cladeWinCount <- strainMatch |> filter(subtype==clade) |> summarise(n_distinct(strain))
    tmpLine <- cbind(tmpLine, cladeCount=cladeWinCount)
    tmpHeader <- c(tmpHeader, clade)
    names(tmpLine) <- tmpHeader
  }
  otherCladeCount <- strainMatch |> filter(! subtype %in% topSubsA) |> summarise(n_distinct(strain))
  tmpLine <- cbind(tmpLine, cladeCount=otherCladeCount)
  tmpHeader <- c(tmpHeader, "allOthers")
  names(tmpLine) <- tmpHeader
  genomeAdf <- rbind(genomeAdf, tmpLine)
}

genomeAdf$propTotal <- genomeAdf$strainCount/nrow(rsvAout)
genomeAdf$propUniq <- genomeAdf$strainCount/length(unique(rsvAout$strain))
dfHeader <- colnames(genomeAdf)

for (clade in topSubsA) {
  cladeProp <- genomeAdf[,which(dfHeader==clade)]/cladeCountA$uniStrainCount[which(cladeCountA$subtype==clade)]
  genomeAdf$propClade <- cladeProp
  dfHeader <- c(dfHeader, paste0("prop_", clade))
  names(genomeAdf) <- dfHeader
}

allOtherCount <- cladeCountA |>
  filter(topSubtype == F) |>
  select(uniStrainCount) |>
  sum()

genomeAdf$prop_allOthers <- genomeAdf$allOthers/allOtherCount
write_tsv(genomeAdf, "RSV-A_strainCountPerWindow.tsv")

ggplot(genomeAdf, aes(x=winStart)) + 
  geom_line(aes(y=propTotal), color="grey50") + 
  geom_line(aes(y=propUniq), color="black") + 
  geom_line(aes(y=prop_allOthers), color="red") + 
  theme_bw() 


#And the same with RSV-b
rsvBmeta <- read_tsv("RSV-B_nextclade_metadata_20230206.tsv")
dim(rsvBmeta)
length(unique(rsvBmeta$strain))
length(unique(rsvBmeta$clade))

rsvBmeta <- rsvBmeta |>
  separate(genbank_accession_rev, c("genbankID", "genbankV"), remove=F)
unique(rsvBmeta$genbankV)

rsvBout <- rsvBmeta |>
  mutate(length = alignmentEnd - alignmentStart) |>
  mutate(spp = "RSV-B") |>
  mutate(segment = 1) |>
  select(accession, length, alignmentStart, alignmentEnd, genbank_accession_rev, 
         spp, segment, clade, strain, country, host, date) |>
  rename(genbank_version = genbank_accession_rev) |>
  rename(subtype = clade)

#write_tsv(rsvBout, "RSV-B_metadata_parsed.txt")

ggplot(rsvBout) + geom_bar(aes(x=subtype, fill=subtype)) + theme_bw()

cladeCountB <- rsvBout %>%
  group_by(subtype) %>%
  summarise(subCount=n()) %>%
  arrange(desc(subCount)) 

cladeUniCountB <- rsvBout |>
  group_by(subtype) |>
  summarise(uniStrainCount = n_distinct(strain))

cladeCountB <- left_join(cladeCountB, cladeUniCountB)

cladeCountB <- cladeCountB |>
  mutate(totalCount=nrow(rsvBmeta)) |>
  mutate(propTot=subCount/totalCount) |>
  mutate(propUni = uniStrainCount/subCount)

cutoffB <- nrow(rsvBmeta)/length(unique(rsvBmeta$clade))

cladeCountB <- cladeCountB |>
  mutate(topSubtype = ifelse(subCount>cutoffB, T, F))

head(cladeCountB)
cladeCountB |>
  filter(topSubtype == F) |>
  select(subCount) |>
  sum()

topSubsB <- cladeCountB$subtype[which(cladeCountB$topSubtype==T)]

rsvBout <- rsvBout |>
  mutate(topSubtype = ifelse(subtype %in% topSubsB, T, F))

rsvBout <- rsvBout |>
  mutate(subtypeFolder = ifelse(subtype %in% topSubsB, paste0("subtype_", subtype), "subtype_allOthers"))

write_tsv(rsvBout, "RSV-B_strains_metadata.txt")
write_tsv(cladeCountB, "RSV-B_subtype_strain_counts.txt")

ggplot(rsvBout, aes(x=subtypeFolder, y=length, color=subtypeFolder)) + geom_violin() + theme_bw()

rsvBrefStrain <- "OP975389"
refBline <- rsvBout |> filter(accession==rsvBrefStrain)
refBstart <- refBline$alignmentStart
refBend <- refBline$alignmentEnd

genomeBdf <- tibble()
for (i in (refBstart:refBend)) {
  winEnd <- i + 19
  if (winEnd > refBend) {winEnd <- refBend} 
  strainMatch <- rsvBout |>
    filter(alignmentStart<=i & alignmentEnd>=winEnd)
  tmpLine <- tibble(winStart=i, winEnd=winEnd, 
                    strainCount=nrow(strainMatch), strainUni=length(unique(strainMatch$strain)))
  tmpHeader <- c("winStart", "winEnd", "strainCount", "strainUni")
  for (clade in topSubsB) {
    cladeWinCount <- strainMatch |> filter(subtype==clade) |> summarise(n_distinct(strain))
    tmpLine <- cbind(tmpLine, cladeCount=cladeWinCount)
    tmpHeader <- c(tmpHeader, clade)
    names(tmpLine) <- tmpHeader
  }
  otherCladeCount <- strainMatch |> filter(! subtype %in% topSubsB) |> summarise(n_distinct(strain))
  tmpLine <- cbind(tmpLine, cladeCount=otherCladeCount)
  tmpHeader <- c(tmpHeader, "allOthers")
  names(tmpLine) <- tmpHeader
  genomeBdf <- rbind(genomeBdf, tmpLine)
}

genomeBdf$propTotal <- genomeBdf$strainCount/nrow(rsvBout)
genomeBdf$propUniq <- genomeBdf$strainCount/length(unique(rsvBout$strain))
dfHeader <- colnames(genomeBdf)

for (clade in topSubsB) {
  cladeProp <- genomeBdf[,which(dfHeader==clade)]/cladeCountB$uniStrainCount[which(cladeCountB$subtype==clade)]
  genomeBdf$propClade <- cladeProp
  dfHeader <- c(dfHeader, paste0("prop_", clade))
  names(genomeBdf) <- dfHeader
}

allOtherCount <- cladeCountB |>
  filter(topSubtype == F) |>
  select(uniStrainCount) |>
  sum()

genomeBdf$prop_allOthers <- genomeBdf$allOthers/allOtherCount
write_tsv(genomeBdf, "RSV-B_strainCountPerWindow.tsv")

ggplot(genomeBdf, aes(x=winStart)) + 
  geom_line(aes(y=propTotal), color="grey50") + 
  geom_line(aes(y=propUniq), color="black") + 
  geom_line(aes(y=prop_B.D.4.1), color="blue") + 
  theme_bw() 

ggplot(genomeBdf, aes(x=winStart)) + 
  geom_line(aes(y=strainCount), color="grey50") + 
  geom_line(aes(y=strainUni), color="black") + 
  geom_line(aes(y=B.D.4.1), color="blue") + 
  theme_bw() 

guidesOutBfull <- read_tsv("RSV-B_all_guides_hit_summary.tsv")
nrow(guidesOutBfull |> filter(prop_win_strains_hit>1))
nrow(guidesOutBfull |> filter(prop_win_total_hit>1))

guidesOutB <- guidesOutBfull |>
  filter(refStart<refBend-30)
summary(guidesOutB)

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=total_hit)) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=total_hit)) + 
  geom_line(data=genomeBdf, aes(x=winStart, y=strainCount), color="blue") +
  theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit, color=strains_hit)) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_total_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_total_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_strains_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D.4.1.1, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D.4.1.1, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D.4.1, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D.4.1, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D.4, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D.4, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_allOthers, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_allOthers, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

nrow(guidesOutB |> filter(prop_win_B.D.4.1.1>1))
summary(guidesOutB |> filter(prop_win_B.D.4.1.1>1))

nrow(guidesOutB |> filter(prop_win_B.D >1))
summary(guidesOutB |> filter(prop_win_B.D >1))

nrow(guidesOutB |> filter(prop_win_B.D >1))
summary(guidesOutB |> filter(prop_win_B.D >1))

guidesOutB$filterColor <- "dodgerblue"
guidesOutB$filterColor[which(guidesOutB$prop_win_B.D.4.1.1>1 | 
                               guidesOutB$prop_win_B.D>1 | 
                               guidesOutB$prop_win_B.D.4.1>1 | 
                               guidesOutB$prop_win_B.D.4>1 | 
                               guidesOutB$prop_win_allOthers>1)] <- "grey20"
guidesOutB$filterColor[which(guidesOutB$strains_hit<500)] <- "sienna"

dim(guidesOutB)
dim(guidesOutB |> filter(filterColor=="dodgerblue"))

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D.4.1.1), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D.4.1.1), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D.4.1), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D.4.1), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_B.D.4), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_B.D.4), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_allOthers), color=guidesOutB$filterColor) + theme_bw()

ggplot(guidesOutB, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_allOthers), color=guidesOutB$filterColor) + theme_bw()

guidesOutBkeep <- guidesOutB |>
  filter(filterColor=="dodgerblue") |>
  select(-filterColor)

write_tsv(guidesOutBkeep, "RSV-B_all_guides_filtered.tsv")

guidesOutBwindows <- guidesOutBkeep |>
  select(accession, refStart, target, spacer, strand, GC_content, A_content) |>
  rename(segment=accession) |>
  rename(start = refStart)

write_tsv(guidesOutBwindows, "RSV-B_all_guides_windows.txt")

guidesOutBspacer <- guidesOutBwindows |>
  select(spacer)
write_tsv(guidesOutBspacer, "RSV-B_all_guides_spacers.txt", col_names = F)

guidesOutBtargets <- c()
for (i in 1:nrow(guidesOutBwindows)) {
  result1 <- paste0(">", guidesOutBwindows$segment[i], "_", guidesOutBwindows$start[i])
  result2 <- guidesOutBwindows$target[i]
  guidesOutBtargets <- c(guidesOutBtargets, result1, result2)
}

writeLines(guidesOutBtargets, "RSV-B_all_guides_targets.fa")


##Checking the output for RSV-A
guidesOutAfull <- read_tsv("RSV-A_all_guides_hit_summary.tsv")
summary(guidesOutAfull)
nrow(guidesOutAfull |> filter(prop_win_strains_hit>1))
nrow(guidesOutAfull |> filter(prop_win_total_hit>1))

ggplot(guidesOutAfull, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit)) + theme_bw()

guidesOutA <- guidesOutAfull |>
  filter(prop_win_strains_hit<1)
summary(guidesOutA)

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=total_hit)) + theme_bw()

#ggplot(guidesOutA, aes(x=refStart)) + 
#  geom_point(aes(y=total_hit)) + 
#  geom_line(data=genomeAdf, aes(x=winStart, y=strainCount), color="blue") +
#  theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit, color=strains_hit)) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_total_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_total_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_strains_hit, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.3.1, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3.1, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D.1, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D.1, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D.3, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D.3, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.3.1.1, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3.1.1, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.3, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D.2.2, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D.2.2, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_allOthers, color=strains_hit)) + 
  scale_color_viridis() + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_allOthers, color=strains_hit)) + 
  scale_color_viridis(direction=-1) + 
  geom_hline(yintercept=1, color="grey50") + theme_bw()

nrow(guidesOutA |> filter(prop_win_A.D>1))

nrow(guidesOutA |> filter(prop_win_A.3.1 >1))
summary(guidesOutA |> filter(prop_win_A.3.1 >1))

guidesOutA$filterColor <- "dodgerblue"
guidesOutA$filterColor[which(guidesOutA$prop_win_A.D>1 | 
                               guidesOutA$prop_win_A.3.1>1 | 
                               guidesOutA$prop_win_A.D.1>1 | 
                               guidesOutA$prop_win_A.3.1.1>1 | 
                               guidesOutA$prop_win_A.D.3>1 |  
                               guidesOutA$prop_win_A.3>1 | 
                               guidesOutA$prop_win_A.D.2.2>1 | 
                               guidesOutA$prop_win_allOthers>1)] <- "grey20"
guidesOutA$filterColor[which(guidesOutA$strains_hit<500)] <- "sienna"

dim(guidesOutA)
dim(guidesOutA |> filter(filterColor=="dodgerblue"))

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_strains_hit), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.3.1), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3.1), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D.1), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3.1.1), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.3.1.1), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3.1.1), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D.3), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D.3), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.3), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.3), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_A.D.2.2), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_A.D.2.2), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=refStart)) + 
  geom_point(aes(y=prop_win_allOthers), color=guidesOutA$filterColor) + theme_bw()

ggplot(guidesOutA, aes(x=prop_win_strains_hit)) + 
  geom_point(aes(y=prop_win_allOthers), color=guidesOutA$filterColor) + theme_bw()

guidesOutAkeep <- guidesOutA |>
  filter(filterColor=="dodgerblue") |>
  select(-filterColor)

write_tsv(guidesOutAkeep, "RSV-A_all_guides_filtered.tsv")

guidesOutAwindows <- guidesOutAkeep |>
  select(accession, refStart, target, spacer, strand, GC_content, A_content) |>
  rename(segment=accession) |>
  rename(start = refStart)

write_tsv(guidesOutAwindows, "RSV-A_all_guides_windows.txt")

guidesOutAspacer <- guidesOutAwindows |>
  select(spacer)
write_tsv(guidesOutAspacer, "RSV-A_all_guides_spacers.txt", col_names = F)

guidesOutAtargets <- c()
for (i in 1:nrow(guidesOutAwindows)) {
  result1 <- paste0(">", guidesOutAwindows$segment[i], "_", guidesOutAwindows$start[i])
  result2 <- guidesOutAwindows$target[i]
  guidesOutAtargets <- c(guidesOutAtargets, result1, result2)
}

writeLines(guidesOutAtargets, "RSV-A_all_guides_targets.fa")

