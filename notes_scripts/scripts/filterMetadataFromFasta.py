__author__ = 'Quinn'

#import os
import argparse
import pandas as pd
from Bio import SeqIO

################################################################
# This script splits a fasta file by subtype based on matching accession 
# to pre-processed metadata file, where a column has been added for subytpeFolder
#
#Input: metadata file for species AND input fasta sequence file
#Output: seperate fasta files for the subtypes specified in subtypeFolder of metadata file 
#
################################################################

parser = argparse.ArgumentParser(description="Split fasta file into subtype fastas based on metadata table")
parser.add_argument('--metadata', help="metadata, required", required=True)
args = parser.parse_args()

metaFile = args.metadata
metaPrefix = metaFile.split(".")[0]
#load genome dataframe and generate genomeID-to-segment lookup dictionary
metaDF = pd.read_csv(metaFile, sep="\t")
metaDict = metaDF.set_index('accession')['subtypeFolder'].to_dict()
subtypeFullIDs = metaDF['subtypeFolder'].unique()

accessionList = []
for subtypeID in subtypeFullIDs:
    subFasta = subtypeID + ".fasta"
    print(subtypeID)
    for seq_record in SeqIO.parse(subFasta, "fasta"):
        seqID = seq_record.id
        accessionList.append(seqID)
    print(len(accessionList))

filterMeta = metaDF[metaDF['accession'].isin(accessionList)]
filterMeta.to_csv(metaPrefix + "_inFasta.tsv", sep="\t", index=False)

strainCount = len(set(filterMeta['strain']))
subtypeCountDict = {}
for subtypeID in subtypeFullIDs:
    subtypeCountDict[subtypeID] = {}
    count = len(set(filterMeta[filterMeta['subtypeFolder']==subtypeID]['strain']))
    subtypeCountDict[subtypeID]['count'] = count
    subtypeCountDict[subtypeID]['prop'] = count/strainCount

subtypeCountDF = pd.DataFrame(data=subtypeCountDict).T
subtypeCountDF.to_csv(metaPrefix + "_inFasta_counts.tsv", sep="\t")
