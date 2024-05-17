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
parser.add_argument('--fasta', required=True,
                    help='The fasta file to split')
#parser.add_argument('--filter', help="to omit writting some sequences, True or False", default=False)
args = parser.parse_args()

fastaName = args.fasta
metaFile = args.metadata
#load genome dataframe and generate genomeID-to-segment lookup dictionary
metaDF = pd.read_csv(metaFile, sep="\t")
metaDict = metaDF.set_index('accession')['subtypeFolder'].to_dict()
folderIDs = metaDF['subtypeFolder'].unique()

fastaDicts = {}
for id in folderIDs:
    fastaDicts[id] = ""

fasta = open(fastaName, 'r')
for seq_record in SeqIO.parse(fasta, "fasta"):
    seqID = seq_record.id
    if seqID in metaDict:
        seqFolder = metaDict[seqID]
        fastaDicts[seqFolder] = fastaDicts[seqFolder] + ">" + seqID + "\n" + seq_record.seq + "\n"

fasta.close() 

for id in folderIDs:
    fastaOut = open(id+".fasta", 'w')
    fastaOut.write(str(fastaDicts[id]))
    fastaOut.close()

