__author__ = 'Quinn'

import Bio
from Bio import SeqIO
from Bio.Seq import MutableSeq
import sys, re, argparse

parser = argparse.ArgumentParser(description="Pick sequences wanted and rename from FASTA (with strain key)")
parser.add_argument('--out', help="Output suffix")
parser.add_argument('--fasta', help="fasta to use, required", required=True)
parser.add_argument('--strainKey', help="Key to strains wanted and how to rename them", required=True)
args = parser.parse_args()

#print(len(sys.argv))
fastaName = args.fasta
if args.out:
    outID = args.out
else:
    outID = "seqsPicked_renamed"

FASTAoutputFile = fastaName + "_" + outID + ".fa"
keyOutputFile = fastaName + "_" + outID + "_outputKey.txt"

strainKeyFile = args.strainKey
renameDict = {}
with open(strainKeyFile) as strainsRename:
    next(strainsRename)
    for line in strainsRename:
        currentLine = line.strip('\n').split('\t')
        accession = currentLine[0]
        newName = currentLine[1]
        renameDict[accession] = newName

keyOutput = open(keyOutputFile, 'w')
keyOutput.write("GenBankAccession\tMetadataID\tfnaDescription\tseqLen\n")
FASTAoutput = open(FASTAoutputFile, 'w')
fasta = open(fastaName, 'r')
for seq_record in SeqIO.parse(fasta, "fasta"):
    seqID = seq_record.id
    if seqID in renameDict.keys():
        newName = renameDict[seqID]
        FASTAoutput.write(">" + seqID + " " + newName + "\n" + str(seq_record.seq) + "\n")
        keyOutput.write(seqID + "\t" + newName + "\t" + seq_record.description + "\t" + str(len(seq_record.seq)) + "\n")

keyOutput.close()
FASTAoutput.close()
