__author__ = 'Quinn'

#import os
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO

################################################################
# This script splits a fasta file by subtype based on matching accession or chosen column 
# to pre-processed metadata file, where a column has been added for subytpeFolder
#
#Input: metadata file for species AND input fasta sequence file, 
# optionally you can choose the columns that contain the unique ids and that you want to group by
#Output: seperate fasta files for the subtypes specified in subtypeFolder of metadata file 
#
################################################################

parser = argparse.ArgumentParser(description="Split fasta file into subtype fastas based on metadata table")
parser.add_argument('--metadata', help="metadata file, required", required=True)
parser.add_argument('--fasta', required=True,
                    help='The fasta file to split')
parser.add_argument('--idColumn', help="column in the metadata with the unique identifer that will match metdata to fasta", 
                    default='accession')
parser.add_argument('--subtypeColumn', help="column in the metadata with the information to split on", 
                    default='subtypeFolder')
#parser.add_argument('--filter', help="to omit writting some sequences, True or False", default=False)
args = parser.parse_args()

fastaName = args.fasta
metaFile = args.metadata
idCol = args.idColumn
subCol = args.subtypeColumn

def splitByGroup(metaFile, idCol, subCol, fastaName):
    """get a list of folders that contains the prefix
    Args:
        the intput metadata file and fasta. And the column to split by and the column witht the unique identifier
    Returns:
        list: a list of groups that the data is split into
    """
    #load genome dataframe and generate genomeID-to-segment lookup dictionary
    metaDF = pd.read_csv(metaFile, sep="\t")
    folderIDs = metaDF[subCol].unique()
    ##The id names having / in them breaks later steps so fixing the names here
    if metaDF[idCol].str.contains("/").any():
        metaDF.insert(0, 'renamedIDs', metaDF[idCol].str.replace('/', '_'))
        idCol = 'renamedIDs'
    ##Get a count of what's in the metadata file by the choosen groups
    metaCount = metaDF.groupby(subCol).count()[idCol].rename('idCount').to_frame()
    metaCount['propTotal'] = metaCount['idCount']/len(metaDF)
    if 'strain' in metaDF.columns:
        metaCount = metaCount.merge(metaDF.groupby(subCol)['strain'].nunique().rename('uniqueStrains'), left_index=True, right_index=True)
        metaCount['propStrains'] = metaCount['uniqueStrains']/metaDF['strain'].nunique()
    ##Don't want to create a ton of fastas so want to figure out which groups have the most then lump the others together
    ##Currently the cutoff is #groups / #total entries
    fullList = metaDF[subCol].unique()
    if len(folderIDs) > 20:
        cutoff = len(metaDF)/len(fullList)
        metaCount['topGroup'] = metaCount['idCount'] >= cutoff
    else:
        metaCount['topGroup'] = True
    #Adding new info about if the groups are getting reshuffled
    topGroups = metaCount.index[metaCount['topGroup']==True].tolist()
    metaCount['assignedGroup'] = np.where(metaCount['topGroup']==True, metaCount.index, 'allOthers')
    metaDF['topGroup'] = metaDF[subCol].isin(topGroups)
    metaDF['assignedGroup'] = np.where(metaDF['topGroup']==True, metaDF[subCol], 'allOthers')
    #Writing this summary out here
    metaCountOutBase = metaFile.split(".")[0]
    metaCount.to_csv(metaCountOutBase + "_strainCounts.tsv", sep="\t")
    ##Now making a dict that will store the ids that need to be with each group
    metaDict = metaDF.set_index(idCol)['assignedGroup'].to_dict()
    fastaDicts = {}
    for id in folderIDs:
        fastaDicts[id] = ""
    accInFasta = []
    fasta = open(fastaName, 'r')
    for seq_record in SeqIO.parse(fasta, "fasta"):
        seqID = str(seq_record.id)
        if "/" in seqID:
            seqID = seqID.replace('/', '_')
        accInFasta.append(str(seqID))
        if seqID in metaDict:
            seqFolder = metaDict[seqID]
            fastaDicts[seqFolder] = fastaDicts[seqFolder] + ">" + seqID + "\n" + seq_record.seq + "\n"
    fasta.close() 
    fastaBase = fastaName.split(".f")[0]
    for id in folderIDs:
        fastaOut = open(fastaBase+"_assignedGroup_"+id+".fasta", 'w')
        fastaOut.write(str(fastaDicts[id]))
        fastaOut.close()
    ##Now check which ids were actually in the fasta file
    metaDF['inFASTA'] = metaDF[idCol].isin(accInFasta)
    metaCount = metaCount.merge(metaDF[metaDF['inFASTA'] == True].groupby(subCol).count()[idCol].rename('idInFASTA'), left_index=True, right_index=True)
    metaCount['propInFASTA'] = metaCount['idInFASTA']/len(accInFasta)
    if 'strain' in metaDF.columns:
        metaCount = metaCount.merge(metaDF[metaDF['inFASTA'] == True].groupby(subCol)['strain'].nunique().rename('strainsInFASTA'), left_index=True, right_index=True)
        metaCount['propStrainsInFASTA'] = metaCount['strainsInFASTA']/metaDF[metaDF['inFASTA'] == True]['strain'].nunique()
    metaCount.to_csv(metaCountOutBase + "_strainCounts.tsv", sep="\t")
    metaDF.to_csv(metaCountOutBase + "_postParsing.tsv", sep="\t")
    return(topGroups)


groupsList = splitByGroup(metaFile, idCol, subCol, fastaName)

groupListFile = open("assignedGroupsList.txt", 'w')
for group in groupsList:
    groupListFile.write(group+"\n")

groupListFile.close()


