__author__ = 'Quinn'

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import os
import shutil
#import time

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
parser.add_argument('--fasta', required=True, help='The fasta file to split')
parser.add_argument('--idColumn', help="column in the metadata with the unique identifer that will match metdata to fasta", 
                    default='accession')
parser.add_argument('--subtypeColumn', help="column in the metadata with the information to split on", 
                    default='subtypeFolder')
parser.add_argument('--segmented', help="if true will count both unique ids in the idColumn and assume the Strain column contains the full strain name", 
                    action='store_true')
parser.add_argument('--strainColumn', help="column in the metadata with the identifer that corresponds to strain or isolate", 
                    default='strain')
parser.add_argument("-w", "--window", type=int, default=20,
                    help="window size")
parser.add_argument("-e", "--enzyme", type=str, help="Cas enzyme type", default='Cas13a')  # "Cas13a" or "Cas12"
parser.add_argument("-g", "--genome_strand", type=str, default="+",
                    help="strandedness of viral genome")
parser.add_argument("-s", "--strand", type=str, default="+",
                    help="strand to target")
parser.add_argument("-o", "--out", type=str, default="windowed_genomes",
                    help="output directory")

#parser.add_argument('--filter', help="to omit writting some sequences, True or False", default=False)
args = parser.parse_args()

fastaName = args.fasta
metaFile = args.metadata
idCol = args.idColumn
subCol = args.subtypeColumn
strainColumn = args.strainColumn

#fastaName = 'h5n1_cattle-outbreak_concat_segKeyed.fasta'
#metaFile = 'nextstrain_avian-flu_h5n1-cattle-outbreak_genome_metadata_longForPipeline_original.txt'
#idCol = 'accession'
#subCol = 'host'

def splitByGroup(metaFile, idCol, subCol, fastaName):
    """get a list of folders that contains the prefix
    Args:
        the intput metadata file and fasta. And the column to split by and the column witht the unique identifier
    Returns:
        list: a list of groups that the data is split into
    """
    #load genome dataframe and generate genomeID-to-segment lookup dictionary
    metaDF = pd.read_csv(metaFile, sep="\t")
    ##The id names having / in them breaks later steps so fixing the names here
    if metaDF[idCol].str.contains("/").any():
        metaDF.insert(0, 'renamedIDs', metaDF[idCol].str.replace('/', '_'))
        idCol = 'renamedIDs'
    ##Get a count of what's in the metadata file by the choosen groups
    ###Filling missing Clade info
    if metaDF[subCol].isnull().sum() != 0:
        metaDF = metaDF.fillna({subCol: "null"})
        folderIDs = metaDF[subCol].unique()
    folderIDs = metaDF[subCol].unique()
    metaCount = metaDF.groupby(subCol).count()[idCol].rename('idCount').to_frame()
    metaCount['propTotal'] = metaCount['idCount']/len(metaDF)
    if args.segmented == True:
        metaCount = metaCount.merge(metaDF.groupby(subCol)[strainColumn].nunique().rename('uniqueStrains'), left_index=True, right_index=True)
        metaCount['propStrains'] = metaCount['uniqueStrains']/metaDF[strainColumn].nunique()   
    ##Don't want to create a ton of fastas so want to figure out which groups have the most then lump the others together
    if len(folderIDs) > 5:
        ##Maybe lump every group that's less than 1% of the data first
        metaCount['below1'] = metaCount['propTotal'] <= 0.01
        metaCount['newGroups'] = np.where(metaCount['below1']==False, metaCount.index, 'smallClades')
        newCount = metaCount.groupby('newGroups').sum()['idCount'].rename('newSum').to_frame()
        newCount['propTotal'] = newCount['newSum']/newCount['newSum'].sum()
        newClades = metaCount['newGroups'].unique()
        smallCladeList =  metaCount.index[metaCount['below1']==True].unique()
        metaDF['smallClade'] = metaDF[subCol].isin(smallCladeList)
        metaDF['assignedGroup'] = np.where(metaDF['smallClade']==False, metaDF[subCol], 'smallClades')
    else:
        metaCount['newGroups'] = metaCount.index
        newClades = metaDF[subCol].unique()
        metaDF['assignedGroup'] = metaDF[subCol]
    #Still want to collapse further if there's a lot of groups
    ##Currently the cutoff is #groups / #total entries
    if len(newClades) > 20:
        cutoff = newCount['newSum'].sum()/len(newCount)
        metaCount['topGroup'] = metaCount['idCount'] >= cutoff
    else:
        metaCount['topGroup'] = True
    topGroups = metaCount['newGroups'][metaCount['topGroup']==True].unique()
    metaCount['assignedGroup'] = np.where(metaCount['topGroup']==True, metaCount['newGroups'], 'allOthers')
    #Adding new info about if the groups are getting reshuffled
    metaDF['topGroup'] = metaDF['assignedGroup'].isin(topGroups)
    metaDF['assignedGroup'] = np.where(metaDF['topGroup']==True, metaDF['assignedGroup'], 'allOthers')
    #Writing this summary out here
    metaCountOutBase = metaFile.split(".")[0]
    metaCount.to_csv(metaCountOutBase + "_uniqueCounts.tsv", sep="\t")
    ##Now making a dict that will store the ids that need to be with each group
    folderIDs = metaDF['assignedGroup'].unique()
    metaDict = metaDF.set_index(idCol)['assignedGroup'].to_dict()
    fastaDicts = {}
    for id in folderIDs:
        fastaDicts[id] = ""
    accInFasta = []
    fasta = open(fastaName, 'r')
    for seq_record in SeqIO.parse(fasta, "fasta"):
        seqID = str(seq_record.id)
        seqShort = seqID.split(".")[0]
        if len(seq_record.seq) > 23:
            if "/" in seqID:
                seqID = seqID.replace('/', '_')
            if seqID in metaDict:
                seqFolder = metaDict[seqID]
                fastaDicts[seqFolder] = fastaDicts[seqFolder] + ">" + seqID + "\n" + seq_record.seq + "\n"
                accInFasta.append(str(seqID))
            elif seqShort in metaDict:
                seqFolder = metaDict[seqShort]
                fastaDicts[seqFolder] = fastaDicts[seqFolder] + ">" + seqShort + "\n" + seq_record.seq + "\n"
                accInFasta.append(str(seqShort))
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
    if args.segmented == True:
        metaCount = metaCount.merge(metaDF[metaDF['inFASTA'] == True].groupby(subCol)[strainColumn].nunique().rename('strainsInFASTA'), left_index=True, right_index=True)
        metaCount['propStrainsInFASTA'] = metaCount['strainsInFASTA']/metaDF[metaDF['inFASTA'] == True][strainColumn].nunique()
    metaCount.to_csv(metaCountOutBase + "_uniqueCounts.tsv", sep="\t")
    metaDF.to_csv(metaCountOutBase + "_postParsing.tsv", sep="\t")
    topGroups = metaCount.sort_values('idCount', ascending=True).index.tolist()
    print("Sorted " + fastaName + " into " + str(len(folderIDs)) + " groups")
    return(folderIDs)

groupsList = splitByGroup(metaFile, idCol, subCol, fastaName)

groupListFile = open("assignedGroupsList.txt", 'w')
for group in groupsList:
    groupListFile.write(group+"\n")

groupListFile.close()


##From the wrapper for breaking genomes
# make output directory
#groupTest = "otherMammal"
fastaBase = fastaName.split(".f")[0]
#inputFastaTest = fastaBase+"_assignedGroup_"+groupTest+".fasta" #args.input
outDirBase = "windowed_genomes" #args.out
outDirBase = args.out

def prepForBreak(groupID, fastaBase, outDirBase):
    inputGroupFasta = fastaBase+"_assignedGroup_"+groupID+".fasta" #args.input
    outDir = outDirBase + "/" + groupID
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    else:
        print(f"Purging output directory {outDir}")
        shutil.rmtree(outDir)
        os.makedirs(outDir)
    #read in genomes and write to individual fasta files
    print("Preparing genomes for break_genome.py...", flush=True)
    list_of_genomes = []
    genome_count = 0
    for record in SeqIO.parse(inputGroupFasta, "fasta"):
        if not os.path.exists(f"{outDir}/{record.id}"):
            os.makedirs(f"{outDir}/{record.id}")
        with open(f"{outDir}/{record.id}/{record.id}.fa", "w") as handle:
            handle.write(f">{record.description}\n{record.seq}\n")
        list_of_genomes.append(record.id)
        genome_count += 1
        print(f"{genome_count} genomes prepared", flush=True) if genome_count % 10000 == 0 else None
    print(f"Done. {genome_count} genomes prepared")
    return(list_of_genomes)

groupDict = {}
for group in groupsList:
    idList = prepForBreak(group, fastaBase, outDirBase)
    groupDict[group] = idList
    print(group + " Done")

##For breaking one genome
enzyme = args.enzyme
gStrand = args.genome_strand
strand = args.strand
windowSize = args.window

#enzyme = 'Cas13a' #args.enzyme
#gStrand = "+" #args.genome_strand
#strand = "+" #args.strand
#windowSize = 20 #args.window

#focalSample = 'A_TIGER_USA_24-038268-001_2024__ha'

def splitIndividualGenome(outDir, focalSample, enzyme, gStrand, strand, windowSize):
    focalOut = f"{outDir}/{focalSample}/" #args.out
    focalFasta = f"{outDir}/{focalSample}/{focalSample}.fa" #args.input
    genome_seq = SeqIO.to_dict(SeqIO.parse(focalFasta, "fasta"))
    windows = []
    for seq_id, seq in genome_seq.items():
        seq_str = str(seq.seq)
        num_windows = len(seq_str) - windowSize - 3
        for i in range(num_windows):
            windows.append({
                'segment': seq_id.split(" ")[0],
                'start': i + 1,
                'target': seq_str[i: i + windowSize]
            })
    df = pd.DataFrame(windows)
    num_rows, num_columns = df.shape
    if num_rows == 0:
        print(f"No windows found in genome {focalFasta}. Exiting...")
        exit(1)
    if enzyme == "Cas13a":
        if strand == gStrand:
            df['spacer'] = df['target'].apply(lambda x: Seq(x).reverse_complement().transcribe())
            df['strand'] = gStrand
        else:
            df['spacer'] = df['target'].apply(lambda x: Seq(x).transcribe())
            df['target'] = df['target'].apply(lambda x: Seq(x).reverse_complement())
            df['strand'] = strand
    else:
        df_minusStrand = df.copy()
        df_minusStrand['target'] = df['spacer'].apply(lambda x: Seq(x).transcribe())
        df_minusStrand['spacer'] = df['target'].apply(lambda x: Seq(x).transcribe())
        df_minusStrand['strand'] = '-'
        df = pd.concat([df, df_minusStrand], ignore_index=True)
    df['GC_content'] = df['spacer'].apply(lambda x: round(gc_fraction(Seq(x)),2))
    df['A_content'] = df['spacer'].apply(lambda x: x.upper().count('A')/windowSize)
    df.to_csv(f"{focalOut}/windows.txt", sep='\t', index=False)
    df['target'].to_csv(f"{focalOut}/targets.txt", index=False)
    df['spacer'].to_csv(f"{focalOut}/spacers.txt", index=False)
    with open(f"{focalOut}/targets.fa", 'w') as f:
        for idx, row in df.iterrows():
            f.write(f">{row['segment']}_{row['start']}\n")
            f.write(f"{row['target']}\n")


for groupID in groupsList:
    #start = time.time()
    outDir = outDirBase + "/" + groupID
    idList = groupDict[groupID]
    for focalSample in idList:
        splitIndividualGenome(outDir, focalSample, enzyme, gStrand, strand, windowSize)
    print(groupID + " Done")
    #end = time.time()
    #length = end - start
    #print("It took", length, "seconds!")



#start = time.time()
#end = time.time()
#length = end - start

# Show the results : this can be altered however you like
#print("It took", length, "seconds!")


