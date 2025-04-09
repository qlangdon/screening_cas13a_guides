__author__ = 'Quinn'

#import os
import argparse
#from re import split
import pandas as pd
import numpy as np


################################################################
# This script goes through tiled 20bp widows for all strains of a species to make a dataframe with counts of that 
# target across all genomes and matches to the metadata.
#
#Input: metadata file for species, path to the windowed genomes in chunks
#Output: full dataframe of all unique targets and the number of times they're found 
#
################################################################

parser = argparse.ArgumentParser(description="Go through windowed genomes and count target frequency")
parser.add_argument('--metadata', help="metadata, required", required=True)
#parser.add_argument('--pathPrefix', required=True,
#                    help='The prefix of the folders with the chunked genomes')
#parser.add_argument('--out', help="Output prefix, required", required=True)
args = parser.parse_args()

#folderPathPrefix=args.pathPrefix
#folderPathPrefix='tiled_genomes/tiled_subtype_'
#metaDF = pd.read_csv('IAV_complete_strains_metadata.txt', sep="\t")

#load genome dataframe and generate genomeID-to-segment lookup dictionary
metaDF = pd.read_csv(args.metadata, sep="\t")
metaDict = metaDF.set_index('accession').T.to_dict('list')
topClades = list(set(metaDF[metaDF['topSubtype']]['subtype']))
topClades.append("allOthers")
topCladesDict = {}
for clade in topClades:
     count = len(set(metaDF[metaDF['subtype']==clade]['strain']))
     topCladesDict[clade] = count

topCladesDict['allOthers'] = len(set(metaDF[metaDF['subtypeFolder']=='subtype_allOthers']['strain']))

topCladesDictSorted = dict(sorted(topCladesDict.items(), key=lambda x:x[1], reverse=True))

segmentListFull = list(set(metaDF['segment']))
segmentList = []
for segID in segmentListFull:
    if not np.isnan(segID):
        #print(segID)
        segmentList.append(segID)

segmentList.append(float('nan'))

segmentDict = {}
for seg in segmentList:
    if not np.isnan(seg):
         count = len(set(metaDF[metaDF['segment']==seg]['strain']))
         segmentDict[seg] = count
    else:
        naCount = sum(np.isnan(metaDF['segment']))
        segmentDict[seg] = naCount

allCladesDict = {}
for curClade in topClades:
    print(curClade)
    sumCountDF = pd.read_csv("tiled_subtype_" + curClade + "_guides_hit_summary.tsv", sep="\t")
    print(len(sumCountDF))
    sumCountDict = sumCountDF.set_index('target').T.to_dict('list')
    for targetSeq in sumCountDF['target']:
        curCounts = sumCountDict[targetSeq]
        curAcc = curCounts[0]
        curStart = curCounts[1]
        curSeg = curCounts[2]
        curTot = curCounts[3]
        curStrain = curCounts[4]
        curSpacer = curCounts[5]
        curStrand = curCounts[6]
        curGC = curCounts[7]
        curA = curCounts[8]
        curCladeCount = curCounts[9]
        cur1 = curCounts[10]
        cur2 = curCounts[11]
        cur3 = curCounts[12]
        cur4 = curCounts[13]
        cur5 = curCounts[14]
        cur6 = curCounts[15]
        cur7 = curCounts[16]
        cur8 = curCounts[17]
        if targetSeq in allCladesDict.keys():
            preAcc = allCladesDict[targetSeq]['accession']
            if curAcc < preAcc:
                allCladesDict[targetSeq]['accession'] = curAcc
                allCladesDict[targetSeq]['start'] = curStart
                allCladesDict[targetSeq]['segment'] = curSeg
            allCladesDict[targetSeq]['total_hit'] += curTot
            allCladesDict[targetSeq]['strains_hit'] += curStrain
            allCladesDict[targetSeq][1.0] += cur1
            allCladesDict[targetSeq][2.0] += cur2
            allCladesDict[targetSeq][3.0] += cur3
            allCladesDict[targetSeq][4.0] += cur4
            allCladesDict[targetSeq][5.0] += cur5
            allCladesDict[targetSeq][6.0] += cur6
            allCladesDict[targetSeq][7.0] += cur7
            allCladesDict[targetSeq][8.0] += cur8
            for cladeID in topClades:
                if cladeID == curClade:
                    countToAdd = curCladeCount
                else:
                    countToAdd = 0
                allCladesDict[targetSeq][cladeID] += countToAdd 
        else:
            targetDict = {'accession':curAcc, 'start':curStart, 'target':targetSeq, 'segment':curSeg, 'total_hit':curTot, 
                                    'strains_hit':curStrain, 'spacer':curSpacer, 'strand':curStrand, 'GC_content':curGC, 'A_content':curA, 
                                    1.0:cur1, 2.0:cur2, 3.0:cur3, 4.0:cur4, 5.0:cur5, 6.0:cur6, 7.0:cur7, 8.0:cur8 }
            for cladeID in topClades:
                if cladeID == curClade:
                    countToAdd = curCladeCount
                else:
                    countToAdd = 0
                targetDict[cladeID] = countToAdd 
            allCladesDict[targetSeq] = targetDict

allCladesDF = pd.DataFrame(data=allCladesDict).T
allCladesDFsorted = allCladesDF.sort_values(by=['total_hit'], ascending = False)
allCladesDFsorted['prop_strains_hit'] = allCladesDFsorted['strains_hit']/len(set(metaDF['strain']))
for clade in topCladesDict.keys():
    allCladesDFsorted['prop_'+clade] = allCladesDFsorted[clade]/topCladesDict[clade]

allCladesDFsortedReorder = allCladesDFsorted.iloc[:,[0,1,2,3,4,5,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37]]
allCladesDFsortedReorder.to_csv("IAV_all_guides_hit_summary.tsv", sep="\t", index=False)


        






   
