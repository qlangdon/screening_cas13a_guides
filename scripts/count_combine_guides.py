__author__ = 'Quinn'

import os
import argparse
import pandas as pd
from pathlib import Path

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
parser.add_argument('--pathPrefix', required=True,
                    help='The prefix of the folders with the chunked genomes')
#parser.add_argument('--out', help="Output prefix, required", required=True)
args = parser.parse_args()

folderPathPrefix=args.pathPrefix
#folderPathPrefix='tiled_genomes/tiled_subtype_'
metaDF = pd.read_csv(args.metadata, sep="\t")

#load genome dataframe and generate genomeID-to-segment lookup dictionary
metaDF = pd.read_csv(args.metadata, sep="\t")
#metaDict = metaDF.set_index('accession').T.to_dict('list')
metaDict = metaDF.set_index('accession').T.to_dict()
topClades = list(set(metaDF[metaDF['topSubtype']]['subtype']))
topCladesDict = {}
for clade in topClades:
     count = len(set(metaDF[metaDF['subtype']==clade]['strain']))
     topCladesDict[clade] = count

topCladesDictSorted = dict(sorted(topCladesDict.items(), key=lambda x:x[1], reverse=True))

def get_folder_list(folderPathPrefix):
    """get a list of folders that contains the prefix
    Args:
        folderPathPrefix (_type_): path to the parent folder to be scanned
    Returns:
        list: a list of folder paths
    """
    folder_list = []
    parent_dir = os.path.dirname(folderPathPrefix)
    prefix = os.path.basename(folderPathPrefix)
    for lv1_dir in os.listdir(parent_dir):
        if prefix in lv1_dir and os.path.isdir(os.path.join(parent_dir,lv1_dir)): # go through all folders that contains prefix
            folder_list.append(os.path.join(parent_dir, lv1_dir))
    return folder_list

def getChunkDict(folderBase):
    """go through each chunk folder and get a dictionary of all the guides and their hits
    Args:
        folderBase (_type_): base folder to search
    Returns:
        dictionary: a dictionary of all the guides in that base folder
    """
    pathList = os.listdir(folderBase)
    guideCountDict = {} #dictionary of guides where the key is the target seq and value is a dictionary of extra info
    for pathStep in pathList:
        curAcc = pathStep
        curMeta = metaDict[curAcc]
        windowsPath =  os.path.join(folderBase, pathStep, "windows.txt")
        windowsFile = open(windowsPath, 'r')
        windowLines = windowsFile.readlines()[1:]
        for line in windowLines:
            winLine = line.strip().split("\t")
            winAcc = winLine[0]
            winStart = winLine[1]
            winTarget = winLine[2]
            winSpacer = winLine[3]
            winStrand = winLine[4]
            winGC = winLine[5]
            winA = winLine[6]
            if winTarget in guideCountDict.keys():
                guideCountDict[winTarget]['accW'].append(winAcc)
                guideCountDict[winTarget]['start'].append(winStart)
                guideCountDict[winTarget]['target'].append(winTarget)
                guideCountDict[winTarget]['spacer'].append(winSpacer)
                guideCountDict[winTarget]['strand'].append(winStrand)
                guideCountDict[winTarget]['GC'].append(winGC)
                guideCountDict[winTarget]['A'].append(winA)
                guideCountDict[winTarget]['accM'].append(curAcc)
                for metaKey in curMeta.keys():
                    guideCountDict[winTarget][metaKey].append(curMeta[metaKey])
            else:
                winDict = {'accW':[winAcc], 'start':[winStart], 'target':[winTarget], 'spacer':[winSpacer], 'strand':[winStrand], 'GC':[winGC], 'A':[winA], 'accM':[curAcc]}
                for metaKey in curMeta.keys():
                    winDict[metaKey] = [curMeta[metaKey]]
                guideCountDict[winTarget] = winDict
        windowsFile.close()
    return guideCountDict

def sumCountDict(guideCountDict, accList):
    """go full hits dictionary and make a summary dataframe
    Args:
        guideCountDict (_type_): guide hit dictionary
        accList : list of the accessions that went into that dictionary (to properly count the metadata)
    Returns:
        dataframe: summary dataframe of all guides hits
    """
    chunkMetaDF = metaDF[metaDF['accession'].isin(accList)]
    chunkTopCladesDict = {}
    strainsCladesDict = {}
    for clade in topCladesDictSorted.keys():
        count = len(set(chunkMetaDF[chunkMetaDF['subtype']==clade]['strain']))
        chunkTopCladesDict[clade] = count
        strainsCladesDict[clade] = set(chunkMetaDF[chunkMetaDF['subtype']==clade]['strain'])
    guideChunkSumDict = {'accession':[], 'start':[], 'target':[], 'total_hit':[], 'strains_hit':[], 'subtypes_hit':[], 'spacer':[], 'strand':[], 'GC_content':[], "A_content":[]}
    for clade in topCladesDictSorted.keys():
        guideChunkSumDict[clade] = []
    for target in guideCountDict.keys():
        targetDict = guideCountDict[target]
        firstAcc = sorted(targetDict['accW'])[0]
        guideChunkSumDict['accession'].append(firstAcc)
        firstIndex = targetDict['accW'].index(firstAcc)
        guideChunkSumDict['start'].append(targetDict['start'][firstIndex])
        guideChunkSumDict['target'].append(targetDict['target'][firstIndex])
        guideChunkSumDict['spacer'].append(targetDict['spacer'][firstIndex])
        guideChunkSumDict['strand'].append(targetDict['strand'][firstIndex])
        guideChunkSumDict['A_content'].append(targetDict['A'][firstIndex])
        guideChunkSumDict['GC_content'].append(targetDict['GC'][firstIndex])
        guideChunkSumDict['total_hit'].append(len(set(targetDict['accM'])))
        guideChunkSumDict['strains_hit'].append(len(set(targetDict['strain'])))
        guideChunkSumDict['subtypes_hit'].append(len(set(targetDict['subtype'])))
        #winTopCladesDict = {}
        for clade in topCladesDictSorted.keys():
            cladeCount = len(strainsCladesDict[clade] & set(targetDict['strain']))
            #winTopCladesDict[clade] = targetDict['clade'].count(clade) ###need to chancge this to somehow be len(set())
            guideChunkSumDict[clade].append(cladeCount)
    guideChunkSumDF = pd.DataFrame(data=guideChunkSumDict)
    guideChunkSumDFsorted = guideChunkSumDF.sort_values(by=['total_hit'], ascending = False)
    guideChunkSumDFsorted['prop_total_hit'] = guideChunkSumDFsorted['total_hit']/len(chunkMetaDF)
    guideChunkSumDFsorted['prop_strains_hit'] = guideChunkSumDFsorted['strains_hit']/len(set(chunkMetaDF['strain']))
    guideChunkSumDFsorted['prop_subtypes_hit'] = guideChunkSumDFsorted['subtypes_hit']/len(set(chunkMetaDF['subtype']))
    for clade in chunkTopCladesDict.keys():
        guideChunkSumDFsorted['prop_'+clade] = guideChunkSumDFsorted[clade]/chunkTopCladesDict[clade]
    return guideChunkSumDFsorted


chunkFolderList = get_folder_list(folderPathPrefix)
fullGuideDict = {}
for folderPath in chunkFolderList:
    folderName = os.path.basename(folderPath)
    accList = os.listdir(folderPath)
    print("Getting guides from " + folderName)
    chunkDict = getChunkDict(folderPath)
    fullGuideDict[folderName] = chunkDict
    print("Summing guides from " + folderName)
    guideChunkSumDFsorted = sumCountDict(chunkDict, accList)
    guideChunkSumDFsorted.to_csv(folderName+"_guide_hit_summary.tsv", sep="\t", index=False)
    print("Done with folder " + folderName)

print("Done parsing chunk folders")

comboGuideDict = {}
for folder in fullGuideDict.keys():
    print(folder)
    for target in fullGuideDict[folder].keys():
        if target in comboGuideDict.keys():
            for key in comboGuideDict[target].keys():
                comboGuideDict[target][key] = comboGuideDict[target][key] + fullGuideDict[folder][target][key]
        else: 
            comboGuideDict[target] = fullGuideDict[folder][target].copy()

print("Done building combined dictionary")

species = metaDF['spp'][1]
fullAcc = metaDF['accession']
guideSumDFsorted = sumCountDict(comboGuideDict, fullAcc)
species = metaDF['spp'][1]
guideSumDFsorted.to_csv(species + "_all_guides_hit_summary.tsv", sep="\t", index=False)

Path("top100guides").mkdir(parents=True, exist_ok=True)
topTargets = guideSumDFsorted['target'][0:100].to_list()
for target in topTargets:
    targetDict = comboGuideDict[target]
    targetDF = pd.DataFrame(data=targetDict)
    targetAcc = guideSumDFsorted[guideSumDFsorted['target']==target]['accession'].to_string(index=False)
    targetStart = guideSumDFsorted[guideSumDFsorted['target']==target]['start'].to_string(index=False)
    targetDF.to_csv('top100guides/' + targetAcc + "_" + targetStart + "_" + target + "_guide_full_hits.tsv", sep="\t", index=False)
