__author__ = 'Quinn'

import os
import argparse
#from re import split
import pandas as pd
import numpy

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
parser.add_argument('--winCounts', required=True, 
                    help="counts of strains per each 20bp window based on ref alignment, required")
parser.add_argument('--pathPrefix', required=True,
                    help='The prefix of the folders with the chunked genomes')
#parser.add_argument('--out', help="Output prefix, required", required=True)
args = parser.parse_args()

folderPathPrefix=args.pathPrefix
#folderPathPrefix='tiled_genomes/tiled_subtype_'
#metaDF = pd.read_csv('RSV-B_strains_metadata.txt', sep="\t")
#winCountsDF = pd.read_csv('RSV-B_strainCountPerWindow.tsv', sep="\t")

#load genome dataframe and generate genomeID-to-segment lookup dictionary
metaDF = pd.read_csv(args.metadata, sep="\t")
metaDict = metaDF.set_index('accession').T.to_dict('list')
startDict = metaDF.set_index('accession')['alignmentStart'].to_dict()
topClades = list(set(metaDF[metaDF['topSubtype']]['subtype']))
topCladesDict = {}
for clade in topClades:
     count = len(set(metaDF[metaDF['subtype']==clade]['strain']))
     topCladesDict[clade] = count

topCladesDictSorted = dict(sorted(topCladesDict.items(), key=lambda x:x[1], reverse=True))

winCountsDF = pd.read_csv(args.winCounts, sep="\t")
countDict = {}
countDict['totalCount'] = winCountsDF.set_index('winStart')['strainCount'].to_dict()
countDict['uniqueCount'] = winCountsDF.set_index('winStart')['strainUni'].to_dict()
countDict['allOthers'] = winCountsDF.set_index('winStart')['allOthers'].to_dict()
for clade in topClades:
     countDict[clade] = winCountsDF.set_index('winStart')[clade].to_dict()


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
        curLen = curMeta[0]
        curStart = curMeta[1]
        curEnd = curMeta[2]
        curGB = curMeta[3]
        curSpp = curMeta[4]
        curSeg = curMeta[5]
        curClade = curMeta[6]
        curStrain = curMeta[7]
        curCountry = curMeta[8]
        curHost = curMeta[9]
        curDate = curMeta[10]
        curTop = curMeta[11]
        curFolder = curMeta[12]
        windowsPath =  os.path.join(folderBase, pathStep, "windows.txt")
        windowsFile = open(windowsPath, 'r')
        windowLines = windowsFile.readlines()[1:]
        for line in windowLines:
            winLine = line.strip().split("\t")
            winAcc = winLine[0]
            winStart = winLine[1]
            refStart =  int(winStart) + (curStart - 1)
            winTarget = winLine[2]
            winSpacer = winLine[3]
            winStrand = winLine[4]
            winGC = winLine[5]
            winA = winLine[6]
            if float(winGC)>0.25 and float(winGC)<0.75:
                if winTarget in guideCountDict.keys():
                    guideCountDict[winTarget]['accW'].append(winAcc)
                    guideCountDict[winTarget]['start'].append(winStart)
                    guideCountDict[winTarget]['refStart'].append(refStart)
                    guideCountDict[winTarget]['target'].append(winTarget)
                    guideCountDict[winTarget]['spacer'].append(winSpacer)
                    guideCountDict[winTarget]['strand'].append(winStrand)
                    guideCountDict[winTarget]['GC'].append(winGC)
                    guideCountDict[winTarget]['A'].append(winA)
                    guideCountDict[winTarget]['accM'].append(curAcc)
                    guideCountDict[winTarget]['length'].append(curLen)
                    guideCountDict[winTarget]['aliStart'].append(curStart)
                    guideCountDict[winTarget]['aliEnd'].append(curEnd)
                    guideCountDict[winTarget]['gb'].append(curGB)
                    guideCountDict[winTarget]['spp'].append(curSpp)
                    guideCountDict[winTarget]['seg'].append(curSeg)
                    guideCountDict[winTarget]['clade'].append(curClade)
                    guideCountDict[winTarget]['strain'].append(curStrain)
                    guideCountDict[winTarget]['country'].append(curCountry)
                    guideCountDict[winTarget]['host'].append(curHost)
                    guideCountDict[winTarget]['date'].append(curDate)
                    guideCountDict[winTarget]['topSub'].append(curTop)
                    guideCountDict[winTarget]['subFolder'].append(curFolder)
                else:
                    winDict = {'accW':[winAcc], 'start':[winStart], 'refStart':[refStart], 'target':[winTarget], 
                            'spacer':[winSpacer], 'strand':[winStrand], 'GC':[winGC], 'A':[winA], 'accM':[curAcc], 
                            'length':[curLen], 'aliStart':[curStart], 'aliEnd':[curEnd], 'gb':[curGB], 'spp':[curSpp], 
                            'seg':[curSeg], 'clade':[curClade], 'strain':[curStrain], 'country':[curCountry], 
                            'host':[curHost], 'date':[curDate], 'topSub':[curTop], 'subFolder':[curFolder]}
                    guideCountDict[winTarget] = winDict
        windowsFile.close()
    return guideCountDict

def sumCountDict(guideCountDict, accList, subPop):
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
    strainsCladesDict['allOthers'] = set(chunkMetaDF[chunkMetaDF['subtypeFolder']=='subtype_allOthers']['strain'])
    chunkTopCladesDict['allOthers'] = len(set(set(chunkMetaDF[chunkMetaDF['subtypeFolder']=='subtype_allOthers']['strain'])))
    guideChunkSumDict = {'accession':[], 'refStart':[], 'target':[], 'total_hit':[], 'prop_win_total_hit':[],
                          'strains_hit':[], 'prop_win_strains_hit':[], 'spacer':[], 'strand':[], 'GC_content':[], "A_content":[]}    
    if subPop == "all":
        for clade in topCladesDictSorted.keys():
            guideChunkSumDict[clade] = []
            guideChunkSumDict['prop_win_'+clade] = []
        guideChunkSumDict['allOthers'] = []
        guideChunkSumDict['prop_win_allOthers'] = []
        guideChunkSumDict['subtypes_hit'] = []
    else:
        guideChunkSumDict[subPop] = []
        guideChunkSumDict['prop_win_'+subPop] = []
    for target in guideCountDict.keys():
        #print(target)
        targetDict = guideCountDict[target]
        firstAcc = sorted(targetDict['accW'])[0]
        guideChunkSumDict['accession'].append(firstAcc)
        firstIndex = targetDict['accW'].index(firstAcc)
        refPos = targetDict['refStart'][firstIndex]
        if refPos > max(countDict['totalCount'].keys()):
            refPos = max(countDict['totalCount'].keys()) ##Some fa are longer than the ref...
        elif numpy.isnan(refPos):
            refPos = 1
        winTotCount = countDict['totalCount'][refPos]
        winStrainCount = countDict['uniqueCount'][refPos]
        guideChunkSumDict['refStart'].append(refPos)
        guideChunkSumDict['target'].append(targetDict['target'][firstIndex])
        guideChunkSumDict['spacer'].append(targetDict['spacer'][firstIndex])
        guideChunkSumDict['strand'].append(targetDict['strand'][firstIndex])
        guideChunkSumDict['A_content'].append(targetDict['A'][firstIndex])
        guideChunkSumDict['GC_content'].append(targetDict['GC'][firstIndex])
        guideChunkSumDict['total_hit'].append(len(set(targetDict['accM'])))
        guideChunkSumDict['strains_hit'].append(len(set(targetDict['strain'])))
        guideChunkSumDict['prop_win_total_hit'].append(len(set(targetDict['accM']))/winTotCount)
        guideChunkSumDict['prop_win_strains_hit'].append(len(set(targetDict['strain']))/winStrainCount)
        #winTopCladesDict = {}
        if subPop == "all":
            for clade in topCladesDictSorted.keys():
                cladeCount = len(strainsCladesDict[clade] & set(targetDict['strain']))
                guideChunkSumDict[clade].append(cladeCount)
                cladeWinCount = countDict[clade][refPos]
                guideChunkSumDict['prop_win_'+clade].append(cladeCount/cladeWinCount)
            cladeCount = len(strainsCladesDict['allOthers'] & set(targetDict['strain']))
            guideChunkSumDict['allOthers'].append(cladeCount)
            cladeWinCount = countDict['allOthers'][refPos]
            guideChunkSumDict['prop_win_allOthers'].append(cladeCount/cladeWinCount)
            guideChunkSumDict['subtypes_hit'].append(len(set(targetDict['clade'])))
        else:
            clade = subPop
            cladeCount = len(strainsCladesDict[clade] & set(targetDict['strain']))
            guideChunkSumDict[clade].append(cladeCount)
            cladeWinCount = countDict[clade][refPos]
            guideChunkSumDict['prop_win_'+clade].append(cladeCount/cladeWinCount)
    guideChunkSumDF = pd.DataFrame(data=guideChunkSumDict)
    guideChunkSumDFsorted = guideChunkSumDF.sort_values(by=['totaql_hit'], ascending = False)
    guideChunkSumDFsorted['prop_total_hit'] = guideChunkSumDFsorted['total_hit']/len(chunkMetaDF)
    guideChunkSumDFsorted['prop_strains_hit'] = guideChunkSumDFsorted['strains_hit']/len(set(chunkMetaDF['strain']))
    if subPop == "all":
        for clade in chunkTopCladesDict.keys():
            guideChunkSumDFsorted['prop_'+clade] = guideChunkSumDFsorted[clade]/chunkTopCladesDict[clade]
        guideChunkSumDFsorted['prop_subtypes_hit'] = guideChunkSumDFsorted['subtypes_hit']/len(set(chunkMetaDF['subtype']))
    else:
        guideChunkSumDFsorted['prop_'+subPop] = guideChunkSumDFsorted[subPop]/chunkTopCladesDict[subPop]
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
    subPop = folderName.split("_")[-1]
    guideChunkSumDFsorted = sumCountDict(chunkDict, accList, subPop)
    guideChunkSumDFsorted.to_csv(folderName+"_guides_hit_summary.tsv", sep="\t", index=False)
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
guideSumDFsorted = sumCountDict(comboGuideDict, fullAcc, 'all')
guideSumDFsorted.to_csv(species + "_all_guides_hit_summary.tsv", sep="\t", index=False)
winSortDF = guideSumDFsorted.sort_values(by=['prop_win_strains_hit'], ascending = False)

topTargets = guideSumDFsorted['target'][0:100].to_list()
for target in topTargets:
    targetDict = comboGuideDict[target]
    targetDF = pd.DataFrame(data=targetDict)
    targetAcc = guideSumDFsorted[guideSumDFsorted['target']==target]['accession'].to_string(index=False)
    targetStart = guideSumDFsorted[guideSumDFsorted['target']==target]['refStart'].to_string(index=False)
    targetDF.to_csv('top100guides/' + targetAcc + "_" + targetStart + "_" + target + "_guide_full_hits.tsv", sep="\t", index=False)
