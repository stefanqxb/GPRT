import numpy as np
import data_manager as dm

kyteDoolittleIndex = { "A" : 1.8, "C" : 2.5, "D" : -3.5, "E" : -3.5, "F" : 2.8,
   "G" : -0.4, "H" : -3.2, "I" : 4.5, "K" : -3.9, "L" : 3.8,
   "M" : 1.9, "N" : -3.5, "P" : -1.6, "Q" : -3.5, "R" : -4.5,
   "S" : -0.8, "T" : -0.7, "V" : 4.2, "W" : -0.9, "Y" : -1.3};

bulkinessIndex = { "A" : 11.5, "C" : 13.46, "D" : 11.68, "E" : 13.57, "F" : 19.80,
   "G" : 3.40, "H" : 13.69, "I" : 21.40, "K" : 15.71, "L" : 21.40,
   "M" : 16.25, "N" : 12.82, "P" : 17.43, "Q" : 14.45, "R" : 14.28,
   "S" : 9.47, "T" : 15.77, "V" : 21.57, "W" : 21.67, "Y" : 18.03};

defaultAlphabet = set(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])

percentageAa = 0.25

def buildRetentionIndex(aaAlphabet, psmDescriptions, 
      normalizeRetentionTimes):
  featureMatrix = computeRetentionIndexFeatureMatrix(aaAlphabet, psmDescriptions)
  normalizeFeatures(featureMatrix)
  
  retentionTimes = [psmd.retentionTime for psmd in psmDescriptions]
  if normalizeRetentionTimes:
    retentionTimes -= np.mean(retentionTimes)
    retentionTimes /= np.std(retentionTimes)
  
  customIndex = np.linalg.lstsq(featureMatrix, retentionTimes)[0]
  return customIndex

def computeRetentionIndexFeatureMatrix(aaAlphabet, psmDescriptions):
  featureMatrix = np.zeros((len(psmDescriptions), len(aaAlphabet)))
  for i, psmd in enumerate(psmDescriptions):
    featureMatrix[i] = computeRetentionIndexFeatures(aaAlphabet, psmd.peptide)
  return featureMatrix
  
def computeRetentionIndexFeatures(aaAlphabet, peptide):
  aas = dm.getAminoAcidList(peptide)
  featureVector = np.zeros((1, len(aaAlphabet)))
  for aa in aas:
    featureVector[0][aaAlphabet.index(aa)] += 1
  return featureVector

def normalizeFeatures(featureMatrix):
  rows, cols = featureMatrix.shape
  colMean = list()
  colStd = list()
  for i in range(cols):
    featureMatrix[:,i] -= np.mean(featureMatrix[:,i])
    featureMatrix[:,i] /= np.std(featureMatrix[:,i])

def hasPtms(aaAlphabet):
  return sum([1 for aa in aaAlphabet if aa not in defaultAlphabet]) > 0
  
def computeRetentionFeatureMatrix(aaAlphabet, psmDescriptions, customIndex):
  ptmsPresent = False
  if hasPtms(aaAlphabet):
    numFeatures = 20 + 20 + 2 + len(aaAlphabet)
    ptmsPresent = True
  else:
    numFeatures = 20 + 1 + len(aaAlphabet)
  
  polarAa, hydrophobicAa = getExtremeRetentionAA(customIndex)
  
  featureMatrix = np.zeros((len(psmDescriptions), 20))
  for i, psmd in enumerate(psmDescriptions):
    featureMatrix[i] = computeIndexFeatures(aaAlphabet, psmd.peptide, customIndex, polarAa, hydrophobicAa)
  
  if not ptmsPresent:
    polarAa, hydrophobicAa = getExtremeRetentionAA(kyteDoolittleIndex)    
    kyteDoolittleFeatureMatrix = np.zeros((len(psmDescriptions), 20))
    for i, psmd in enumerate(psmDescriptions):
      kyteDoolittleFeatureMatrix[i] = computeIndexFeatures(aaAlphabet, psmd.peptide, kyteDoolittleIndex, polarAa, hydrophobicAa)
    featureMatrix = np.concatenate((featureMatrix, kyteDoolittleFeatureMatrix), axis = 1)
  
  return featureMatrix

def getExtremeRetentionAA(index):
  numAa = int(np.ceil(percentageAa * len(index)))
  sortedIndex = sorted(index.items(), key = lambda x : x[1])
  return [x[0] for x in sortedIndex[:numAa]], [x[0] for x in sortedIndex[-1*numAa:]]
  
def computeIndexFeatures(aaAlphabet, peptide, index, polarAa, hydrophobicAa):
  features = []
  aas = dm.getAminoAcidList(peptide)
  features.append(indexSum(aas, index))
  features.append(indexAvg(aas, index))
  features.append(indexN(peptide, index))
  features.append(indexC(peptide, index))
  features.append(indexNearestNeigbour(peptide, index, polarAa))
  features.append(indexMaxPartialSum(peptide, index, 5))
  features.append(indexMaxPartialSum(peptide, index, 2))
  features.append(indexMinPartialSum(peptide, index, 5))
  features.append(indexMinPartialSum(peptide, index, 2))
  features.append(indexMaxHydrophobicSideHelix(peptide, index))
  features.append(indexMinHydrophobicSideHelix(peptide, index))
  features.append(indexMaxHydrophobicMoment(peptide, index, 100, 11))
  features.append(indexMaxHydrophobicMoment(peptide, index, 180, 11))
  features.append(indexMinHydrophobicMoment(peptide, index, 100, 11))
  features.append(indexMinHydrophobicMoment(peptide, index, 180, 11))
  features.append(indexSumSquaredDiff(peptide, index))
  features.append(numberTypeAA(peptide, polarAa))
  features.append(numberConsecTypeAA(peptide, polarAa))
  features.append(numberTypeAA(peptide, hydrophobicAa))
  features.append(numberConsecTypeAA(peptide, hydrophobicAa))
  return features

def indexSum(aas, index):
  return sum([index[aa] for aa in aas])

def indexAvg(peptide, index):
  return 0
  
def indexN(peptide, index):
  return 0

def indexC(peptide, index):
  return 0

def indexNearestNeigbour(peptide, index, polarAa):
  return 0

def indexMaxPartialSum(peptide, index, window):
  return 0

def indexMinPartialSum(peptide, index, window):
  return 0

def indexMaxHydrophobicSideHelix(peptide, index):
  return 0

def indexMinHydrophobicSideHelix(peptide, index):
  return 0

def indexMaxHydrophobicMoment(peptide, index, angle, window):
  return 0

def indexMinHydrophobicMoment(peptide, index, angle, window):
  return 0

def indexSumSquaredDiff(peptide, index):
  return 0

def numberTypeAA(peptide, aas):
  return 0

def numberConsecTypeAA(peptide, aas):
  return 0
