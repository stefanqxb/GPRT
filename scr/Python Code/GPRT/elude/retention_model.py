import numpy as np
import data_manager as dm

def buildRetentionIndex(aaAlphabet, psmDescriptions, 
      normalizeRetentionTimes):
  if normalizeRetentionTimes:
    print "normalizeRetentionTimes is not implemented yet"
  aaAlphabet = sorted(list(aaAlphabet))
  featureMatrix = computeRetentionFeatureMatrix(aaAlphabet, psmDescriptions)
  normalizeFeatures(featureMatrix)
  
  retentionTimes = [psmd.retentionTime for psmd in psmDescriptions]
  
  print featureMatrix[0]
  
  b = np.linalg.lstsq(featureMatrix, retentionTimes)[0]
  return b

def computeRetentionFeatureMatrix(aaAlphabet, psmDescriptions):
  featureMatrix = np.zeros((len(psmDescriptions), len(aaAlphabet)))
  for i, psmd in enumerate(psmDescriptions):
    featureMatrix[i] = computeRetentionFeatures(aaAlphabet, psmd)
  return featureMatrix
  
def computeRetentionFeatures(aaAlphabet, psmd):
  aas = dm.getAminoAcidList(psmd.peptide)
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
    
