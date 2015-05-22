import data_manager as dm
import retention_model as rm
import random
import numpy as np

removeDuplicates = True
removeCommonPeptides = False # not implemented yet
removeInSourceFragments = False # not implemented yet
removeNonEnzymatic = False # not implemented yet
normalizeRetentionTimes = False # not implemented yet

def processTrainData(trainFile, testFile = ""):
  trainPsms, trainAaAlphabet = dm.loadPeptides(trainFile)
  if removeDuplicates:
    trainPsms = dm.removeDuplicates(trainPsms)
  
  if len(testFile) > 0:
    testPsms, testAaAlphabet = dm.loadPeptides(testFile)
    if removeDuplicates:
      testPsms = dm.removeDuplicates(testPsms)
    if removeCommonPeptides:
      print "removeCommonPeptides is not implemented yet"
  
  if removeInSourceFragments:
    print "removeInSourceFragments is not implemented yet"
  
  if removeNonEnzymatic:
    print "removeNonEnzymatic is not implemented yet"
  
  rnd = random.Random(1)
  rnd.shuffle(trainPsms)
  return trainPsms, trainAaAlphabet

def processTrainDataUnitTest():
  psmDescriptions, aaAlphabet = processTrainData('../data/retention_time_peptide.csv')
  print len(psmDescriptions), len(aaAlphabet), psmDescriptions[0].peptide, psmDescriptions[0].retentionTime

def trainRetentionModel(aaAlphabet, psmDescriptions):
  featureMatrix = rm.buildRetentionIndex(aaAlphabet, psmDescriptions, normalizeRetentionTimes)
  return featureMatrix

def trainRetentionModelUnitTest():
  psmDescriptions, aaAlphabet = processTrainData('../data/retention_time_peptide.csv')
  print len(psmDescriptions), len(aaAlphabet), psmDescriptions[0].peptide, psmDescriptions[0].retentionTime
  
  trainedIndex = trainRetentionModel(aaAlphabet, psmDescriptions)
  print trainedIndex
