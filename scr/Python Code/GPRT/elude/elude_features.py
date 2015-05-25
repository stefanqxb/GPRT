import data_manager as dm
import retention_model as rm
import random
import numpy as np

removeDuplicates = True
removeCommonPeptides = False # not implemented yet
removeInSourceFragments = False # not implemented yet
removeNonEnzymatic = False # not implemented yet
normalizeRetentionTimes = True

def processTrainData(trainFile, testFile = ""):
  trainPsms, trainAaAlphabet = dm.loadPeptides(trainFile)
  if removeDuplicates:
    trainPsms = dm.removeDuplicates(trainPsms)
  
  if rm.hasPtms(trainAaAlphabet):
    trainAaAlphabet = rm.defaultAlphabet
  
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
  customIndex = rm.buildRetentionIndex(aaAlphabet, psmDescriptions, normalizeRetentionTimes)
  return dict(zip(aaAlphabet, customIndex))

def trainRetentionModelUnitTest():
  psmDescriptions, aaAlphabet = processTrainData('../data/retention_time_peptide.csv')
  print len(psmDescriptions), len(aaAlphabet), psmDescriptions[0].peptide, psmDescriptions[0].retentionTime
  
  trainedIndex = trainRetentionModel(aaAlphabet, psmDescriptions)
  
  for aa, index in trainedIndex:
    print aa,':',round(index,6)

def getFeatures(trainFile):
  psmDescriptions, aaAlphabet = processTrainData(trainFile)
  customIndex = trainRetentionModel(aaAlphabet, psmDescriptions)
  featureMatrix = rm.computeRetentionFeatureMatrix(aaAlphabet, psmDescriptions, customIndex)
  return psmDescriptions, featureMatrix

def getFeaturesUnitTest():
  psmDescriptions, featureMatrix = getFeatures('../data/retention_time_peptide.csv')
  
  for idx in [100,101,102]:
    print ""
    if idx == 100:
      # feature vector for K.LCNNQEENDAVSSAK.K
      eludeFeatureVector = [0.55743, 0.392982, 0.922222, 0.0666667, 0.0974441, 0.571038, 0.798507, 0.0337838, 0.116279, 0.397471, 0.21625, 0.38153, 0.176623, 0.0741797, 0.00524934, 0.119266, 0.296296, 0.25, 0.2, 0.111111, 0.188252, 0.263656, 1, 0.126391, 0.0666641, 0.221891, 0.496811, 0.356832, 0.456211, 0.354441, 0.289677, 0.214665, 0.274368, 0.123323, 0.0191546, 0.0831153, 0.181818, 0.0666667, 0.142857, 0, 0.193677, 0.204545, 0.105263, 0.142857, 0.0555556, 0.111111, 0, 0, 0, 0, 0.333333, 0.0909091, 0, 0.375, 0, 0.0909091, 0, 0.153846, 0, 0.125, 0, 0]
    elif idx == 101:
      # feature vector for R.LGTPALTSR.G
      eludeFeatureVector = [0.690763, 0.592593, 0.922222, 0, 0, 0.480874, 0.746269, 0.611486, 0.313953, 0.477686, 0.655993, 0, 0.0127273, 0.00215584, 0.0128609, 0.0563497, 0.037037, 0, 0.133333, 0, 0.206105, 0.402284, 1, 0.140665, 0.110162, 0.308878, 0.601676, 0.547744, 0.373248, 0.322974, 0.489505, 0, 0.0285554, 0.00837665, 0.0398314, 0.142634, 0.0909091, 0, 0.142857, 0, 0.0848419, 0.0681818, 0.0526316, 0, 0, 0, 0, 0.0526316, 0, 0, 0, 0.181818, 0, 0, 0.0909091, 0, 0.333333, 0.0769231, 0.222222, 0, 0, 0]
    else:
      # feature vector for K.TVIVAALDGTFQR.K
      eludeFeatureVector = [0.763855, 0.711201, 0.422222, 0, 0.105431, 0.852459, 0.977612, 0.412162, 4.13106e-16, 0.840601, 0.473735, 0.528567, 0.25974, 0.317542, 0.152231, 0.132353, 0.111111, 0.0416667, 0.333333, 0.222222, 0.284283, 0.438798, 0.336514, 0.140665, 0.0294914, 0.506548, 0.623036, 0.400291, 0.35311, 0.661824, 0.276678, 0.564349, 0.250851, 0.449059, 0.306076, 0.115925, 0.0909091, 0, 0.357143, 0.285714, 0.191194, 0.159091, 0.105263, 0, 0.0555556, 0, 0.2, 0.0526316, 0, 0.125, 0, 0.0909091, 0, 0, 0, 0.0909091, 0.333333, 0, 0.222222, 0.25, 0, 0]
      
    print "Peptide:", psmDescriptions[idx].peptide
    # print psmDescriptions[idx].peptide, len(featureMatrix[idx]), featureMatrix[idx][:20]
    
    i = 1
    errors = 0
    for eludeFeature, ourFeature in zip(eludeFeatureVector, featureMatrix[idx]):
      if abs(eludeFeature - ourFeature) > 1e-5 and i not in range(21,41):
        print "Wrong feature:", i, ", EludeFeature:", eludeFeature, ", OurFeature:", ourFeature
        errors += 1
      i += 1
    if errors == 0:
      print "Unit test succeeded"

