#!/usr/bin/python3

import glob
import json
import math
import os
from os import path

runProgram="DBSCAN"
checkProgram="DBSCANCheck"
dataDir = "../datasets/"

def MakePBBSTest(dataFile, numprocs, eps, minpts, rho = 0.01):
  return [numprocs, dataFile + ".pbbs","-eps " + str(eps) + " -minpts " + str(minpts) + " -rho " + str(rho),"-eps " + str(eps) + " -minpts " + str(minpts) + " -rho " + str(rho)]

def myround(x, base=1):
    return base * round(x/base)

def ScaleRange(a,b,n):
  m = n-1
  jump = (b-a)/m
  myRet = []
  for idx in range(0,n-1):
    myRet.append(myround(a + jump * idx))
  myRet.append(myround(b))
  return myRet

# given a dictionary of eps range, generate a dic of eps list values
def GenEpsRange(myDict):
  myNewDict = dict()
  for k in myDict:
    myRange = myDict[k]
    if myRange[0] == myRange[1]:
      myNewDict[k] = ScaleRange(max(10, myRange[0] / 2), min(3000,myRange[0] * 2), 6)
    elif "UniformFill" in k:
      myNewDict[k] = ScaleRange(10, 3000, 6)
    else:
      myDiff = myRange[1] - myRange[0]
      myNewDict[k] = ScaleRange(max(10,int(myRange[0] - myDiff)), min(3000,myRange[1] + myDiff), 6)
  return myNewDict

############################################ settings
dataFileIO = open("dataFiles.json","r")
dataFiles = json.load(dataFileIO)
dataFileIO.close()

epsRangeIO = open("epsRange.json", "r")
epsRange = json.load(epsRangeIO)
defaultEpsRange = GenEpsRange(epsRange)

# hard code params for real datasets
defaultEpsRange["7D_HouseHold_2M"] = ScaleRange(max(10, 2000 / 2), min(3000,2000 * 2), 6)
defaultEpsRange["3D_GeoLife_24M"] = [20, 40, 80, 160]
defaultEpsRange["3D_Cosmo50_321M"] = [0.01, 0.02, 0.04, 0.08]
defaultEpsRange["2D_OpenStreetMap_2B"] = [0.01, 0.02, 0.04, 0.08]
defaultEpsRange["13D_TeraClickLog_4B"] = [1500, 3000, 6000, 12000]

epsRangeIO.close()

defaultMinptsRange = [10,100,500,1000,5000,10000]

from myParams import myParams
defaultParam = dict()
for k in myParams:
  myParams[k] = sorted(myParams[k])
  defaultParam[k] = myParams[k][int(math.floor(len(myParams[k]) / 2))]

# hard code params for real datasets
defaultParam["7D_HouseHold_2M"] = (2000, 100)
defaultParam["3D_GeoLife_24M"] = (40, 100)
defaultParam["3D_Cosmo50_321M"] = (0.02, 100)
defaultParam["2D_OpenStreetMap_2B"] = (0.02, 100)
defaultParam["13D_TeraClickLog_4B"] = (3000, 100)

defaultRho = 0.01
defaultRhoRange = [0.001, 0.01, 0.02, 0.05, 0.1]

defaultThreadRange = [1,4,8,16,24,36,72]

fconfig = open("myconfig.json")
config = json.load(fconfig)
defaultThreadCount = config["maxnumprocs"]

if config["syntheticdata"] == "yes":
  dataFiles = dataFiles["synthetic"]
else:
  dataFiles = dataFiles["real"]

############################################

if __name__ == "__main__":
  defaultThreadRange.reverse()

  header = ""
  header += "config-id: " + config["config-id"] + "\n"
  header += "configuration: " + config["configuration"] + "\n"
  header += "machine: " + config["machine"] + "\n"

  tests = []

  def ParseTestsThreadCount(test, defaultValRange):
    if "whitelist" in test:
      for dataset in test["whitelist"]:
        for val in defaultValRange:
          tests.append(MakePBBSTest(dataset, val, defaultParam[dataset][0], defaultParam[dataset][1], defaultRho))
    elif "blacklist" in test:
      for dataset in dataFiles:
        if (config["syntheticdata"] == "yes" and "10M" not in dataset) or (config["ndalgo"] == "2d" and "2D" not in dataset):
          continue
        if dataset in test["blacklist"]:
          continue
        for val in defaultValRange:
          tests.append(MakePBBSTest(dataset, val, defaultParam[dataset][0], defaultParam[dataset][1], defaultRho))

  def ParseTestsEps(test, defaultValRange):
    if "whitelist" in test:
      for dataset in test["whitelist"]:
        for val in defaultValRange[dataset]:
          tests.append(MakePBBSTest(dataset, defaultThreadCount, val, defaultParam[dataset][1], defaultRho))
    elif "blacklist" in test:
      for dataset in dataFiles:
        if (config["syntheticdata"] == "yes" and "10M" not in dataset) or (config["ndalgo"] == "2d" and "2D" not in dataset):
          continue
        if dataset in test["blacklist"]:
          continue
        for val in defaultValRange[dataset]:
          tests.append(MakePBBSTest(dataset, defaultThreadCount, val, defaultParam[dataset][1], defaultRho))

  def ParseTestsMinpts(test, defaultValRange):
    if "whitelist" in test:
      for dataset in test["whitelist"]:
        for val in defaultValRange:
          tests.append(MakePBBSTest(dataset, defaultThreadCount, defaultParam[dataset][0], val, defaultRho))
    elif "blacklist" in test:
      for dataset in dataFiles:
        if (config["syntheticdata"] == "yes" and "10M" not in dataset) or (config["ndalgo"] == "2d" and "2D" not in dataset):
          continue
        if dataset in test["blacklist"]:
          continue
        for val in defaultValRange:
          tests.append(MakePBBSTest(dataset, defaultThreadCount, defaultParam[dataset][0], val, defaultRho))

  def ParseTestsRho(test, defaultValRange):
    if "whitelist" in test:
      for dataset in test["whitelist"]:
        for val in defaultValRange:
          tests.append(MakePBBSTest(dataset, defaultThreadCount, defaultParam[dataset][0], defaultParam[dataset][1], val))
    elif "blacklist" in test:
      for dataset in dataFiles:
        if (config["syntheticdata"] == "yes" and "10M" not in dataset) or (config["ndalgo"] == "2d" and "2D" not in dataset):
          continue
        if dataset in test["blacklist"]:
          continue
        for val in defaultValRange:
          tests.append(MakePBBSTest(dataset, defaultThreadCount, defaultParam[dataset][0], defaultParam[dataset][1], val))
          
  def ParseTestsNumPts(test):
    if "whitelist" in test:
      print("warning: numpts tests do not support whitelist")
    elif "blacklist" in test:
      for dataset in dataFiles:
        if (dataset in test["blacklist"]) or (config["ndalgo"] == "2d" and "2D" not in dataset):
          continue
        tests.append(MakePBBSTest(dataset, defaultThreadCount, defaultParam[dataset][0], defaultParam[dataset][1], defaultRho))

  if "testnumpts" in config:
    ParseTestsNumPts(config["testnumpts"])

  if "testeps" in config:
    ParseTestsEps(config["testeps"], defaultEpsRange)

  if "testminpts" in config:
    ParseTestsMinpts(config["testminpts"], defaultMinptsRange)
    
  if "testrho" in config:
    ParseTestsRho(config["testrho"], defaultRhoRange)
    
  if "testthreadcount" in config:
    ParseTestsThreadCount(config["testthreadcount"], defaultThreadRange)

  #todo eliminate repeated tests

  b4prune = len(tests)
  nottested = list()
  pruned = list()
  for test in tests:
    if test not in pruned:
      if not path.exists(dataDir + test[1]):
        print("warning: " + dataDir + test[1] + " does not exist")
        if test not in nottested:
          nottested.append(test)
      else:
        pruned.append(test)
  print("tests cannot be done: ")
  for test in nottested:
    print(test)
  print("tests that will be done: ")
  for test in pruned:
    print(test)
  print("num tests before pruning and checking " + str(b4prune))
  print("num tests remains " + str(len(pruned)))

  import runDBSCANConfig
  runDBSCANConfig.timeAllArgs(runProgram, checkProgram, dataDir, pruned, header)

