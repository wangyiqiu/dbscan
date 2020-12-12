#!/usr/bin/python3

import glob
import json
import math
import os
from os import path
from dbscanParse import ParseTestFile
from genTest import *
import pprint

# simply find first found
def FindTest(t_serialalgo, t_ndalgo, t_algo, t_tests, t_dataset, t_eps, t_minpts, t_threadcount):
  for test in t_tests:
    if test.GetDatasetName().split(".")[0]==t_dataset and test.m_eps==t_eps and test.m_minpts==t_minpts and test.m_numprocs==t_threadcount and t_serialalgo==test.m_serialalgo and t_ndalgo==test.m_ndalgo and t_algo==test.m_algo:
      # if "7D_UniformFill_10M" in t_dataset and "7D_UniformFill_10M" in test.GetDatasetName():
      #   print(t_dataset)
      #   print(test.GetDatasetName(),test.m_eps,test.m_minpts,test.m_numprocs)
      #   print("match!")
      if test.m_totaltime < 0:
        #print("totaltime",test.m_totaltime)
        #return None
        continue
      else:
        return test
    # else:
    #   if "7D_UniformFill_10M" in t_dataset and "7D_UniformFill_10M" in test.GetDatasetName():
    #     print(t_dataset)
    #     print(test.GetDatasetName(),test.m_eps,test.m_minpts,test.m_numprocs)
    #     print("no match")
  return None

def FindTestRho(t_serialalgo, t_ndalgo, t_algo, t_tests, t_dataset, t_eps, t_minpts, t_rho, t_threadcount):
  for test in t_tests:
    if test.GetDatasetName().split(".")[0]==t_dataset and test.m_eps==t_eps and test.m_minpts==t_minpts and test.m_numprocs==t_threadcount and t_serialalgo==test.m_serialalgo and t_ndalgo==test.m_ndalgo and t_algo==test.m_algo and t_rho==test.m_rho:
      if test.m_totaltime < 0:
        return None
      else:
        return test
  return None

def FindTestByDimName(t_serialalgo, t_ndalgo, t_algo, t_tests, t_dimname, t_eps, t_minpts, t_threadcount, t_numpts):
  for test in t_tests:
    if test.GetDimName()==t_dimname and test.m_eps==t_eps and test.m_minpts==t_minpts and test.m_numprocs==t_threadcount and test.m_numpts==t_numpts and t_serialalgo==test.m_serialalgo and t_ndalgo==test.m_ndalgo and t_algo==test.m_algo:
      if test.m_totaltime < 0:
        return None
      else:
        return test
  return None

def FilterTest(t_tests, t_dataset):
  return [i for i in t_tests if i.GetDatasetName()==t_dataset]

# a snapshot is an aggregation of all tests on one dataset
def CreateTestSnapshot(t_dataset, t_varAttribute, config):
  snapshot = dict()
  snapshot["config-id"] = config["config-id"]
  snapshot["configuration"] = config["configuration"]
  snapshot["machine"] = config["machine"]
  snapshot["dataset"] = t_dataset
  snapshot["dimname"] = None
  snapshot["testcount"] = 0
  # variable attribute
  snapshot["variable"] = t_varAttribute
  snapshot["x"] = list()
  #fixed attributes
  snapshot["serialalgo"] = set()
  snapshot["ndalgo"] = set()
  snapshot["algo"] = set()
  snapshot["datafile"] = set()
  snapshot["eps"] = set()
  snapshot["rounds"] = set()
  snapshot["minpts"] = set()
  snapshot["rho"] = set()
  snapshot["dim"] = set()
  snapshot["numprocs"] = set()
  snapshot["numpts_attr"] = set()
  # test attributes
  snapshot["assigngridtime"] = list()
  snapshot["markcoretime"] = list()
  snapshot["buildcoreptquadtree"] = list()
  snapshot["buildgetnbrkdt"] = list()
  snapshot["clustercoretime"] = list()
  snapshot["assignbordertime"] = list()
  snapshot["totaltime"] = list()
  snapshot["numgrids"] = list()
  snapshot["numcore"] = list()
  snapshot["numborder"] = list()
  snapshot["numclusterpts"] = list()
  snapshot["numpts"] = list()
  snapshot["numclusters"] = list()
  snapshot["numnoise"] = list()
  return snapshot

# it is user's responsibility to aggregate the correct test
def AggregateTest(snapshot, t_varAttribute, t_test):
  # results
  snapshot["assigngridtime"].append(t_test.m_assigngridtime)
  snapshot["markcoretime"].append(t_test.m_markcoretime)
  snapshot["buildcoreptquadtree"].append(t_test.m_buildcoreptquadtree)
  snapshot["buildgetnbrkdt"].append(t_test.m_buildgetnbrkdt)
  snapshot["clustercoretime"].append(t_test.m_clustercoretime)
  snapshot["assignbordertime"].append(t_test.m_assignbordertime)
  snapshot["totaltime"].append(t_test.m_totaltime)
  snapshot["numgrids"].append(t_test.m_numgrids)
  snapshot["numcore"].append(t_test.m_numcore)
  snapshot["numborder"].append(t_test.m_numborder)
  snapshot["numclusterpts"].append(t_test.m_numclusterpts)
  snapshot["numpts"].append(t_test.m_numpts)
  snapshot["numclusters"].append(t_test.m_numclusters)
  snapshot["numnoise"].append(t_test.m_numnoise)
  # attributes
  snapshot["serialalgo"].add(t_test.m_serialalgo)
  snapshot["ndalgo"].add(t_test.m_ndalgo)
  snapshot["algo"].add(t_test.m_algo)
  snapshot["datafile"].add(t_test.m_datafile)
  snapshot["dimname"] = t_test.GetDimName()
  snapshot["eps"].add(t_test.m_eps)
  snapshot["rounds"].add(t_test.m_rounds)
  snapshot["minpts"].add(t_test.m_minpts)
  snapshot["rho"].add(t_test.m_rho)
  snapshot["dim"].add(t_test.m_dim)
  snapshot["numprocs"].add(t_test.m_numprocs)
  snapshot["numpts_attr"].add(t_test.m_numpts)
  # variable
  if t_varAttribute=="eps":
    snapshot["x"].append(t_test.m_eps)
  elif t_varAttribute=="minpts":
    snapshot["x"].append(t_test.m_minpts)
  elif t_varAttribute=="rho":
    snapshot["x"].append(t_test.m_rho)
  elif t_varAttribute=="numprocs":
    snapshot["x"].append(t_test.m_numprocs)
  elif t_varAttribute=="numpts_attr":
    snapshot["x"].append(t_test.m_numpts)
  else:
    print("warning: var attribute illegal")
    return snapshot
  snapshot["testcount"] += 1
  return snapshot

def CleanupSnapshot(snapshot):
  return snapshot

#####################################################3

# read files, for pbbs output
textFiles = glob.glob("./outputs/*.txt")
print(textFiles)
tests = list()
for textFile in textFiles:
  tests += ParseTestFile(textFile)
  
'''
# read files, for hpdbscan output
from hpParser import *
for test in tests:
  print(test.m_datafile + " " + test.GetMyBasicInfo())

# read files, for pdsdbscan output
from pdsParser import *
for test in tests:
  print(test.m_datafile + " " + test.GetMyBasicInfo())
'''

# indicate what test to parse
serialalgo="parallel"
ndalgo="nd"
algo="exact"
#####################################################3

plotsData = dict()
plotsData["numprocsPlots"] = dict()
plotsData["epsPlots"] = dict()
plotsData["minptsPlots"] = dict()
plotsData["numptsPlots"] = dict()
plotsData["rhoPlots"] = dict()

for targetDataname in dataFiles:
  if (config["syntheticdata"] == "yes" and "10M" not in targetDataname) or (config["ndalgo"] == "2d" and "2D" not in targetDataname):
    continue
  targetDataset = targetDataname + ".pbbs"
  targetTests = FilterTest(tests, targetDataset)

  mysnapshot = CreateTestSnapshot(targetDataset, "numprocs", config)
  for tcount in defaultThreadRange:
    # if "7D_UniformFill_10M" not in targetDataname:
    #   continue
    # else:
    #print(targetDataname,"!!!!!!!!!!!!!!!!!")
    #print("find",targetDataname,defaultParam[targetDataname][0], defaultParam[targetDataname][1], tcount);
    found =FindTest(serialalgo, ndalgo, algo, tests, targetDataname, defaultParam[targetDataname][0], defaultParam[targetDataname][1], tcount)
    if found==None:
      print("warning: threadcount did not find " + targetDataset + " tcount " + str(tcount) + " eps " + str(defaultParam[targetDataname][0]) + " minpts " + str(defaultParam[targetDataname][1]))
    else:
      AggregateTest(mysnapshot, "numprocs", found)
  plotsData["numprocsPlots"][targetDataname] = CleanupSnapshot(mysnapshot)
  #continue

  mysnapshot = CreateTestSnapshot(targetDataset, "eps", config)
  for eps in defaultEpsRange[targetDataname]:
    found =FindTest(serialalgo, ndalgo, algo, tests, targetDataname, eps, defaultParam[targetDataname][1], defaultThreadCount)
    if found==None:
      print("warning: eps did not find " + targetDataset + " tcount " + str(defaultThreadCount) + " eps " + str(eps) + " minpts " + str(defaultParam[targetDataname][1]))
    else:
      AggregateTest(mysnapshot, "eps", found)
  plotsData["epsPlots"][targetDataname] = CleanupSnapshot(mysnapshot)

  mysnapshot = CreateTestSnapshot(targetDataset, "minpts", config)
  for val in defaultMinptsRange:
    found =FindTest(serialalgo, ndalgo, algo, tests, targetDataname, defaultParam[targetDataname][0], val, defaultThreadCount)
    if found==None:
      print("warning: minpts did not find " + targetDataset + " tcount " + str(defaultThreadCount) + " eps " + str(defaultParam[targetDataname][0]) + " minpts " + str(val))
    else:
      AggregateTest(mysnapshot, "minpts", found)
  plotsData["minptsPlots"][targetDataname] = CleanupSnapshot(mysnapshot)

  mysnapshot = CreateTestSnapshot(targetDataset, "rho", config)
  for val in defaultRhoRange:
    found =FindTestRho(serialalgo, ndalgo, algo, tests, targetDataname, defaultParam[targetDataname][0], defaultParam[targetDataname][1], val, defaultThreadCount)
    if found==None:
      print("warning: rho did not find " + targetDataset + " tcount " + str(defaultThreadCount) + " eps " + str(defaultParam[targetDataname][0]) + " minpts " + str(defaultParam[targetDataname][1]) + " rho " + str(val))
    else:
      AggregateTest(mysnapshot, "rho", found)
  plotsData["rhoPlots"][targetDataname] = CleanupSnapshot(mysnapshot)

# parse by dataset size, slightly different from above
dimNames = set()
for targetDataname in dataFiles:
  if config["ndalgo"] == "2d" and "2D" not in targetDataname:
    continue
  toRemove = targetDataname.split('_')[-1]
  dimName = targetDataname.replace(toRemove, "")
  dimNames.add(dimName)

defaultDataSizes = [100000,500000,1000000,2000000,5000000,10000000]
defaultDataSizeStr = {"100000":"100K","500000":"500K","1000000":"1M","2000000":"2M","5000000":"5M","10000000":"10M"}

if "testnumpts" in config:
  plotsData["numptsPlots"] = dict()
  for dimName in dimNames:
    mysnapshot = CreateTestSnapshot(dimName, "numpts_attr", config)
    for dataSize in defaultDataSizes:
      targetDataname = dimName+defaultDataSizeStr[str(dataSize)]
      found =FindTestByDimName(serialalgo, ndalgo, algo, tests, dimName, defaultParam[targetDataname][0], defaultParam[targetDataname][1], defaultThreadCount, dataSize)
      if found==None:
        print("warning: datasize did not find " + dimName + " tcount " + str(defaultThreadCount) + " eps " + str(defaultParam[targetDataname][0]) + " minpts " + str(defaultParam[targetDataname][1]) + " size " + str(dataSize))
      else:
        print("ok " + dimName + str(found.m_numpts))
        AggregateTest(mysnapshot, "numpts_attr", found)

    plotsData["numptsPlots"][dimName] = CleanupSnapshot(mysnapshot)

print("#numprocsPlots " + str(len(plotsData["numprocsPlots"])))
print("#epsPlots " + str(len(plotsData["epsPlots"])))
print("#minptsPlots " + str(len(plotsData["minptsPlots"])))
print("#numptsPlots " + str(len(plotsData["numptsPlots"])))
plotsPretty = pprint.pformat(plotsData)
#print(plotsPretty)

fout = open("reproduced-dbscan.py", "w")
fout.write("plotData = " + str(plotsPretty))
fout.close()
