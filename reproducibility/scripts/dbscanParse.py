from tokenDefs import *
import tokenDefs
import glob

class DTest:

  m_assigngridtime = None
  m_markcoretime = None
  m_buildcoreptquadtree = None
  m_buildgetnbrkdt = None
  m_clustercoretime = None
  m_assignbordertime = None
  m_totaltime = None
  m_numgrids = None
  m_numcore = None
  m_numborder = None
  m_numnoise = None
  m_numclusterpts = None
  m_numpts = None
  m_numprocs = None
  m_numclusters = None
  m_serialalgo = None
  m_ndalgo = None
  m_algo = None
  m_datafile = None
  m_eps = None
  m_rounds = None
  m_minpts = None
  m_rho = None
  m_dim = None

  m_timer1 = None
  m_timer2 = None
  m_timer3 = None
  m_timer4 = None
  m_timer5 = None
  m_timer6 = None
  m_timer7 = None
  m_timer8 = None
  m_timer9 = None
  m_timer10 = None


  def __init__(self, paragraph="", eps=-1, minpts=-1, totaltime=-1, numprocs=-1, datafile="", program="", numpts=-1, dimName=""):
    if program == "hpdbscan":
      self.m_serialalgo = "parallel"
      self.m_ndalgo = "nd"
      self.m_algo = "hpdbscan"
      self.m_eps = eps
      self.m_minpts = minpts
      self.m_totaltime = totaltime
      self.m_numprocs = numprocs
      self.m_datafile = datafile
      self.m_numpts = numpts
      self.m_dimName = dimName
      return
    
    if program == "pdsdbscan":
      self.m_serialalgo = "parallel"
      self.m_ndalgo = "nd"
      self.m_algo = "pdsdbscan"
      self.m_eps = eps
      self.m_minpts = minpts
      self.m_totaltime = totaltime
      self.m_numprocs = numprocs
      self.m_datafile = datafile
      self.m_numpts = numpts
      self.m_dimName = dimName
      return
    
    self.m_assigngridtime = ExtractNumber(paragraph, tokenDefs.assigngridtime_time)
    self.m_markcoretime = ExtractNumber(paragraph, tokenDefs.markcore_time)
    self.m_buildcoreptquadtree = ExtractNumber(paragraph, tokenDefs.buildcoreptquadtree_time)
    self.m_buildgetnbrkdt = ExtractNumber(paragraph, tokenDefs.buildgetnbrkdt_time)
    self.m_clustercoretime = ExtractNumber(paragraph, tokenDefs.clustercore_time)
    self.m_assignbordertime = ExtractNumber(paragraph, tokenDefs.assignborder_time)
    self.m_cleanuptime = ExtractNumber(paragraph, tokenDefs.cleanup_time)
    self.m_totaltime = ExtractNumber(paragraph, tokenDefs.total_time)
    self.m_numgrids = ExtractNumber(paragraph, tokenDefs.num_grids)
    self.m_numcore = ExtractNumber(paragraph, tokenDefs.num_core)
    self.m_numborder = ExtractNumber(paragraph, tokenDefs.num_border)
    self.m_numclusterpts = ExtractNumber(paragraph, tokenDefs.num_clusterpts)
    self.m_numpts = ExtractNumber(paragraph, tokenDefs.num_pts)
    self.m_numprocs = ExtractNumber(paragraph, tokenDefs.num_procs)
    self.m_numclusters = ExtractNumber(paragraph, tokenDefs.num_clusters)
    self.m_numnoise = ExtractNumber(paragraph, tokenDefs.num_noise)
    self.m_serialalgo = ExtractString(paragraph, tokenDefs.serial_algo)
    self.m_ndalgo = ExtractString(paragraph, tokenDefs.nd_algo)
    self.m_algo = ExtractString(paragraph, tokenDefs.algo)
    self.m_datafile = ExtractString(paragraph, tokenDefs.param_datafile)
    self.m_eps = ExtractNumber(paragraph, tokenDefs.param_eps)
    self.m_rounds = ExtractNumber(paragraph, tokenDefs.param_rounds)
    self.m_minpts = ExtractNumber(paragraph, tokenDefs.param_minpts)
    self.m_rho = ExtractNumber(paragraph, tokenDefs.param_rho)
    self.m_dim = ExtractNumber(paragraph, tokenDefs.param_dim)
    self.m_timer1 = ExtractNumber(paragraph, tokenDefs.timer1)
    self.m_timer2 = ExtractNumber(paragraph, tokenDefs.timer2)
    self.m_timer3 = ExtractNumber(paragraph, tokenDefs.timer3)
    self.m_timer4 = ExtractNumber(paragraph, tokenDefs.timer4)
    self.m_timer5 = ExtractNumber(paragraph, tokenDefs.timer5)
    self.m_timer6 = ExtractNumber(paragraph, tokenDefs.timer6)
    self.m_timer7 = ExtractNumber(paragraph, tokenDefs.timer7)
    self.m_timer8 = ExtractNumber(paragraph, tokenDefs.timer8)
    self.m_timer9 = ExtractNumber(paragraph, tokenDefs.timer9)
    self.m_timer10 = ExtractNumber(paragraph, tokenDefs.timer10)

  def GetMySpeedUp(self, other):
    # assigngrid
    # (kdtree)
    # markcore
    # (quadtree)
    # clustercore
    # assignborder
    # total
    nontimer = [other.m_assigngridtime / self.m_assigngridtime, other.m_buildgetnbrkdt / self.m_buildgetnbrkdt, other.m_markcoretime / self.m_markcoretime, other.m_buildcoreptquadtree / self.m_buildcoreptquadtree, other.m_clustercoretime / self.m_clustercoretime, other.m_assignbordertime / self.m_assignbordertime, other.m_cleanuptime / self.m_cleanuptime, other.m_totaltime / self.m_totaltime]
    timer = [other.m_timer1 / self.m_timer1, other.m_timer2 / self.m_timer2, other.m_timer3 / self.m_timer3, other.m_timer4 / self.m_timer4, other.m_timer5 / self.m_timer5, other.m_timer6 / self.m_timer6, other.m_timer7 / self.m_timer8, other.m_timer8 / self.m_timer8, other.m_timer9 / self.m_timer9, other.m_timer10 / self.m_timer10]
    return timer + nontimer
    
  def IsSingleThread(self):
    return self.m_numprocs <= 1

  def IsSerial(self):
    return self.m_serialalgo == "serial"

  def GetMyTiming(self):
    return [self.m_timer1, self.m_timer2, self.m_timer3, self.m_timer4, self.m_timer5, self.m_timer6, self.m_timer7, self.m_timer8, self.m_timer9, self.m_timer10, self.m_assigngridtime, self.m_buildgetnbrkdt, self.m_markcoretime, self.m_buildcoreptquadtree, self.m_clustercoretime, self.m_assignbordertime, self.m_cleanuptime, self.m_totaltime]

  def GetDataset(self):
    return self.m_datafile

  def GetDatasetName(self):
    datafilename = self.m_datafile.split('/')
    return datafilename[-1]

  # dimension + name, no pts count
  def GetDimName(self):
    targetDataname = self.GetDatasetName()
    toRemove = targetDataname.split('_')[-1]
    dimName = targetDataname.replace(toRemove, "")
    return dimName
  
  def GetResultInfo(self):
    return "{:1.0f} pts, {:1.0f} grids, {:1.0f} core, {:1.0f} border, {:1.0f} noise, {:1.0f} clusters, {:1.0f} clustered pts".format(self.m_numpts, self.m_numgrids, self.m_numcore, self.m_numborder, self.m_numnoise, self.m_numclusters, self.m_numclusterpts)

  def GetMyInfo(self):
    return "{serialalgo:}-{algo:}-{ndalgo:} {numprocs:1.0f}-threaded {rounds:1.0f}-rounds [eps={eps:}, minpts={minpts}, rho={rho}]".format(algo=self.m_algo, serialalgo=self.m_serialalgo, ndalgo=self.m_ndalgo, numprocs=self.m_numprocs, eps=self.m_eps, minpts=self.m_minpts, rho=self.m_rho, rounds=self.m_rounds)

  def GetMyBasicInfo(self):
    return "{serialalgo:}-{algo:}-{ndalgo:} {numprocs:1.0f}-threaded [eps={eps:}, minpts={minpts}]".format(algo=self.m_algo, serialalgo=self.m_serialalgo, ndalgo=self.m_ndalgo, numprocs=self.m_numprocs, eps=self.m_eps, minpts=self.m_minpts)

  def StringMySpeedUp2(self, other):
    # self will be run2, the faster one
    outputStr = ""
    dataset1 = other.GetDataset()
    dataset2 = self.GetDataset()
    if dataset1 == dataset2:
      outputStr += "dataset: " + dataset1 + "\n"
    else:
      outputStr += "dataset1: " + dataset1 + "\n"
      outputStr += "dataset2: " + dataset2 + "\n"

    outputStr += "result1: " + other.GetResultInfo() + "\n"
    outputStr += "result2: " + self.GetResultInfo() + "\n"

    if (other.GetResultInfo() != self.GetResultInfo()):
      outputStr += "warning: results differ\n"
    
    outputStr += "run1: " + other.GetMyInfo() + "\n"
    outputStr += "run2: " + self.GetMyInfo() + "\n"

    speedups = self.GetMySpeedUp(other)
    mytiming = self.GetMyTiming()
    mypercent = [100*i / mytiming[17] for i in mytiming]
    keywords = ["timer1", "timer2", "timer3", "timer4", "timer5", "timer6", "timer7", "timer8", "timer9", "timer10", "assign-grid: ", "kdtree-build", "mark-core: ", "quadtree-build: ", "cluster-core: ", "assign-border: ", "clean-up", "total: "]
    guyshetoldmenottoworryabout = other.GetMyTiming()
    outputStr += "{0:20} {1:10} {2:10} {3:9} {4:10}\n".format("<stage>", "<timing1>", "<timing2>", "<percent2>", "<speedup2>")
    for i, speedup in enumerate(speedups):
      if mytiming[i] < 0:
        continue
      outputStr += "{0:20} {1:<10.6f} {2:<10.6f} {3:>7.2f}%     {4:>5.2f}x\n".format(keywords[i], guyshetoldmenottoworryabout[i], mytiming[i], mypercent[i], speedup)
    return outputStr

def ParseTestFile(fileName):
  f1 = open(fileName)
  data1 = f1.read()
  splittedString = test_split.split(data1)
  tests = list()
  for paragraph in splittedString:
    finalString = ""
    if (ExtractString(paragraph, tokenDefs.param_datafile) != ""): # if is a valid test
      splittedStringRounds = round_split.split(paragraph)
      fastest = -1
      minTiming = -1
      for i,subpara in enumerate(splittedStringRounds):
        if (ExtractString(subpara, tokenDefs.param_datafile) != ""): # if is a valid header
          finalString += subpara
          continue
        timing = ExtractNumber(subpara, total_time)
        if timing >= 0: # if is a valid round
          if fastest == -1:
            fastest = i
            minTiming = timing
          else:
            if timing < minTiming:
              fastest = i
              minTiming = timing
          continue
        if (ExtractString(subpara, tokenDefs.num_clusters) != -1): # if is a valid clusterinfo
          finalString += subpara
          continue
      finalString += splittedStringRounds[fastest]
      tests.append(DTest(paragraph=finalString))

  print("[ParseTestFile] " + fileName + " parsed " + str(len(tests)) + " tests")
  return tests

def CompareTwo():

  textFiles = glob.glob("./*.txt")
  print("Found raw output files: ", textFiles)
  if len(textFiles) != 2:
    print("Test outputs not found")
    exit(1)

  print(textFiles)
  testList = list()
  for textFile in textFiles:
    testList.append(ParseTestFile(textFile))

  if testList[0][0].IsSerial():
    fastTests = testList[1]
    slowTests = testList[0]
  elif testList[1][0].IsSerial():
    fastTests = testList[0]
    slowTests = testList[1]
  elif testList[0][0].IsSingleThread():
    fastTests = testList[1]
    slowTests = testList[0]
  else:
    fastTests = testList[0]
    slowTests = testList[1]

  output = ""
  count = 0

  # by component
  speedupVec = [0.0] * 18
  avgTimeVec = [0.0] * 18
  timePercentVec = [0.0] * 18

  # by dataset
  datasetBest = dict()
  datasetWorst = dict()
  datasetSum = dict()
  datasetCount = dict()
  datasetTime = dict()

  for i,fastTest in enumerate(fastTests):
    if "GaussianDisc" in fastTest.m_datafile: # or "UniformFill" in fastTest.m_datafile:
      continue
    output += fastTest.StringMySpeedUp2(slowTests[i])
    output += "\n"

    # by component
    tmp1 = fastTest.GetMySpeedUp(slowTests[i])
    time1 = fastTest.GetMyTiming()
    timePercentVec = [s + time1[i]/time1[17] for i,s in enumerate(timePercentVec)]
    avgTimeVec = [s + time1[i] for i,s in enumerate(avgTimeVec)]
    speedupVec = [s + tmp1[i] for i,s in enumerate(speedupVec)]

    # by dataset
    if (fastTest.m_datafile not in datasetBest) or (datasetBest[fastTest.m_datafile][17] < tmp1[17]):
      datasetBest[fastTest.m_datafile] = tmp1
    if (fastTest.m_datafile not in datasetWorst) or (datasetWorst[fastTest.m_datafile][17] > tmp1[17]):
      datasetWorst[fastTest.m_datafile] = tmp1
    if (fastTest.m_datafile not in datasetSum):
      datasetSum[fastTest.m_datafile] = 0.0
      datasetTime[fastTest.m_datafile] = 0.0
      datasetCount[fastTest.m_datafile] = 0.0
    datasetSum[fastTest.m_datafile] += tmp1[17]
    datasetTime[fastTest.m_datafile] += time1[17]
    datasetCount[fastTest.m_datafile] += 1

    count += 1.0

  output += ">>> By Dataset\n"
  worstAvg = 0.0
  bestAvg = 0.0
  datasetSortArray = []
  for key in datasetBest.keys():
    datasetSortArray.append( (datasetSum[key] / datasetCount[key], key) )
  datasetSortArray.sort(reverse = 0)
  for datasetAvg, dataset in datasetSortArray:
    output += "{0:50} {1:<10.6f} {2:>5.2f}x - {3:>5.2f}x, Avg {4:>5.2f}x\n".format(dataset, datasetTime[dataset]/datasetCount[dataset], datasetWorst[dataset][17], datasetBest[dataset][17], datasetAvg)#datasetSum[dataset] / datasetCount[dataset])
    worstAvg += datasetWorst[dataset][17]
    bestAvg += datasetBest[dataset][17]
  output += "{0:50}            {1:>5.2f}x - {2:>5.2f}x\n".format("Dataset average", worstAvg/len(datasetBest.keys()), bestAvg/len(datasetBest.keys()) )

  keywords = ["timer1", "timer2", "timer3", "timer4", "timer5", "timer6", "timer7", "timer8", "timer9", "timer10", "assign-grid: ", "kdtree-build", "mark-core: ", "quadtree-build: ", "cluster-core: ", "assign-border: ", "clean-up", "total: "]
  output  += "\n>>> By Component\n"
  output += "{0:20} {1:10} {2:10} {3:9}\n".format("<stage>", "<time>", "<percent>", "<speedup>")
  for i, speedup in enumerate(speedupVec):
    if timePercentVec[i] < 0:
      continue
    output += "{0:20} {1:<10.6f} {2:>6.2f}%      {3:>5.2f}x\n".format(keywords[i], avgTimeVec[i] / count, 100*timePercentVec[i] / count, speedup / count)

  outFile = open("paramfind.out", "w")
  outFile.write(output)
  outFile.close()

  print(output)

if __name__ == "__main__":
  CompareTwo()
