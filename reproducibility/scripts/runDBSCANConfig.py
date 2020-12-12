import subprocess
import sys
import random
import os
import time
from os import path

def onPprocessors(command,p) :
  return "CILK_NWORKERS="+str(p)+" " + command

def shellGetOutput(commandStr) :
  process = subprocess.Popen(commandStr,shell=True,stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
  output, err = process.communicate()
  strOut = output.decode('utf-8') + '\n' + err.decode('utf-8')
  return strOut

def stripFloat(val) :
  trunc = float(int(val*1000))/1000
  return str(trunc).rstrip('0')    

def runSingle(runProgram, options, ifile, procs) :
  comString = "./"+runProgram+" "+options+" "+ifile
  myProgOutput = comString + "\n"
  currentTime = "Test start time: " + time.strftime("%Y%m%d-%H:%M:%S") + "\n"
  myProgOutput += currentTime
  print(currentTime)
  comString = onPprocessors(comString,procs)
  out = shellGetOutput(comString)
  myProgOutput += out
  return myProgOutput

def runTest(runProgram, checkProgram, dataDir, test, rounds, procs, noOutput) :
  random.seed()
  myOut = ""
  outFile="/tmp/ofile%d_%d" %(random.randint(0, 1000000), random.randint(0, 1000000)) 
  [specifiedProcs, inputFileNames, runOptions, checkOptions] = test
  if type(inputFileNames) is str :
    inputFileNames = [inputFileNames]
    shortInputNames = " ".join(inputFileNames)
    if len(dataDir)>0:
      #out = shellGetOutput("cd " + dataDir + "; make " + shortInputNames)
      longInputNames = " ".join(dataDir + "/" + name for name in inputFileNames)
      runOptions = runOptions + " -r " + str(rounds)
    if (noOutput == 0) :
      runOptions = runOptions + " -o " + outFile
    if specifiedProcs <= 0:
      finalProcs = procs
    else:
      finalProcs = specifiedProcs
    singleOut = runSingle(runProgram, runOptions, longInputNames, finalProcs)
    myOut += singleOut
    if (noOutput == 0) :
      checkString = ("./" + checkProgram + " " + checkOptions + " "
                     + longInputNames + " " + outFile)
      checkOut = shellGetOutput(checkString)
      myOut += checkOut
      if path.exists(outFile):
        os.remove(outFile)
    print(myOut)
    return myOut
  
def averageTime(times) :
    return sum(times)/len(times)

def timeAll(name, runProgram, checkProgram, dataDir, tests, rounds, procs, message, noOutput, header) :
  #out = shellGetOutput("rm ./DBSCAN")
  #out = shellGetOutput("ln -s " + runProgram + " .")
  try:
    results = list()
    i = 0
    for i,test in enumerate(tests):
      i += 1
      results.append(runTest(runProgram, checkProgram, dataDir, test, rounds, procs, noOutput))
      if i % 50 == 0:
        timestr = time.strftime("%Y%m%d-%H%M%S")
        allTestOutputs = ""
        for testOut in results:
          allTestOutputs += testOut
        myTextFile = open("./outputs/" + message + "-" + str(i) + "-" + timestr + ".txt", "a")
        myTextFile.write(message + "\n" + header + "\n" + allTestOutputs)
        myTextFile.close()
        results = list()

    timestr = time.strftime("%Y%m%d-%H%M%S")
    allTestOutputs = ""
    for testOut in results:
      allTestOutputs += testOut
    myTextFile = open("./outputs/" + message + "-" + str(i) + "-" + timestr + ".txt", "a")
    myTextFile.write(message + "\n" + header + "\n" + allTestOutputs)
    myTextFile.close()
  except KeyboardInterrupt:
    return 1

def getOption(str) :
  a = sys.argv
  l = len(a)
  for i in range(1, l) :
    if (a[i] == str) :
      return True
  return False

def getArg(str, default) :
  a = sys.argv
  l = len(a)
  for i in range(1, l) :
    if (a[i] == str and  (i+1 != l)) :
        return sys.argv[i+1]
  return default

def getArgs() :
  noOutput = getOption("-x")
  processors = int(getArg("-p", 0)) # needs more fine grained control than this
  rounds = int(getArg("-r", 1))
  message = getArg("-m", "")
  return (noOutput, rounds, processors, message)

def timeAllArgs(runProgram, checkProgram, dataDir, tests, header) :
  (noOutput, rounds, procs, message) = getArgs()
  name = os.path.basename(os.getcwd())
  timeAll(name, runProgram, checkProgram, dataDir, tests, rounds, procs, message, noOutput, header)
