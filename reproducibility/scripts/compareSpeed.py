from genTest import *
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.transforms import Bbox

# synthetic dataset testing
import synt1_2dboxbcp_new as synt1_2dboxbcp
import synt1_2dboxdt_mod as synt1_2dboxdt
import synt1_2dboxusec_new as synt1_2dboxusec
import synt1_2dgridbcp_new as synt1_2dgridbcp
import synt1_2dgriddt_mod as synt1_2dgriddt
import synt1_2dgridusec_new as synt1_2dgridusec
import synt1_ndexactbcp
import synt1_ndexactquadtree
import synt1_ndapprox_long as synt1_ndapprox
import synt1_ndapprox_long_rho
import synt2_ndapprox_quadtree as synt1_ndapproxquadtree # includes rho
import synt1_gan # also includes real
import synt1_gan_approx # also includes real
import synt1_hpdbscan_all as synt1_hpdbscan
import synt1_pdsdbscan_all as synt1_pdsdbscan

import synt1_23uf_serialbucketing_new as synt1_23uf_serialbucketing # new version does not use jemalloc
import synt3_ndexact_bucketing_all as synt3_ndexact_bucketing
import synt4_quadtreemc_quadraticbcp_all as synt3_quadtreemc_quadraticbcp
import synt4_quadtreemc_quadraticbcp_bucketing as synt_quadtreemc_quadraticbcp_bucketing
import synt4_scanmc_quadtreebcp_all as synt3_scanmc_quadtreebcp

# _mod suffix: modified description in py files
# _new suffix: new testing after 2d changes to make sure 2d is at least as fast as nd

# small real dataset testing
import real2_ndapprox_hh as real1_ndapprox
import real2_ndapprox_hh_quadtree as real1_ndapproxquadtree
import real2_ndexactbcp_hh as real1_ndexactbcp
import real2_ndexactbcpbucketing_hh as real1_ndexactbcpbucketing
import real2_ndexactquadtree_hh as real1_ndexactquadtree
import real4_quadtreemc_quadraticbcp_all as real_quadtreemc_quadraticbcp
import real4_scanmc_quadtreebcp_all as real_scanmc_quadtreebcp
import real4_quadtreemc_quadraticbcp_bucketing as real_quadtreemc_quadraticbcp_bucketing

import real1_pdsdbscan as real_pdsdbscan
# none of hpdbscan completed on two small real datasets

syntDataFiles = list()
dims = ["3D", "5D", "7D"]
dataNames = ["VisualSim","VisualVar","UniformFill"]
for dataName in dataNames:
  for dim in dims:
    syntDataFiles.append(dim + "_" + dataName + "_10M")
realDataFiles = ["7D_HouseHold_2M", "3D_GeoLife_24M"]

synt_our_exact = synt1_ndexactbcp.plotData
synt_our_exact_bucket = synt3_ndexact_bucketing.plotData
synt_our_exact_qt = synt3_quadtreemc_quadraticbcp.plotData
synt_our_exact_qt_bucket = synt_quadtreemc_quadraticbcp_bucketing.plotData
synt_our_approx = synt1_ndapprox.plotData
synt_our_approx_qt = synt1_ndapproxquadtree.plotData
synt_hpdbscan = synt1_hpdbscan.plotData
synt_pdsdbscan = synt1_pdsdbscan.plotData
#old
synt_our_exact_all_qt = synt1_ndexactquadtree.plotData

real_our_exact = real1_ndexactbcp.plotData
real_our_exact_bucket = real1_ndexactbcpbucketing.plotData
real_our_exact_qt = real_quadtreemc_quadraticbcp.plotData
real_our_exact_qt_bucket = real_quadtreemc_quadraticbcp_bucketing.plotData
real_our_approx = real1_ndapprox.plotData
real_our_approx_qt = real1_ndapproxquadtree.plotData
real_hpdbscan = None # all failed
real_pdsdbscan = real_pdsdbscan.plotData
#old
real_our_exact_all_qt = real1_ndexactquadtree.plotData

synt_real_gan = synt1_gan.plotData
synt_real_gan_approx = synt1_gan_approx.plotData

# new data on 0221
import synt4_approx_bucketing as synt_approx_bucketing
import synt4_approx_qt_bucketing as synt_approx_qt_bucketing
import real4_approx_bucketing as real_approx_bucketing
import real4_approx_qt_bucketing as real_approx_qt_bucketing

synt_our_approx_bucketing = synt_approx_bucketing.plotData
synt_our_approx_qt_bucketing = synt_approx_qt_bucketing.plotData
real_our_approx_bucketing = real_approx_bucketing.plotData
real_our_approx_qt_bucketing = real_approx_qt_bucketing.plotData

def GetParallelRange(method, dataset):
  best = 10000000
  worst = -1
  bestStr = ""
  if dataset in method["numprocsPlots"] and len(method["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = method["numprocsPlots"][dataset]["totaltime"]
    for timing in timings:
      if timing < best:
        best = timing
      if timing > worst:
        worst = timing
  return best, worst

def GetFastestSerial(dataset):
  best = 10000000
  if dataset in synt_real_gan:
    timing = synt_real_gan[dataset]["totaltime"]
    if timing < best:
      best = timing
      methodStr = "gan-tao-v2"
  if dataset in synt_our_exact["numprocsPlots"] and len(synt_our_exact["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact"
  if dataset in synt_our_exact_bucket["numprocsPlots"] and len(synt_our_exact_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-bucket"
  if dataset in synt_our_exact_qt["numprocsPlots"] and len(synt_our_exact_qt["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact_qt["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact_qt["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt"
  if dataset in synt_our_exact_qt_bucket["numprocsPlots"] and len(synt_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact_qt_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt-bucket"
  if dataset in synt_hpdbscan["numprocsPlots"] and len(synt_hpdbscan["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_hpdbscan["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_hpdbscan["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "hpdbscan"
  if dataset in synt_pdsdbscan["numprocsPlots"] and len(synt_pdsdbscan["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_pdsdbscan["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_pdsdbscan["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "pdsdbscan"

  if dataset in real_our_exact["numprocsPlots"] and len(real_our_exact["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact"
  if dataset in real_our_exact_bucket["numprocsPlots"] and len(real_our_exact_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-bucket"
  if dataset in real_our_exact_qt["numprocsPlots"] and len(real_our_exact_qt["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact_qt["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact_qt["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt"
  if dataset in real_our_exact_qt_bucket["numprocsPlots"] and len(real_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact_qt_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt-bucket"
  if dataset in real_pdsdbscan["numprocsPlots"] and len(real_pdsdbscan["numprocsPlots"][dataset]["totaltime"]) > 0 and real_pdsdbscan["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_pdsdbscan["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "pdsdbscan"
  return best, methodStr

def GetFastestParallelBaseline(dataset):
  best = 100000.0
  methodStr = "noData"
  if dataset in synt_hpdbscan["numprocsPlots"] and len(synt_hpdbscan["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = synt_hpdbscan["numprocsPlots"][dataset]["totaltime"]
    for timing in timings:
      if timing < best:
        best = timing
        methodStr = "hpdbscan"
  if dataset in synt_pdsdbscan["numprocsPlots"] and len(synt_pdsdbscan["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = synt_pdsdbscan["numprocsPlots"][dataset]["totaltime"]
    for timing in timings:
      if timing < best:
        best = timing
        methodStr = "pdsdbscan"
  if dataset in real_pdsdbscan["numprocsPlots"] and len(real_pdsdbscan["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = real_pdsdbscan["numprocsPlots"][dataset]["totaltime"]
    for timing in timings:
      if timing < best:
        best = timing
        methodStr = "pdsdbscan"
  if methodStr == "noData":
    return -1, methodStr
  else:
    return best, methodStr

# on all cores
def GetSlowestParallelBaseline(dataset):
  worst = -1
  methodStr = "noData"
  if dataset in synt_hpdbscan["numprocsPlots"] and len(synt_hpdbscan["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = synt_hpdbscan["numprocsPlots"][dataset]["totaltime"]
    cores = synt_hpdbscan["numprocsPlots"][dataset]["x"]
    if len(cores)>0 and cores[-1] == 72:
      if timings[-1] > worst:
        worst = timings[-1]
        methodStr = "hpdbscan"
  if dataset in synt_pdsdbscan["numprocsPlots"] and len(synt_pdsdbscan["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = synt_pdsdbscan["numprocsPlots"][dataset]["totaltime"]
    cores = synt_pdsdbscan["numprocsPlots"][dataset]["x"]
    if len(cores)>0 and cores[-1] == 72:
      if timings[-1] > worst:
        worst = timings[-1]
        methodStr = "pdsdbscan"
  if dataset in real_pdsdbscan["numprocsPlots"] and len(real_pdsdbscan["numprocsPlots"][dataset]["totaltime"]) > 0:
    timings = real_pdsdbscan["numprocsPlots"][dataset]["totaltime"]
    cores = real_pdsdbscan["numprocsPlots"][dataset]["x"]
    if len(cores)>0 and cores[-1] == 72:
      if timings[-1] > worst:
        worst = timings[-1]
        methodStr = "pdsdbscan"
  if methodStr == "noData":
    return -1, methodStr
  else:
    return worst, methodStr

def Analyze(syntMethod, realMethod):
  selfRelativeLo = 100000
  selfRelativeHi = -1
  selfRelativeTotal = 0
  selfRelativeCount = 0
  overBestLo = 100000
  overBestHi = -1
  overBestTotal = 0
  overBestCount = 0
  overParLo = 100000
  overParHi = -1
  for dataset in syntDataFiles:
    # print("dataset",dataset)
    bestSerial,name = GetFastestSerial(dataset)
    # print("our best serial", name, bestSerial)
    bestParBaseline,name = GetSlowestParallelBaseline(dataset)
    #print("slowest par", bestParBaseline, name)
    bestPar, worstPar = GetParallelRange(syntMethod, dataset)
    # print("our par range", bestPar, worstPar)
    selfRelative = worstPar / bestPar
    if selfRelative < selfRelativeLo:
      selfRelativeLo = selfRelative
    if selfRelative > selfRelativeHi:
      selfRelativeHi = selfRelative
    selfRelativeTotal += selfRelative
    selfRelativeCount += 1
    overBest = bestSerial / bestPar
    if overBest < overBestLo:
      overBestLo = overBest
    if overBest > overBestHi:
      overBestHi = overBest
    overBestTotal += overBest
    overBestCount += 1
    if bestParBaseline < 0:
        continue
    overPar = bestParBaseline / bestPar
    #print("overPar",overPar)
    if overPar < overParLo:
      overParLo = overPar
    if overPar > overParHi:
      overParHi = overPar
  for dataset in realDataFiles:
    # print("dataset",dataset)
    bestSerial,name = GetFastestSerial(dataset)
    bestParBaseline,name = GetSlowestParallelBaseline(dataset)
    #print("best par baseline",bestParBaseline,name)
    bestPar, worstPar = GetParallelRange(realMethod, dataset)
    # print("our par range",bestPar, worstPar)
    selfRelative = worstPar / bestPar
    if selfRelative < selfRelativeLo:
      selfRelativeLo = selfRelative
    if selfRelative > selfRelativeHi:
      selfRelativeHi = selfRelative
    selfRelativeTotal += selfRelative
    selfRelativeCount += 1      
    overBest = bestSerial / bestPar
    if overBest < overBestLo:
      overBestLo = overBest
    if overBest > overBestHi:
      overBestHi = overBest
    overBestTotal += overBest
    overBestCount += 1      
    if bestParBaseline < 0:
        continue
    overPar = bestParBaseline / bestPar
    if overPar < overParLo:
      overParLo = overPar
    if overPar > overParHi:
      overParHi = overPar
  print("self-relative", "{:>4.2f}".format(selfRelativeLo),"{:>4.2f}".format(selfRelativeHi),\
        "{:>4.2f}".format(selfRelativeTotal/selfRelativeCount))
  print("over-best", "{:>4.2f}".format(overBestLo), "{:>4.2f}".format(overBestHi), "{:>4.2f}".format(overBestTotal/overBestCount))
  print("over-par-baseline", "{:>4.2f}".format(overParLo), "{:>4.2f}".format(overParHi))
  print("")
  
print("our-exact")
Analyze(synt_our_exact, real_our_exact)
print("our-exact-bucket")
Analyze(synt_our_exact_bucket, real_our_exact_bucket)
print("our-exact-qt")
Analyze(synt_our_exact_qt, real_our_exact_qt)
print("our-exact-qt-bucket")
Analyze(synt_our_exact_qt_bucket, real_our_exact_qt_bucket)
# print("our-exact-all-qt")
# Analyze(synt_our_exact_all_qt, real_our_exact_all_qt
print("")
print("our-approx")
Analyze(synt_our_approx, real_our_approx)
print("our-approx-bucketing")
Analyze(synt_our_approx_bucketing, real_our_approx_bucketing)
print("our-approx-qt")
Analyze(synt_our_approx_qt, real_our_approx_qt)
print("our-approx-qt-bucketing")
Analyze(synt_our_approx_qt_bucketing, real_our_approx_qt_bucketing)
print("")
print("our-exact-bucket", "geolife")
print("eps", real_our_exact_bucket["epsPlots"]["3D_GeoLife_24M"]["x"])
print("time", real_our_exact_bucket["epsPlots"]["3D_GeoLife_24M"]["totaltime"])
print("")

def GetFastestSerialBaseline(dataset):
  best = 10000000
  methodStr = "noData"
  if dataset in synt_real_gan:
    timing = synt_real_gan[dataset]["totaltime"]
    if timing < best:
      best = timing
      methodStr = "gan-tao-v2"
  if methodStr == "noData":
    return -1, methodStr
  return best, methodStr

def GetFastestSerialApproxBaseline(dataset):
  best = 10000000
  methodStr = "noData"
  if dataset in synt_real_gan_approx:
    timing = synt_real_gan_approx[dataset]["totaltime"]
    if timing < best:
      best = timing
      methodStr = "gan-tao-v2-approx"
  if methodStr == "noData":
    return -1, methodStr
  return best, methodStr

def GetOurFastestSerialExact(dataset):
  best = 10000000
  if dataset in synt_our_exact["numprocsPlots"] and len(synt_our_exact["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact"
  if dataset in synt_our_exact_bucket["numprocsPlots"] and len(synt_our_exact_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-bucket"
  if dataset in synt_our_exact_qt["numprocsPlots"] and len(synt_our_exact_qt["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact_qt["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact_qt["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt"
  if dataset in synt_our_exact_qt_bucket["numprocsPlots"] and len(synt_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_exact_qt_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt-bucket"
  if dataset in real_our_exact["numprocsPlots"] and len(real_our_exact["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact"
  if dataset in real_our_exact_bucket["numprocsPlots"] and len(real_our_exact_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-bucket"
  if dataset in real_our_exact_qt["numprocsPlots"] and len(real_our_exact_qt["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact_qt["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact_qt["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt"
  if dataset in real_our_exact_qt_bucket["numprocsPlots"] and len(real_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_exact_qt_bucket["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt-bucket"
  return best, methodStr

def GetOurFastestParallelExact(dataset):
  best = 10000000
  if dataset in synt_our_exact["numprocsPlots"] and len(synt_our_exact["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_exact["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact"
  if dataset in synt_our_exact_bucket["numprocsPlots"] and len(synt_our_exact_bucket["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_exact_bucket["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact-bucket"
  if dataset in synt_our_exact_qt["numprocsPlots"] and len(synt_our_exact_qt["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_exact_qt["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt"
  if dataset in synt_our_exact_qt_bucket["numprocsPlots"] and len(synt_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt-bucket"
  if dataset in real_our_exact["numprocsPlots"] and len(real_our_exact["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_exact["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact"
  if dataset in real_our_exact_bucket["numprocsPlots"] and len(real_our_exact_bucket["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_exact_bucket["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact-bucket"
  if dataset in real_our_exact_qt["numprocsPlots"] and len(real_our_exact_qt["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_exact_qt["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt"
  if dataset in real_our_exact_qt_bucket["numprocsPlots"] and len(real_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_exact_qt_bucket["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-exact-qt-bucket"
  return best, methodStr

def GetOurFastestSerialApprox(dataset):
  best = 10000000
  if dataset in synt_our_approx["numprocsPlots"] and len(synt_our_approx["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_approx["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_approx["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx"
  if dataset in synt_our_approx_bucketing["numprocsPlots"] and len(synt_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_approx_bucketing["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx-bucketing"
  if dataset in synt_our_approx_qt["numprocsPlots"] and len(synt_our_approx_qt["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_approx_qt["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_approx_qt["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt"
  if dataset in synt_our_approx_qt_bucketing["numprocsPlots"] and len(synt_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0 and synt_our_approx_qt_bucketing["numprocsPlots"][dataset]["x"][0] == 1:
    timing = synt_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt-bucketing"
  if dataset in real_our_approx["numprocsPlots"] and len(real_our_approx["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_approx["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_approx["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx"
  if dataset in real_our_approx_bucketing["numprocsPlots"] and len(real_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_approx_bucketing["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx-bucketing"
  if dataset in real_our_approx_qt["numprocsPlots"] and len(real_our_approx_qt["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_approx_qt["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_approx_qt["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt"
  if dataset in real_our_approx_qt_bucketing["numprocsPlots"] and len(real_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0 and real_our_approx_qt_bucketing["numprocsPlots"][dataset]["x"][0] == 1:
    timing = real_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"][0]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt-bucketing"
  return best, methodStr

def GetOurFastestParallelApprox(dataset):
  best = 10000000
  if dataset in synt_our_approx["numprocsPlots"] and len(synt_our_approx["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_approx["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx"
  if dataset in synt_our_approx_bucketing["numprocsPlots"] and len(synt_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx-bucketing"      
  if dataset in synt_our_approx_qt["numprocsPlots"] and len(synt_our_approx_qt["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_approx_qt["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt"
  if dataset in synt_our_approx_qt_bucketing["numprocsPlots"] and len(synt_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = synt_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt-bucketing"
  if dataset in real_our_approx["numprocsPlots"] and len(real_our_approx["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_approx["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx"
  if dataset in real_our_approx_bucketing["numprocsPlots"] and len(real_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_approx_bucketing["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx-bucketing"
  if dataset in real_our_approx_qt["numprocsPlots"] and len(real_our_approx_qt["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_approx_qt["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt"
  if dataset in real_our_approx_qt_bucketing["numprocsPlots"] and len(real_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"]) > 0:
    timing = real_our_approx_qt_bucketing["numprocsPlots"][dataset]["totaltime"][-1]
    if timing < best:
      best = timing
      methodStr = "our-approx-qt-bucketing"
  return best, methodStr

avgExact = 0
exactCounter = 0
for dataset in syntDataFiles+realDataFiles:
  blTime, name = GetFastestSerialBaseline(dataset)
  ourTime, nameour = GetOurFastestSerialExact(dataset)
  if blTime > 0:
    avgExact += blTime/ourTime
    exactCounter += 1
print("speedup over exact gantao-v2", avgExact/exactCounter)

avgApprox = 0
approxCounter = 0
for dataset in syntDataFiles+realDataFiles:
  blTime, name = GetFastestSerialApproxBaseline(dataset)
  ourTime, nameour = GetOurFastestSerialApprox(dataset)
  if blTime > 0:
    avgApprox += blTime/ourTime
    approxCounter += 1
print("speedup over approx gantao-v2", avgApprox/approxCounter)

avgApprox = 0
approxCounter = 0
for dataset in syntDataFiles+realDataFiles:
  ourATime, nameour = GetOurFastestSerialApprox(dataset)
  ourETime, nameour = GetOurFastestSerialExact(dataset)
  avgApprox += ourATime/ourETime
  approxCounter += 1
print("serial exact over approx", avgApprox/approxCounter)

avgApprox = 0
approxCounter = 0
for dataset in syntDataFiles+realDataFiles:
  ourATime, nameour = GetOurFastestParallelApprox(dataset)
  ourETime, nameour = GetOurFastestParallelExact(dataset)
  avgApprox += ourATime/ourETime
  approxCounter += 1
print("parallel exact over approx", avgApprox/approxCounter)

