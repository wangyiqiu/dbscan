from genTest import *
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['text.usetex'] = True


from matplotlib.transforms import Bbox

# synthetic dataset testing
import synt1_2dboxbcp_new as synt1_2dboxbcp
import synt1_2dboxdt_mod as synt1_2dboxdt
import synt1_2dboxusec_new as synt1_2dboxusec
import synt1_2dgridbcp_new as synt1_2dgridbcp
import synt1_2dgriddt_mod as synt1_2dgriddt
import synt1_2dgridusec_new as synt1_2dgridusec
import synt1_ndexactbcp
# import synt1_ndexactquadtree
import synt1_ndapprox_long as synt1_ndapprox
import synt1_ndapprox_long_rho
import synt2_ndapprox_quadtree as synt1_ndapproxquadtree # includes rho
import synt1_gan # also includes real
import synt1_hpdbscan_all as synt1_hpdbscan
import synt1_pdsdbscan_all as synt1_pdsdbscan

import synt1_23uf_serialbucketing_new as synt1_23uf_serialbucketing # new version does not use jemalloc
import synt3_ndexact_bucketing_all as synt3_ndexact_bucketing
# import synt4_quadtreemc_quadraticbcp_all as synt3_quadtreemc_quadraticbcp
import synt4_quadtreemc_quadraticbcp_all as synt1_ndexactquadtree
import synt4_quadtreemc_quadraticbcp_bucketing as synt_quadtreemc_quadraticbcp_bucketing
import synt4_scanmc_quadtreebcp_all as synt3_scanmc_quadtreebcp
import synt4_approx_bucketing as synt_approx_bucketing
import synt4_approx_qt_bucketing as synt_approx_qt_bucketing

# _mod suffix: modified description in py files
# _new suffix: new testing after 2d changes to make sure 2d is at least as fast as nd

# small real dataset testing
import real2_ndapprox_hh as real1_ndapprox
import real2_ndapprox_hh_quadtree as real1_ndapproxquadtree
import real2_ndexactbcp_hh as real1_ndexactbcp
import real2_ndexactbcpbucketing_hh as real1_ndexactbcpbucketing
# import real2_ndexactquadtree_hh as real1_ndexactquadtree
import real4_quadtreemc_quadraticbcp_all as real1_ndexactquadtree
# import real4_quadtreemc_quadraticbcp_all as real_quadtreemc_quadraticbcp
import real4_scanmc_quadtreebcp_all as real_scanmc_quadtreebcp
import real4_quadtreemc_quadraticbcp_bucketing as real_quadtreemc_quadraticbcp_bucketing

import real1_pdsdbscan as real_pdsdbscan
import real4_approx_bucketing as real_approx_bucketing
import real4_approx_qt_bucketing as real_approx_qt_bucketing

# none of hpdbscan completed on two small real datasets

plotDir = "./plots/"
speedupDir = "./plots/speedup/"
epsDir = "./plots/eps/"
minptsDir = "./plots/minpts/"
numptsDir = "./plots/numpts/"

plotIndividual=False

plot2D = True
if False:
  plotExactBcp = False
  plotExactQuadtree = False
  plotApprox = False
  plotGan = False
  plotHp = False
  plotPds = False
  plotSBucketing = False
else:
  plotExactBcp = True
  plotExactQuadtree = True
  plotApprox = True
  plotGan = True
  plotHp = True
  plotPds = True
  plotSBucketing = True

saveFileName = "fullplots"

datasetOrder = list()
dims = ["2D", "3D", "5D", "7D"]
dataNames = ["VisualSim","VisualVar","UniformDisc","GaussianDisc","UniformFill","UniformStrip","GaussianStrip"]
plotLetters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"]
def GetLetter(idx):
  return r"$\bf{(" + plotLetters[idx] + ")}$ "

def FixName(myname):
  myname = myname.replace("VisualSim","SS-simden")
  myname = myname.replace("VisualVar","SS-varden")
  myname = myname.replace("7D_HouseHold_2M","7D_Household_2.05M")
  myname = myname.replace("3D_GeoLife_24M","3D_GeoLife_24.9M")
  myname = myname.replace("_","-")
  return myname

realDataFiles = ["7D_HouseHold_2M", "3D_GeoLife_24M"]

def GetName(name):
  # todo
  return name

for dataName in dataNames:
  if "Strip" in dataName:
    datasetOrder.append("2D_" + dataName + "_")
  else:
    for dim in dims:
      datasetOrder.append(dim + "_" + dataName + "_")

datasetOrder2D = list()
for dataName in datasetOrder:
  if "2D" in dataName:
    datasetOrder2D.append(dataName)

matplotlib.rcParams.update({'font.size': 6})

def PlotDecomp(ax, data):
    axs[i,j].plot(data["x"], data["totaltime"], c='k', marker='^', mfc="None",linewidth=1, label="total")
    axs[i,j].plot(data["x"], data["assigngridtime"], c='k', marker='o', mfc="None", ms=4, linewidth=0.3,linestyle="--", label="build-grid")
    axs[i,j].plot(data["x"], data["markcoretime"], c='k', marker='v', mfc="None", ms=4, linewidth=0.3,linestyle="--", label="mark-core")
    axs[i,j].plot(data["x"], data["clustercoretime"], c='k', marker='s', mfc="None", ms=4, linewidth=0.3,linestyle="--", label="cluster-core")
    axs[i,j].plot(data["x"], data["assignbordertime"], c='k', marker='X', mfc="None", ms=4, linewidth=0.3,linestyle="--", label="assign-border")

def Plot2D(axs, i, j, toPlot, plotType, plotFunc):
    my2dmarkers = ['s','^','o','s','^','o']
    my2dcolors = ['b','r','b',"tab:brown","tab:pink","tab:gray",'r',"tab:olive","tab:purple","tab:cyan"]
    if "2D" in toPlot:
      data = synt1_2dgriddt.plotData[plotType][toPlot]
      plotFunc(data["x"], data["totaltime"], c=my2dcolors[0], marker=my2dmarkers[1],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="3our-2d-grid-delaunay", linestyle='-',mfc="None")
      data = synt1_2dgridusec.plotData[plotType][toPlot]
      plotFunc(data["x"], data["totaltime"], c=my2dcolors[0], marker=my2dmarkers[2],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="2our-2d-grid-usec", linestyle='-',mfc="None")
      data = synt1_2dgridbcp.plotData[plotType][toPlot]
      plotFunc(data["x"], data["totaltime"], c=my2dcolors[0], marker=my2dmarkers[0],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="1our-2d-grid-bcp", linestyle='-',mfc="None")
      data = synt1_2dboxdt.plotData[plotType][toPlot]
      plotFunc(data["x"], data["totaltime"], c=my2dcolors[1], marker=my2dmarkers[4],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="6our-2d-box-delaunay", linestyle='-.',mfc="None")
      data = synt1_2dboxusec.plotData[plotType][toPlot]
      plotFunc(data["x"], data["totaltime"], c=my2dcolors[1], marker=my2dmarkers[5],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="5our-2d-box-usec", linestyle='-.',mfc="None")
      data = synt1_2dboxbcp.plotData[plotType][toPlot]
      plotFunc(data["x"], data["totaltime"], c=my2dcolors[1], marker=my2dmarkers[3],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="4our-2d-box-bcp", linestyle='-.',mfc="None")

def Plot2DCores(axs, i, j, toPlot):
    my2dmarkers = ['.','.','.','.','.','.']
    my2dmarkers = ['s','^','o','s','^','o']
    my2dcolors = ['b','r','b',"tab:brown","tab:pink","tab:gray",'r',"tab:olive","tab:purple","tab:cyan"]
    if "2D" in toPlot:
      data = synt1_2dgriddt.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [i/2 for i in myx]
      axs[i,j].semilogy(myx, data["totaltime"], c=my2dcolors[0], marker=my2dmarkers[1],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="3our-2d-grid-delaunay", linestyle='-',mfc="None")
      data = synt1_2dgridusec.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [i/2 for i in myx]
      axs[i,j].semilogy(myx, data["totaltime"], c=my2dcolors[0], marker=my2dmarkers[2],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="2our-2d-grid-usec", linestyle='-',mfc="None")
      data = synt1_2dgridbcp.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [i/2 for i in myx]
      axs[i,j].semilogy(myx, data["totaltime"], c=my2dcolors[0], marker=my2dmarkers[0],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="1our-2d-grid-bcp", linestyle='-',mfc="None")
      data = synt1_2dboxdt.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [i/2 for i in myx]
      axs[i,j].semilogy(myx, data["totaltime"], c=my2dcolors[1], marker=my2dmarkers[4],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="6our-2d-box-delaunay", linestyle='-.',mfc="None")
      data = synt1_2dboxusec.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [i/2 for i in myx]
      axs[i,j].semilogy(myx, data["totaltime"], c=my2dcolors[1], marker=my2dmarkers[5],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="5our-2d-box-usec", linestyle='-.',mfc="None")
      data = synt1_2dboxbcp.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [i/2 for i in myx]
      axs[i,j].semilogy(myx, data["totaltime"], c=my2dcolors[1], marker=my2dmarkers[3],linewidth=0.3,markeredgewidth=0.3,markersize=4, label="4our-2d-box-bcp", linestyle='-.',mfc="None")

def Plot2DSpeedup(axs, toPlot):
    my2dmarkers = ['.','.','.','.','.','.']
    my2dmarkers = ['s','^','o','s','^','o']
    my2dcolors = ['b','r','b',"tab:brown","tab:pink","tab:gray",'r',"tab:olive","tab:purple","tab:cyan"]
    if "2D" in toPlot:
      data = synt1_2dgriddt.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [math.ceil(i/2) for i in myx]
      runTimes = data["totaltime"]
      #speedups = [runTimes[0]/i for i in runTimes]
      speedups,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      axs.semilogx(myx, speedups, c=my2dcolors[0], marker=my2dmarkers[1],linewidth=0.3,markeredgewidth=0.3,markersize=4,mfc="None", label="3our-2d-grid-delaunay", linestyle='-',basex=2)
      data = synt1_2dgridusec.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [math.ceil(i/2) for i in myx]
      runTimes = data["totaltime"]
      #speedups = [runTimes[0]/i for i in runTimes]
      speedups,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      axs.semilogx(myx, speedups, c=my2dcolors[0], marker=my2dmarkers[2],linewidth=0.3,markeredgewidth=0.3,markersize=4,mfc="None", label="2our-2d-grid-usec", linestyle='-',basex=2)
      data = synt1_2dgridbcp.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [math.ceil(i/2) for i in myx]
      runTimes = data["totaltime"]
      #speedups = [runTimes[0]/i for i in runTimes]
      speedups,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      axs.semilogx(myx, speedups, c=my2dcolors[0], marker=my2dmarkers[0],linewidth=0.3,markeredgewidth=0.3,markersize=4,mfc="None", label="1our-2d-grid-bcp", linestyle='-',basex=2)

      data = synt1_2dboxdt.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [math.ceil(i/2) for i in myx]
      runTimes = data["totaltime"]
      #speedups = [runTimes[0]/i for i in runTimes]
      speedups,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      axs.semilogx(myx, speedups, c=my2dcolors[1], marker=my2dmarkers[4],linewidth=0.3,markeredgewidth=0.3,markersize=4,mfc="None", label="6our-2d-box-delaunay", linestyle='-.',basex=2)
      data = synt1_2dboxusec.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [math.ceil(i/2) for i in myx]
      runTimes = data["totaltime"]
      #speedups = [runTimes[0]/i for i in runTimes]
      speedups,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      axs.semilogx(myx, speedups, c=my2dcolors[1], marker=my2dmarkers[5],linewidth=0.3,markeredgewidth=0.3,markersize=4,mfc="None", label="5our-2d-box-usec", linestyle='-.',basex=2)
      data = synt1_2dboxbcp.plotData["numprocsPlots"][toPlot]
      myx = data["x"]
      myx = [math.ceil(i/2) for i in myx]
      runTimes = data["totaltime"]
      #speedups = [runTimes[0]/i for i in runTimes]
      speedups,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      axs.semilogx(myx, speedups, c=my2dcolors[1], marker=my2dmarkers[3],linewidth=0.3,markeredgewidth=0.3,markersize=4,mfc="None", label="4our-2d-box-bcp", linestyle='-.',basex=2)

# approx methods are not counted
def GetSpeedups(toPlot, runTimes):
  best = 10000000
  bestStr = ""
  data = synt1_ndexactbcp.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact"
  data = synt3_ndexact_bucketing.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact-bucket " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact-bucket"
  data = synt1_ndexactquadtree.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact-qt " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact-qt"
  data = synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact-qt-bucket " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact-qt-bucket"
  data = real1_ndexactbcp.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact"
  data = real1_ndexactbcpbucketing.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact-bucket"
  data = real1_ndexactquadtree.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact-qt " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact-qt"
  data = real1_ndexactquadtree.plotData["numprocsPlots"]
  if toPlot in data:
    if data[toPlot]["totaltime"][0] < best:
      best = data[toPlot]["totaltime"][0]
      bestStr = "our-exact-qt-bucket " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "our-exact-qt-bucket"
  data = synt1_gan.plotData
  if toPlot in data:
    if data[toPlot]["totaltime"] < best:
      best = data[toPlot]["totaltime"]
      bestStr = "Gan&Tao-v2 " + "({:>4.2f}".format(best) + " sec)"
      methodStr = "Gan&Tao-v2"
  speedUps = [best/i for i in runTimes]
  return speedUps, bestStr, methodStr

def TitleSuffix(eps, minpts):
  return "\n(eps:"+str(int(eps))+",minpts:"+str(int(minpts))+")"
def TitleSuffixEps(eps):
  return "\n(eps:"+str(int(eps))+")"
def TitleSuffixMinpts(minpts):
  return "\n(minpts:"+str(int(minpts))+")"

def unique_everseen(seq, key=None):
    seen = set()
    seen_add = seen.add
    return [x for x,k in zip(seq,key) if not (k in seen or seen_add(k))]

def reorderLegend(handles, labels, order=None, unique=True):
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0])) # sort both labels and handles by labels
    if order is not None: # Sort according to a given list (not necessarily complete)
      keys=dict(zip(order,range(len(order))))
      labels, handles = zip(*sorted(zip(labels, handles), key=lambda t,keys=keys: keys.get(t[0],np.inf)))
    if unique:
      labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels)) # Keep only the first of each handle
    labels = (i[1:] for i in labels)
    return(handles, labels)

def get_ax_size(ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

def set_ax_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)
  
def SaveThisSubplot(ax, name):
  from matplotlib.transforms import TransformedBbox, Affine2D
  #bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  tight_bbox_raw = ax.get_tightbbox(fig.canvas.get_renderer())
  bbox = TransformedBbox(tight_bbox_raw, Affine2D().scale(1./fig.dpi))
  # bbox=bbox.expanded(1.25,1.25)
  # bbox.x1 -= 0.145;
  # bbox.y0 += 0.03;
  print(name)
  print(get_ax_size(ax))
  fig.savefig(name + ".pdf", bbox_inches=bbox, format="pdf")

realDataFiles.reverse()
datasetOrder = []
datasetOrder = list()
dims = ["3D", "5D", "7D"]
dataNames = ["VisualSim","VisualVar","UniformFill"]
for dataName in dataNames:
  for dim in dims:
    datasetOrder.append(dim + "_" + dataName + "_")
#datasetOrder.append("3D_GaussianDisc_")
print(datasetOrder)

# higher-d plots

# SCALABILITY for self-comparison, the baselines

def GetSelfSpeedups(toPlot, runTimes):
  best = runTimes[0]
  bestStr = "self-relative speedup"
  speedUps = [best/i for i in runTimes]
  return speedUps, bestStr, None

# Self-relative speedups

toPlot = "3D_VisualVar_10M"
printDataName = GetName(toPlot)
fig, axs = plt.subplots(1, 1, constrained_layout=False, figsize =(3.8,1.6))
myax = axs
if toPlot in synt1_ndexactquadtree.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, c='b',linewidth=0.3,marker='^',label="0our-exact-qt",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
if toPlot in synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, linestyle='--', c='g',linewidth=0.3,marker='^',label="1our-exact-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
if toPlot in synt1_ndexactbcp.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, c='b',linewidth=0.3,marker='s',label="2our-exact",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
if toPlot in synt3_ndexact_bucketing.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt3_ndexact_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt3_ndexact_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, linestyle='--', c='g',linewidth=0.3,marker='s',label="3our-exact-bucketing",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)      
if toPlot in synt1_ndapproxquadtree.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt1_ndapproxquadtree.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_ndapproxquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, c='r',linewidth=0.3,marker='x',label="4our-approx-qt",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)
if toPlot in synt_approx_qt_bucketing.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt_approx_qt_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt_approx_qt_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, linestyle='--', c='k',linewidth=0.3,marker='x',label="5our-approx-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)
if toPlot in synt1_ndapprox.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt1_ndapprox.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_ndapprox.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, c='r',linewidth=0.3,marker='+',label="6our-approx",markeredgewidth=0.3,markersize=5,mfc="None", basex=2)
if toPlot in synt_approx_bucketing.plotData["numprocsPlots"]:
  numCores = [math.ceil(i/2) for i in synt_approx_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt_approx_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCores, scalability, linestyle='--', c='k',linewidth=0.3,marker='+',label="7our-approx-bucketing",markeredgewidth=0.3,markersize=5,mfc="None", basex=2)  
if toPlot in synt1_hpdbscan.plotData["numprocsPlots"]:
  numCoresTest = [math.ceil(i/2) for i in synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='^',label="8hpdbscan",markeredgewidth=0.3,markersize=4,mfc='k',basex=2)
if toPlot in synt1_pdsdbscan.plotData["numprocsPlots"]:
  numCoresTest = [math.ceil(i/2) for i in synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
  myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='s',label="9pdsdbscan",markeredgewidth=0.3,markersize=3,mfc='k',basex=2)
  
numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
scalability,baselineStr,methodStr = GetSelfSpeedups(toPlot, runTimes)
myax.set_ylabel("self-relative speedup", fontsize=9,labelpad=0)
myax.set_xlabel("num-threads", fontsize=9,labelpad=0)
myax.set_title(FixName(printDataName) + TitleSuffix(list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["eps"])[0], list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["minpts"])[0]),fontsize=9.0, fontstyle="italic")
myax.grid(ls=':',linewidth=0.3,c='k')
myax.set_xticks(numCores)
myax.set_xticklabels(["1","4","8","16","24","36","36h"])
for l in myax.yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in myax.xaxis.get_ticklabels():
  l.set_visible(True)
  l.set_fontsize(9)

handles, labels = axs.get_legend_handles_labels()
handles, labels = reorderLegend(handles, labels)
fig.legend(handles, labels, loc='upper right',ncol=1, bbox_to_anchor=(1, 1),fontsize=9,labelspacing=0.2)
plt.subplots_adjust(left=0.09,right=0.5,bottom=0.19,top=0.8,wspace=0.45,hspace=0.7)
#box = axs.get_position()
#axs.set_position([box.x0+0.01,box.y0+0.05,box.x1-box.x0,box.y1-box.y0])
plt.savefig(plotDir + "paper_nd_speedup_self.pdf", format="pdf")

# MINPTS

fig, axs = plt.subplots(3, 4, constrained_layout=False, figsize =(8.5,11/2.2))
for j in range(3):
  for i in range(3):
    accessIdx = j*3 + i
    accessLetter = i*3+j
    if accessIdx >= len(datasetOrder):
      continue
    toPlot = datasetOrder[accessIdx]
    toPlot = toPlot + "10M"
    printDataName = GetName(toPlot)
    if False:#plot2D:
      Plot2D(axs, i, j, toPlot, "minptsPlots", axs[i,j].loglog)
    if toPlot in synt1_ndexactquadtree.plotData["minptsPlots"]:
      data = synt1_ndexactquadtree.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], c='b', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="0our-exact-qt",mfc="None")
    if toPlot in synt_quadtreemc_quadraticbcp_bucketing.plotData["minptsPlots"]:
      data = synt_quadtreemc_quadraticbcp_bucketing.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='g', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="1our-exact-qt-bucketing",mfc="None")
    if toPlot in synt1_ndexactbcp.plotData["minptsPlots"]:
      data = synt1_ndexactbcp.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], c='b', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="2our-exact",mfc="None")
    if toPlot in synt3_ndexact_bucketing.plotData["minptsPlots"]:
      data = synt3_ndexact_bucketing.plotData["minptsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='g', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="3our-exact-bucketing",mfc="None")
    if toPlot in synt1_ndapproxquadtree.plotData["minptsPlots"]:
      data = synt1_ndapproxquadtree.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], c='r', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="4our-approx-qt",mfc="None")
    if toPlot in synt_approx_qt_bucketing.plotData["minptsPlots"]:
      data = synt_approx_qt_bucketing.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='k', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="5our-approx-qt-bucketing",mfc="None")
    if toPlot in synt1_ndapprox.plotData["minptsPlots"]:
      data = synt1_ndapprox.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], c='r', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="6our-approx",mfc="None")
    if toPlot in synt_approx_bucketing.plotData["minptsPlots"]:
      data = synt_approx_bucketing.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='k', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="7our-approx-bucketing",mfc="None")
    if toPlot in synt1_hpdbscan.plotData["minptsPlots"]:
      data = synt1_hpdbscan.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
    if toPlot in synt1_pdsdbscan.plotData["minptsPlots"]:
      data = synt1_pdsdbscan.plotData["minptsPlots"][toPlot]
      axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
    axs[i,j].set_ylabel("time (sec)", fontsize=9,labelpad=0)
    axs[i,j].set_xlabel("minpts", fontsize=9,labelpad=0)
    axs[i,j].set_title(GetLetter(accessLetter) + FixName(printDataName)+TitleSuffixEps(list(synt1_ndexactbcp.plotData["minptsPlots"][toPlot]["eps"])[0]),fontsize=9.0, fontstyle="italic")
    axs[i,j].grid(ls=':',linewidth=0.3,c='k')
    for l in axs[i,j].yaxis.get_ticklabels():
      l.set_fontsize(9)
    for l in axs[i,j].xaxis.get_ticklabels():
      l.set_fontsize(9)

axs[0,3].axis("off")
positions = [(1,3),(2,3)]
for idx,toPlot in enumerate(realDataFiles):
  i = positions[idx][0]
  j = positions[idx][1]
  accessIdx = j*3+i
  accessLetter = j*3+i-1
  if toPlot in real1_ndexactquadtree.plotData["minptsPlots"]:
    data = real1_ndexactquadtree.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], c='b',linewidth=0.3,marker='^',label="0our-exact-qt",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real_quadtreemc_quadraticbcp_bucketing.plotData["minptsPlots"]:
    data = real_quadtreemc_quadraticbcp_bucketing.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='g',linewidth=0.3,marker='^',label="1our-exact-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndexactbcp.plotData["minptsPlots"]:
    data = real1_ndexactbcp.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], c='b',linewidth=0.3,marker='s',label="2our-exact",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndexactbcpbucketing.plotData["minptsPlots"]:
    data = real1_ndexactbcpbucketing.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='g',linewidth=0.3,marker='s',label="3our-exact-bucketing",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndapproxquadtree.plotData["minptsPlots"]:
    data = real1_ndapproxquadtree.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], c='r',linewidth=0.3,marker='x',label="4our-approx-qt",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real_approx_qt_bucketing.plotData["minptsPlots"]:
    data = real_approx_qt_bucketing.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='k',linewidth=0.3,marker='x',label="5our-approx-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndapprox.plotData["minptsPlots"]:
    data = real1_ndapprox.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], c='r',linewidth=0.3,marker='+',label="6our-approx",markeredgewidth=0.3,markersize=5,mfc="None")
  if toPlot in real_approx_bucketing.plotData["minptsPlots"]:
    data = real_approx_bucketing.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], linestyle='--', c='k',linewidth=0.3,marker='+',label="7our-approx-bucketing",markeredgewidth=0.3,markersize=5,mfc="None")
  if toPlot in real_pdsdbscan.plotData["minptsPlots"]:
    data = real_pdsdbscan.plotData["minptsPlots"][toPlot]
    axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
  axs[i,j].set_ylabel("time (sec)", fontsize=9,labelpad=0)
  axs[i,j].set_xlabel("minpts", fontsize=9, labelpad=0)
  axs[i,j].set_title(GetLetter(accessLetter) + FixName(toPlot)+TitleSuffixEps(list(real1_ndexactbcp.plotData["minptsPlots"][toPlot]["eps"])[0]),fontsize=9.0, fontstyle="italic")
  axs[i,j].grid(ls=':',linewidth=0.3,c='k')
  for l in axs[i,j].yaxis.get_ticklabels():
    l.set_fontsize(9)
  for l in axs[i,j].xaxis.get_ticklabels():
    l.set_fontsize(9)

handles, labels = axs[0,0].get_legend_handles_labels()
for i in range(3):
  for j in range(4):
    handles2, labels2 = axs[i,j].get_legend_handles_labels()
    handles += handles2
    labels += labels2
handles, labels = reorderLegend(handles, labels)
plt.subplots_adjust(left=0.08,right=0.98,bottom=0.07,top=0.93,wspace=0.45,hspace=0.7)
for idx,toPlot in enumerate(realDataFiles):
  i = positions[idx][0]
  j = positions[idx][1]
  box = axs[i,j].get_position()
  #axs[i,j].set_position([box.x0+0.01,box.y0+0.05,box.x1-box.x0,box.y1-box.y0])
  axs[i,j].set_position([box.x0+0.01,box.y0,box.x1-box.x0,box.y1-box.y0])
fig.legend(handles, labels, loc='upper right',ncol=1, bbox_to_anchor=(1, 1),fontsize=9,labelspacing=0.2)
plt.savefig(plotDir + "paper_nd_minpts.pdf", format="pdf")

# EPS
fig, axs = plt.subplots(3, 4, constrained_layout=False, figsize =(8.5,11/2.2))
for j in range(3):
  for i in range(3):
    accessIdx = j*3 + i
    accessLetter = i*3+j
    if accessIdx >= len(datasetOrder):
      axs[i,j].axis("off")
      continue
    toPlot = datasetOrder[accessIdx]
    toPlot = toPlot + "10M"
    printDataName = GetName(toPlot)
    if toPlot in synt1_ndexactquadtree.plotData["epsPlots"]:
      data = synt1_ndexactquadtree.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], c='b', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="0our-exact-qt",mfc="None")
    if toPlot in synt_quadtreemc_quadraticbcp_bucketing.plotData["epsPlots"]:
      data = synt_quadtreemc_quadraticbcp_bucketing.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='g', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="1our-exact-qt-bucketing",mfc="None")
    if toPlot in synt1_ndexactbcp.plotData["epsPlots"]:
      data = synt1_ndexactbcp.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], c='b', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="2our-exact",mfc="None")
    if toPlot in synt3_ndexact_bucketing.plotData["epsPlots"]:
      data = synt3_ndexact_bucketing.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='g', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="3our-exact-bucketing",mfc="None")
    if toPlot in synt1_ndapproxquadtree.plotData["epsPlots"]:
      data = synt1_ndapproxquadtree.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], c='r', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="4our-approx-qt",mfc="None")
    if toPlot in synt_approx_qt_bucketing.plotData["epsPlots"]:
      data = synt_approx_qt_bucketing.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='k', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="5our-approx-qt-bucketing",mfc="None")
    if toPlot in synt1_ndapprox.plotData["epsPlots"]:
      data = synt1_ndapprox.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], c='r', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="6our-approx",mfc="None")
    if toPlot in synt_approx_bucketing.plotData["epsPlots"]:
      data = synt_approx_bucketing.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='k', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="7our-approx-bucketing",mfc="None")
    if toPlot in synt1_hpdbscan.plotData["epsPlots"]:
      data = synt1_hpdbscan.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
    if toPlot in synt1_pdsdbscan.plotData["epsPlots"]:
      data = synt1_pdsdbscan.plotData["epsPlots"][toPlot]
      axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
    axs[i,j].set_ylabel("time (sec)", fontsize=9,labelpad=0)
    axs[i,j].set_xlabel("epsilon", fontsize=9,labelpad=0)
    axs[i,j].set_title(GetLetter(accessLetter) + FixName(printDataName)+TitleSuffixMinpts(list(synt1_ndexactbcp.plotData["epsPlots"][toPlot]["minpts"])[0]),fontsize=9, fontstyle="italic")
    axs[i,j].grid(ls=':',linewidth=0.3,c='k')
    for l in axs[i,j].yaxis.get_ticklabels():
      l.set_fontsize(9)
    for l in axs[i,j].xaxis.get_ticklabels():
      l.set_fontsize(9)

axs[0,3].axis("off")
positions = [(1,3),(2,3)]
for idx,toPlot in enumerate(realDataFiles):
  i = positions[idx][0]
  j = positions[idx][1]
  accessIdx = j*3 + i
  accessLetter = j*3+i-1
  if toPlot in real1_ndexactquadtree.plotData["epsPlots"]:
    data = real1_ndexactquadtree.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], c='b',linewidth=0.3,marker='^',label="0our-exact-qt",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real_quadtreemc_quadraticbcp_bucketing.plotData["epsPlots"]:
    data = real_quadtreemc_quadraticbcp_bucketing.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='g',linewidth=0.3,marker='^',label="1our-exact-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndexactbcp.plotData["epsPlots"]:
    data = real1_ndexactbcp.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], c='b',linewidth=0.3,marker='s',label="2our-exact",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndexactbcpbucketing.plotData["epsPlots"]:
    data = real1_ndexactbcpbucketing.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='g',linewidth=0.3,marker='s',label="3our-exact-bucketing",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndapproxquadtree.plotData["epsPlots"]:
    data = real1_ndapproxquadtree.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], c='r',linewidth=0.3,marker='x',label="4our-approx-qt",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real_approx_qt_bucketing.plotData["epsPlots"]:
    data = real_approx_qt_bucketing.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='k',linewidth=0.3,marker='x',label="5our-approx-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None")
  if toPlot in real1_ndapprox.plotData["epsPlots"]:
    data = real1_ndapprox.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], c='r',linewidth=0.3,marker='+',label="6our-approx",markeredgewidth=0.3,markersize=5,mfc="None")
  if toPlot in real_approx_bucketing.plotData["epsPlots"]:
    data = real_approx_bucketing.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], linestyle='--', c='k',linewidth=0.3,marker='+',label="7our-approx-bucketing",markeredgewidth=0.3,markersize=5,mfc="None")
  if toPlot in real_pdsdbscan.plotData["epsPlots"]:
    data = real_pdsdbscan.plotData["epsPlots"][toPlot]
    axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
  axs[i,j].set_ylabel("time (sec)", fontsize=9,labelpad=0)
  axs[i,j].set_xlabel("epsilon", fontsize=9,labelpad=0)
  axs[i,j].set_title(GetLetter(accessLetter) + FixName(toPlot)+TitleSuffixMinpts(list(real1_ndexactbcp.plotData["epsPlots"][toPlot]["minpts"])[0]),fontsize=9.0, fontstyle="italic")
  axs[i,j].grid(ls=':',linewidth=0.3,c='k')
  for l in axs[i,j].yaxis.get_ticklabels():
    l.set_fontsize(9)
  for l in axs[i,j].xaxis.get_ticklabels():
    l.set_fontsize(9)

handles, labels = axs[0,0].get_legend_handles_labels()
for i in range(3):
  for j in range(4):
    handles2, labels2 = axs[i,j].get_legend_handles_labels()
    handles += handles2
    labels += labels2
handles, labels = reorderLegend(handles, labels)
plt.subplots_adjust(left=0.08,right=0.98,bottom=0.07,top=0.93,wspace=0.45,hspace=0.7)
for idx,toPlot in enumerate(realDataFiles):
  i = positions[idx][0]
  j = positions[idx][1]
  box = axs[i,j].get_position()
  #axs[i,j].set_position([box.x0+0.01,box.y0+0.05,box.x1-box.x0,box.y1-box.y0])
  axs[i,j].set_position([box.x0+0.01,box.y0,box.x1-box.x0,box.y1-box.y0])
fig.legend(handles, labels, loc='upper right', ncol=1, bbox_to_anchor=(1, 1),fontsize=9,labelspacing=0.2)
plt.savefig(plotDir + "paper_nd_eps.pdf", format="pdf")

# SCALABILITY

fig, axs = plt.subplots(3, 4, constrained_layout=False, figsize =(8.5,11/2.2))
for j in range(3):
  for i in range(3):
    myax = axs[i,j]
    accessIdx = j*3 + i
    accessLetter = i*3+j
    if accessIdx >= len(datasetOrder):
      continue
    toPlot = datasetOrder[accessIdx]
    toPlot = toPlot + "10M"
    printDataName = GetName(toPlot)
    if toPlot in synt1_ndexactquadtree.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, c='b',linewidth=0.3,marker='^',label="0our-exact-qt",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
    if toPlot in synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, linestyle='--', c='g',linewidth=0.3,marker='^',label="1our-exact-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
    if toPlot in synt1_ndexactbcp.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, c='b',linewidth=0.3,marker='s',label="2our-exact",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
    if toPlot in synt3_ndexact_bucketing.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt3_ndexact_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt3_ndexact_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, linestyle='--', c='g',linewidth=0.3,marker='s',label="3our-exact-bucketing",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)      
    if toPlot in synt1_ndapproxquadtree.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt1_ndapproxquadtree.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt1_ndapproxquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, c='r',linewidth=0.3,marker='x',label="4our-approx-qt",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)
    if toPlot in synt_approx_qt_bucketing.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt_approx_qt_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt_approx_qt_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, linestyle='--', c='k',linewidth=0.3,marker='x',label="5our-approx-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)
    if toPlot in synt1_ndapprox.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt1_ndapprox.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt1_ndapprox.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, c='r',linewidth=0.3,marker='+',label="6our-approx",markeredgewidth=0.3,markersize=5,mfc="None", basex=2)
    if toPlot in synt_approx_bucketing.plotData["numprocsPlots"]:
      numCores = [math.ceil(i/2) for i in synt_approx_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt_approx_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCores, scalability, linestyle='--', c='k',linewidth=0.3,marker='+',label="7our-approx-bucketing",markeredgewidth=0.3,markersize=5,mfc="None", basex=2)      
    if toPlot in synt1_hpdbscan.plotData["numprocsPlots"]:
      numCoresTest = [math.ceil(i/2) for i in synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='^',label="8hpdbscan",markeredgewidth=0.3,markersize=4,mfc='k',basex=2)
    if toPlot in synt1_pdsdbscan.plotData["numprocsPlots"]:
      numCoresTest = [math.ceil(i/2) for i in synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["x"]]
      runTimes = synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
      scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
      myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='s',label="9pdsdbscan",markeredgewidth=0.3,markersize=3,mfc='k',basex=2)
      
    numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    myax.set_ylabel("speedup over serial-" + "\n" + baselineStr, fontsize=9,labelpad=0)
    myax.set_xlabel("num-threads", fontsize=9,labelpad=0)
    myax.set_title(GetLetter(accessLetter) + FixName(printDataName) + TitleSuffix(list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["eps"])[0], list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["minpts"])[0]),fontsize=9.0, fontstyle="italic")
    myax.grid(ls=':',linewidth=0.3,c='k')
    myax.set_xticks(numCores)
    myax.set_xticklabels(["1","4","8","16","24","36","36h"])
    for l in myax.yaxis.get_ticklabels():
      l.set_fontsize(9)
    for l in myax.xaxis.get_ticklabels():
      l.set_visible(True)
      l.set_fontsize(9)

axs[0,3].axis("off")
positions = [(1,3),(2,3)]
for idx,toPlot in enumerate(realDataFiles):
  if plotIndividual:
    fig, axs = plt.subplots(1, 1, constrained_layout=True, figsize =(8.5/4,11/6))
    myax = axs
  else:
    i = positions[idx][0]
    j = positions[idx][1]
    myax = axs[i,j]
    accessIdx = j*3 + i
    accessLetter = j*3+i-1
  if toPlot in real1_ndexactquadtree.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, c='b',linewidth=0.3,marker='^',label="0our-exact-qt",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
  if toPlot in real_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real_quadtreemc_quadraticbcp_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, linestyle='--', c='g',linewidth=0.3,marker='^',label="1our-exact-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
  if toPlot in real1_ndexactbcp.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, c='b',linewidth=0.3,marker='s',label="2our-exact",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
  if toPlot in real1_ndexactbcpbucketing.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real1_ndexactbcpbucketing.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real1_ndexactbcpbucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, linestyle='--', c='g',linewidth=0.3,marker='s',label="3our-exact-bucketing",markeredgewidth=0.3,markersize=4,mfc="None",basex=2)
  if toPlot in real1_ndapproxquadtree.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real1_ndapproxquadtree.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real1_ndapproxquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, c='r',linewidth=0.3,marker='x',label="4our-approx-qt",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)
  if toPlot in real_approx_qt_bucketing.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real_approx_qt_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real_approx_qt_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, linestyle='--', c='k',linewidth=0.3,marker='x',label="5our-approx-qt-bucketing",markeredgewidth=0.3,markersize=4,mfc="None", basex=2)
  if toPlot in real1_ndapprox.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real1_ndapprox.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real1_ndapprox.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, c='r',linewidth=0.3,marker='+',label="6our-approx",markeredgewidth=0.3,markersize=5,mfc="None", basex=2)
  if toPlot in real_approx_bucketing.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in real_approx_bucketing.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real_approx_bucketing.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCores, scalability, linestyle='--', c='k',linewidth=0.3,marker='+',label="7our-approx-bucketing",markeredgewidth=0.3,markersize=5,mfc="None", basex=2)
  if toPlot  in real_pdsdbscan.plotData["numprocsPlots"]:
    numCoresTest = [math.ceil(i/2) for i in real_pdsdbscan.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = real_pdsdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
    scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
    axs[i,j].semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='s',label="9pdsdbscan",markeredgewidth=0.3,markersize=3,mfc='k',basex=2)
  axs[i,j].set_ylabel("speedup over serial-" + "\n" + baselineStr, fontsize=9,labelpad=0)
  axs[i,j].set_xlabel("num-threads", fontsize=9,labelpad=0)
  axs[i,j].set_title(GetLetter(accessLetter) + FixName(toPlot) +TitleSuffix(list(real1_ndexactbcpbucketing.plotData["numprocsPlots"][toPlot]["eps"])[0], list(real1_ndexactbcpbucketing.plotData["numprocsPlots"][toPlot]["minpts"])[0]),fontsize=9, fontstyle="italic")
  axs[i,j].grid(ls=':',linewidth=0.3,c='k')
  
  numCores = [math.ceil(i/2) for i in real1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["x"]]
  axs[i,j].set_xticks(numCores)
  axs[i,j].set_xticklabels(["1","4","8","16","24","36","36h"])
  for l in axs[i,j].yaxis.get_ticklabels():
    l.set_fontsize(9)
  for l in axs[i,j].xaxis.get_ticklabels():
    l.set_visible(True)
    l.set_fontsize(9)

if not plotIndividual:
  handles, labels = axs[0,0].get_legend_handles_labels()
  for i in range(3):
    for j in range(4):
      handles2, labels2 = axs[i,j].get_legend_handles_labels()
      handles += handles2
      labels += labels2
  handles, labels = reorderLegend(handles, labels)
  plt.subplots_adjust(left=0.08,right=0.98,bottom=0.07,top=0.93,wspace=0.45,hspace=0.7)
  for idx,toPlot in enumerate(realDataFiles):
    i = positions[idx][0]
    j = positions[idx][1]
    box = axs[i,j].get_position()
    #axs[i,j].set_position([box.x0+0.01,box.y0+0.05,box.x1-box.x0,box.y1-box.y0])
    axs[i,j].set_position([box.x0+0.01,box.y0,box.x1-box.x0,box.y1-box.y0])
  fig.legend(handles, labels, loc='upper right', ncol=1, bbox_to_anchor=(1, 1),fontsize=9,labelspacing=0.2)
  plt.savefig(plotDir + "paper_nd_speedup.pdf", format="pdf")

# RHO small

rhoDataOrder = ["5D_VisualSim_", "5D_VisualVar_"]
fig, axs = plt.subplots(1, 2, constrained_layout=False, figsize =(8.5/2.2,11/5.5))
for j in range(2):
  accessIdx = j
  if accessIdx >= len(rhoDataOrder):
    axs[j].axis("off")
    continue
  toPlot = rhoDataOrder[accessIdx]
  toPlot = toPlot + "10M"
  printDataName = GetName(toPlot)
  # find the fastest parallel timing
  fastestExact = 1000000000
  fastestMethod = ""
  if toPlot in synt1_ndexactquadtree.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in synt1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = synt1_ndexactquadtree.plotData["numprocsPlots"][toPlot]["totaltime"]
    if runTimes[-1] < fastestExact:
      fastestExact = runTimes[-1]
      fastestMethod = "0our-exact-qt"
  if toPlot in synt1_ndexactbcp.plotData["numprocsPlots"]:
    numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
    runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
    if runTimes[-1] < fastestExact:
      fastestExact = runTimes[-1]
      fastestMethod = "2our-exact"
  axs[j].axhline(y=fastestExact,linestyle='-',lw=0.5,c='b',label="9our-best-exact")
  if toPlot in synt1_ndapprox_long_rho.plotData["rhoPlots"] and plotApprox:
    data = synt1_ndapprox_long_rho.plotData["rhoPlots"][toPlot]
    axs[j].semilogx(data["x"], data["totaltime"], c='r', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="6our-approx",mfc="None")
  if toPlot in synt1_ndapproxquadtree.plotData["rhoPlots"] and plotApprox:
    data = synt1_ndapproxquadtree.plotData["rhoPlots"][toPlot]
    axs[j].semilogx(data["x"], data["totaltime"], c='r', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="4our-approx-qt",mfc="None")
  axs[j].set_ylabel("time (sec)", fontsize=9)
  axs[j].set_xlabel("rho", fontsize=9)
  axs[j].set_title(GetLetter(accessIdx) + FixName(printDataName)+TitleSuffixMinpts(list(synt1_ndapproxquadtree.plotData["rhoPlots"][toPlot]["minpts"])[0]),fontsize=9, fontstyle="italic")
  axs[j].grid(ls=':',linewidth=0.3,c='k')
  axs[j].set_ylim(bottom=0,top=2)
  for l in axs[j].yaxis.get_ticklabels():
    l.set_fontsize(9)
  for l in axs[j].xaxis.get_ticklabels():
    l.set_fontsize(9)
handles, labels = axs[0].get_legend_handles_labels()
for j in range(2):
    handles2, labels2 = axs[j].get_legend_handles_labels()
    handles += handles2
    labels += labels2
handles, labels = reorderLegend(handles, labels)
plt.subplots_adjust(left=0.12,right=0.95,bottom=0.2,top=0.68,wspace=0.45,hspace=0.7)
fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 0.99),fontsize=9,labelspacing=0.2)
plt.savefig(plotDir + "paper_nd_rho_small.pdf", format="pdf")

# NUMPTS small
# numptsDataOrder = ["5D_VisualSim_", "5D_VisualVar_"]
# fig, axs = plt.subplots(1, 2, constrained_layout=False, figsize =(8.5/2,11/4.5))
# for j in range(2):
#     accessIdx = j
#     if accessIdx >= len(datasetOrder):
#       axs[i,j].axis('off')
#       continue
#     dimName = numptsDataOrder[accessIdx]
#     toPlot = dimName + "10M"
#     printDataName = GetName(toPlot)
#     if dimName in synt1_ndexactquadtree.plotData["numptsPlots"]:
#       data = synt1_ndexactquadtree.plotData["numptsPlots"][dimName]
#       axs[j].loglog(data["x"], data["totaltime"], c='b', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="0our-exact-qt",mfc="None")
#     if dimName in synt1_ndexactbcp.plotData["numptsPlots"]:
#       data = synt1_ndexactbcp.plotData["numptsPlots"][dimName]
#       axs[j].loglog(data["x"], data["totaltime"], c='b', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="2our-exact",mfc="None")
#     if dimName in synt1_ndapproxquadtree.plotData["numptsPlots"]:
#       data = synt1_ndapproxquadtree.plotData["numptsPlots"][dimName]
#       axs[j].loglog(data["x"], data["totaltime"], c='r', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="4our-approx-qt",mfc="None")
#     if dimName in synt1_ndapprox.plotData["numptsPlots"]:
#       data = synt1_ndapprox.plotData["numptsPlots"][dimName]
#       axs[j].loglog(data["x"], data["totaltime"], c='r', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="6our-approx",mfc="None")
#     if dimName in synt1_hpdbscan.plotData["numptsPlots"]:
#       data = synt1_hpdbscan.plotData["numptsPlots"][dimName]
#       axs[j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
#     if dimName in synt1_pdsdbscan.plotData["numptsPlots"]:
#       data = synt1_pdsdbscan.plotData["numptsPlots"][dimName]
#       axs[j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
#     axs[j].set_ylabel("time (sec)", fontsize=9)
#     axs[j].set_xlabel("num-pts", fontsize=9)
#     axs[j].set_title(GetLetter(accessIdx) + FixName(printDataName)+"\n(various params)",fontsize=9, fontstyle="italic")
#     axs[j].grid(ls=':',linewidth=0.3,c='k')
#     for l in axs[j].yaxis.get_ticklabels():
#       l.set_fontsize(9)
#     for l in axs[j].xaxis.get_ticklabels():
#       l.set_fontsize(9)
# handles, labels = axs[0].get_legend_handles_labels()
# for j in range(2):
#     handles2, labels2 = axs[j].get_legend_handles_labels()
#     handles += handles2
#     labels += labels2
# handles, labels = reorderLegend(handles, labels)
# plt.subplots_adjust(left=0.13,right=0.95,bottom=0.18,top=0.6,wspace=0.45,hspace=0.7)
# fig.legend(handles, labels, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1),fontsize=9,labelspacing=0.2)
# plt.savefig(plotDir + "paper_nd_numpts_small.pdf", format="pdf")

# NUMPTS, will not include real datasets
# fig, axs = plt.subplots(3, 4, constrained_layout=False, figsize =(8.5,11/2))
# for j in range(4):
#   for i in range(3):
#     accessIdx = j*3 + i
#     accessLetter = i*4+j
#     if accessIdx >= len(datasetOrder):
#       axs[i,j].axis('off')
#       continue
#     dimName = datasetOrder[accessIdx]
#     toPlot = dimName + "10M"
#     printDataName = GetName(toPlot)
#     # if False:
#     #   Plot2D(axs, i, j, dimName, "numptsPlots", axs[i,j].loglog)
#     if dimName in synt1_ndexactquadtree.plotData["numptsPlots"]:
#       data = synt1_ndexactquadtree.plotData["numptsPlots"][dimName]
#       axs[i,j].loglog(data["x"], data["totaltime"], c='b', marker='D',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="0our-exact-qt",mfc="None")
#     if dimName in synt1_ndexactbcp.plotData["numptsPlots"]:
#       data = synt1_ndexactbcp.plotData["numptsPlots"][dimName]
#       axs[i,j].loglog(data["x"], data["totaltime"], c='b', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="2our-exact",mfc="None")
#     if dimName in synt3_ndexact_bucketing.plotData["numptsPlots"]:
#       data = synt3_ndexact_bucketing.plotData["numptsPlots"][dimName]
#       axs[i,j].semilogy(data["x"], data["totaltime"], c='g', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="3our-exact-bucketing",mfc="None")      
#     if dimName in synt1_ndapproxquadtree.plotData["numptsPlots"]:
#       data = synt1_ndapproxquadtree.plotData["numptsPlots"][dimName]
#       axs[i,j].loglog(data["x"], data["totaltime"], c='r', marker='x',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="4our-approx-qt",mfc="None")
#     if dimName in synt1_ndapprox.plotData["numptsPlots"]:
#       data = synt1_ndapprox.plotData["numptsPlots"][dimName]
#       axs[i,j].loglog(data["x"], data["totaltime"], c='r', marker='+',linewidth=0.3,markeredgewidth=0.3,markersize=5,label="6our-approx",mfc="None")
#     if dimName in synt1_hpdbscan.plotData["numptsPlots"]:
#       data = synt1_hpdbscan.plotData["numptsPlots"][dimName]
#       axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
#     if dimName in synt1_pdsdbscan.plotData["numptsPlots"]:
#       data = synt1_pdsdbscan.plotData["numptsPlots"][dimName]
#       axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
#     axs[i,j].set_ylabel("time (sec)", fontsize=9)
#     axs[i,j].set_xlabel("num-pts", fontsize=9)
#     axs[i,j].set_title(GetLetter(accessLetter) + FixName(printDataName),fontsize=9, fontstyle="italic")
#     axs[i,j].grid(ls=':',linewidth=0.3,c='k')
#     for l in axs[i,j].yaxis.get_ticklabels():
#       l.set_fontsize(9)
#     for l in axs[i,j].xaxis.get_ticklabels():
#       l.set_fontsize(9)

# handles, labels = axs[0,0].get_legend_handles_labels()
# handles.reverse()
# labels.reverse()
# handles, labels = axs[0,0].get_legend_handles_labels()
# for i in range(3):
#   for j in range(4):
#     handles2, labels2 = axs[i,j].get_legend_handles_labels()
#     handles += handles2
#     labels += labels2
# handles, labels = reorderLegend(handles, labels)
# plt.subplots_adjust(left=0.08,right=0.98,bottom=0.07,top=0.85,wspace=0.45,hspace=0.7)
# fig.legend(handles, labels, loc='upper center',ncol=4, bbox_to_anchor=(0.5, 1),fontsize=9,labelspacing=0.2)
# plt.savefig(plotDir + "paper_nd_numpts.pdf", format="pdf")




# for 2d plots, do ss-simden and ssvarden, all four plots
fig, axs = plt.subplots(2, 4, constrained_layout=False, figsize =(8.5,11/2.8))
toPlot="2D_VisualSim_10M"
dimName = "2D_VisualSim_"
i=0
j=0
accessIdx = i * 4 + j
myax=axs[i,j] #eps
Plot2D(axs, i, j, toPlot, "epsPlots", axs[i,j].semilogy)
if toPlot in synt1_hpdbscan.plotData["epsPlots"] and plotHp:
  data = synt1_hpdbscan.plotData["epsPlots"][toPlot]
  axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
if toPlot in synt1_pdsdbscan.plotData["epsPlots"] and plotPds:
  data = synt1_pdsdbscan.plotData["epsPlots"][toPlot]
  axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
axs[i,j].set_ylabel("time (sec)", fontsize=9)
axs[i,j].set_xlabel("epsilon", fontsize=9)
axs[i,j].set_title(GetLetter(accessIdx) + FixName(toPlot)+TitleSuffixMinpts(list(synt1_2dgridbcp.plotData["epsPlots"][toPlot]["minpts"])[0]),fontsize=9, fontstyle="italic")
axs[i,j].grid(ls=':',linewidth=0.3,c='k')
for l in axs[i,j].yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in axs[i,j].xaxis.get_ticklabels():
  l.set_fontsize(9)

i=0
j=1
accessIdx = i * 4 + j
myax=axs[i,j] #minpts
Plot2D(axs, i, j, toPlot, "minptsPlots", axs[i,j].loglog)
if toPlot in synt1_hpdbscan.plotData["minptsPlots"] and plotHp:
  data = synt1_hpdbscan.plotData["minptsPlots"][toPlot]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
if toPlot in synt1_pdsdbscan.plotData["minptsPlots"] and plotPds:
  data = synt1_pdsdbscan.plotData["minptsPlots"][toPlot]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
axs[i,j].set_ylabel("time (sec)", fontsize=9)
axs[i,j].set_xlabel("minpts", fontsize=9)
axs[i,j].set_title(GetLetter(accessIdx) + FixName(toPlot)+TitleSuffixEps(list(synt1_ndexactbcp.plotData["minptsPlots"][toPlot]["eps"])[0]),fontsize=9, fontstyle="italic")
axs[i,j].grid(ls=':',linewidth=0.3,c='k')
for l in axs[i,j].yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in axs[i,j].xaxis.get_ticklabels():
  l.set_fontsize(9)

i=0
j=2
accessIdx = i * 4 + j
myax=axs[i,j] #numpts
Plot2D(axs, i, j, dimName, "numptsPlots", axs[i,j].loglog)
if dimName in synt1_hpdbscan.plotData["numptsPlots"] and plotHp:
  data = synt1_hpdbscan.plotData["numptsPlots"][dimName]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
if dimName in synt1_pdsdbscan.plotData["numptsPlots"] and plotPds:
  data = synt1_pdsdbscan.plotData["numptsPlots"][dimName]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
axs[i,j].set_ylabel("time (sec)", fontsize=9)
axs[i,j].set_xlabel("num-pts", fontsize=9)
axs[i,j].set_title(GetLetter(accessIdx) + FixName(dimName)+"NumPts\n"+"(various params)",fontsize=9.0, fontstyle="italic")
axs[i,j].grid(ls=':',linewidth=0.3,c='k')
for l in axs[i,j].yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in axs[i,j].xaxis.get_ticklabels():
  l.set_fontsize(9)

i=0
j=3
accessIdx = i * 4 + j
myax=axs[i,j] #speedup
Plot2DSpeedup(myax, toPlot)
if toPlot in synt1_hpdbscan.plotData["numprocsPlots"] and plotHp:
  numCoresTest = [math.ceil(i/2) for i in synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
  myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='^',label="8hpdbscan",markeredgewidth=0.3,markersize=4,mfc='k',basex=2)
if toPlot  in synt1_pdsdbscan.plotData["numprocsPlots"] and plotPds:
  numCoresTest = [math.ceil(i/2) for i in synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
  myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='s',label="9pdsdbscan",markeredgewidth=0.3,markersize=3,mfc='k',basex=2)
numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
myax.set_ylabel("speedup over serial-" + "\n" + baselineStr, fontsize=9)
myax.set_xlabel("num-threads", fontsize=9)
myax.set_title(GetLetter(accessIdx) + FixName(toPlot) + TitleSuffix(list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["eps"])[0], list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["minpts"])[0]),fontsize=9.0, fontstyle="italic")
myax.grid(ls=':',linewidth=0.3,c='k')
myax.set_xticks(numCores)
myax.set_xticklabels(["1","4","8","16","24","36","36h"])
for l in myax.yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in myax.xaxis.get_ticklabels():
  l.set_visible(True)
  l.set_fontsize(9)

toPlot="2D_VisualVar_10M"
dimName = "2D_VisualVar_"
i=1
j=0
accessIdx = i * 4 + j
myax=axs[i,j] #eps
Plot2D(axs, i, j, toPlot, "epsPlots", axs[i,j].semilogy)
if toPlot in synt1_hpdbscan.plotData["epsPlots"] and plotHp:
  data = synt1_hpdbscan.plotData["epsPlots"][toPlot]
  axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
if toPlot in synt1_pdsdbscan.plotData["epsPlots"] and plotPds:
  data = synt1_pdsdbscan.plotData["epsPlots"][toPlot]
  axs[i,j].semilogy(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
axs[i,j].set_ylabel("time (sec)", fontsize=9)
axs[i,j].set_xlabel("epsilon", fontsize=9)
axs[i,j].set_title(GetLetter(accessIdx) + FixName(toPlot)+TitleSuffixMinpts(list(synt1_2dgridbcp.plotData["epsPlots"][toPlot]["minpts"])[0]),fontsize=9, fontstyle="italic")
axs[i,j].grid(ls=':',linewidth=0.3,c='k')
for l in axs[i,j].yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in axs[i,j].xaxis.get_ticklabels():
  l.set_fontsize(9)

i=1
j=1
accessIdx = i * 4 + j
myax=axs[i,j] #minpts
Plot2D(axs, i, j, toPlot, "minptsPlots", axs[i,j].loglog)
if toPlot in synt1_hpdbscan.plotData["minptsPlots"] and plotHp:
  data = synt1_hpdbscan.plotData["minptsPlots"][toPlot]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
if toPlot in synt1_pdsdbscan.plotData["minptsPlots"] and plotPds:
  data = synt1_pdsdbscan.plotData["minptsPlots"][toPlot]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
axs[i,j].set_ylabel("time (sec)", fontsize=9)
axs[i,j].set_xlabel("minpts", fontsize=9)
axs[i,j].set_title(GetLetter(accessIdx) + FixName(toPlot)+TitleSuffixEps(list(synt1_ndexactbcp.plotData["minptsPlots"][toPlot]["eps"])[0]),fontsize=9, fontstyle="italic")
axs[i,j].grid(ls=':',linewidth=0.3,c='k')
for l in axs[i,j].yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in axs[i,j].xaxis.get_ticklabels():
  l.set_fontsize(9)

i=1
j=2
accessIdx = i * 4 + j
myax=axs[i,j] #numpts
Plot2D(axs, i, j, dimName, "numptsPlots", axs[i,j].loglog)
if dimName in synt1_hpdbscan.plotData["numptsPlots"] and plotHp:
  data = synt1_hpdbscan.plotData["numptsPlots"][dimName]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='^',linewidth=0.3,markeredgewidth=0.3,markersize=4,label="8hpdbscan",mfc='k')
if dimName in synt1_pdsdbscan.plotData["numptsPlots"] and plotPds:
  data = synt1_pdsdbscan.plotData["numptsPlots"][dimName]
  axs[i,j].loglog(data["x"], data["totaltime"], c='k', marker='s',linewidth=0.3,markeredgewidth=0.3,markersize=3,label="9pdsdbscan",mfc='k')
axs[i,j].set_ylabel("time (sec)", fontsize=9)
axs[i,j].set_xlabel("num-pts", fontsize=9)
axs[i,j].set_title(GetLetter(accessIdx) + FixName(dimName)+"NumPts\n"+"(various params)",fontsize=9.0, fontstyle="italic")
axs[i,j].grid(ls=':',linewidth=0.3,c='k')
for l in axs[i,j].yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in axs[i,j].xaxis.get_ticklabels():
  l.set_fontsize(9)

i=1
j=3
accessIdx = i * 4 + j
myax=axs[i,j] #speedup
Plot2DSpeedup(myax, toPlot)
if toPlot in synt1_hpdbscan.plotData["numprocsPlots"] and plotHp:
  numCoresTest = [math.ceil(i/2) for i in synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_hpdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
  myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='^',label="8hpdbscan",markeredgewidth=0.3,markersize=4,mfc='k',basex=2)
if toPlot  in synt1_pdsdbscan.plotData["numprocsPlots"] and plotPds:
  numCoresTest = [math.ceil(i/2) for i in synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["x"]]
  runTimes = synt1_pdsdbscan.plotData["numprocsPlots"][toPlot]["totaltime"]
  scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
  myax.semilogx(numCoresTest, scalability, c='k',linewidth=0.3,marker='s',label="9pdsdbscan",markeredgewidth=0.3,markersize=3,mfc='k',basex=2)
numCores = [math.ceil(i/2) for i in synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["x"]]
runTimes = synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["totaltime"]
scalability,baselineStr,methodStr = GetSpeedups(toPlot, runTimes)
myax.set_ylabel("speedup over serial-" + "\n" + baselineStr, fontsize=9)
myax.set_xlabel("num-threads", fontsize=9)
myax.set_title(GetLetter(accessIdx) + FixName(toPlot) + TitleSuffix(list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["eps"])[0], list(synt1_ndexactbcp.plotData["numprocsPlots"][toPlot]["minpts"])[0]),fontsize=9.0, fontstyle="italic")
myax.grid(ls=':',linewidth=0.3,c='k')
myax.set_xticks(numCores)
myax.set_xticklabels(["1","4","8","16","24","36","36h"])
for l in myax.yaxis.get_ticklabels():
  l.set_fontsize(9)
for l in myax.xaxis.get_ticklabels():
  l.set_visible(True)
  l.set_fontsize(9)

handles, labels = axs[0,0].get_legend_handles_labels()
for i in range(2):
  for j in range(4):
    handles2, labels2 = axs[i,j].get_legend_handles_labels()
    handles += handles2
    labels += labels2
handles, labels = reorderLegend(handles, labels)
plt.subplots_adjust(left=0.08,right=0.98,bottom=0.1,top=0.8,wspace=0.45,hspace=0.7)
fig.legend(handles, labels, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1),fontsize=9,labelspacing=0.2)
plt.savefig(plotDir + "paper_2dplots.pdf", format="pdf")

