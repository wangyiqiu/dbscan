plotData = dict()

# serial test only c5.18x
# each data set uses default parameter scaled

plotData["2D_VisualSim_10M"] = {"algo":"algo1", "eps":400, "minpts":100, "totaltime":5.411256075, "numclusters":9}
plotData["3D_VisualSim_10M"] = {"algo":"algo1", "eps":1000, "minpts":10, "totaltime":15.160243034, "numclusters":11}
plotData["5D_VisualSim_10M"] = {"algo":"algo1", "eps":1000, "minpts":100, "totaltime":84.108774185, "numclusters":12}
plotData["7D_VisualSim_10M"] = {"algo":"algo1", "eps":2000, "minpts":10, "totaltime":238.551452160, "numclusters":12}
plotData["2D_VisualVar_10M"] = {"algo":"algo1", "eps":1000, "minpts":100, "totaltime":6.072793961, "numclusters":1}
plotData["3D_VisualVar_10M"] = {"algo":"algo1", "eps":2000, "minpts":100, "totaltime":26.407938957, "numclusters":3}
plotData["5D_VisualVar_10M"] = {"algo":"algo1", "eps":3000, "minpts":10, "totaltime":83.766330004, "numclusters":7}
plotData["7D_VisualVar_10M"] = {"algo":"algo1", "eps":3000, "minpts":10, "totaltime":152.037672043, "numclusters":7}
plotData["2D_UniformDisc_10M"] = {"algo":"algo1", "eps":333.335, "minpts":100, "totaltime":8.764185905, "numclusters":5}
plotData["3D_UniformDisc_10M"] = {"algo":"algo1", "eps":666.67, "minpts":100, "totaltime":72.350967169, "numclusters":5}
plotData["5D_UniformDisc_10M"] = {"algo":"algo1", "eps":6666.7, "minpts":100, "totaltime":153.002969027, "numclusters":5}
plotData["7D_UniformDisc_10M"] = {"algo":"algo1", "eps":13333.4, "minpts":100, "totaltime":550.621845007, "numclusters":5}
plotData["2D_GaussianDisc_10M"] = {"algo":"algo1", "eps":333.335, "minpts":100, "totaltime":10.758620024, "numclusters":5}
plotData["3D_GaussianDisc_10M"] = {"algo":"algo1", "eps":666.67, "minpts":100, "totaltime":61.330317736, "numclusters":5}
plotData["5D_GaussianDisc_10M"] = {"algo":"algo1", "eps":6666.7, "minpts":100, "totaltime":167.281432867, "numclusters":5}
plotData["7D_GaussianDisc_10M"] = {"algo":"algo1", "eps":13333.4, "minpts":100, "totaltime":504.613725901, "numclusters":5}
#plotData["2D_UniformFill_10M"] = {"algo":"algo1", "eps":63245.6, "minpts":10, "totaltime":-1, "numclusters":-1} # floating point exception
#plotData["3D_UniformFill_10M"] = {"algo":"algo1", "eps":63245.6, "minpts":10, "totaltime":-1, "numclusters":-1} # floating point exception
#plotData["5D_UniformFill_10M"] = {"algo":"algo1", "eps":63245.6, "minpts":100, "totaltime":-1, "numclusters":-1} # floating point exception
plotData["7D_UniformFill_10M"] = {"algo":"algo1", "eps":63245.6, "minpts":10, "totaltime":1170.047274828, "numclusters":1}
plotData["2D_UniformStrip_10M"] = {"algo":"algo1", "eps":166.665, "minpts":100, "totaltime":11.445612907, "numclusters":3}
plotData["2D_GaussianStrip_10M"] = {"algo":"algo1", "eps":333.33, "minpts":100, "totaltime":11.087140083, "numclusters":4}
plotData["7D_HouseHold_2M"] = {"algo":"algo1", "eps":2000, "minpts":100, "totaltime":85.829560995, "numclusters":27} # scaled to gan range in each dimension
#plotData[""] = {"algo":"algo1", "eps":, "minpts":, "totaltime":, "numclusters":}
#plotData[""] = {"algo":"algo1", "eps":, "minpts":, "totaltime":, "numclusters":}
