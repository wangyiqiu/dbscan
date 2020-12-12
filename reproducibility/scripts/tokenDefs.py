import re

test_split = re.compile("(>testsplit\*\*)")
round_split = re.compile("(>roundsplit\*\*)")

# timings
assigngridtime_time = re.compile(">assigngridtime\*(.+)\*")
markcore_time = re.compile(">markcoretime\*(.+)\*")
buildcoreptquadtree_time = re.compile(">buildcoreptquadtree\*(.+)\*")
buildgetnbrkdt_time = re.compile(">buildgetnbrkdt\*(.+)\*")
clustercore_time = re.compile(">clustercoretime\*(.+)\*")
assignborder_time = re.compile(">assignbordertime\*(.+)\*")
cleanup_time = re.compile(">cleanuptime\*(.+)\*")
total_time = re.compile(">totaltime\*(.+)\*")

# quantities
num_grids = re.compile(">numgrids\*(.+)\*")
num_core = re.compile(">numcore\*(.+)\*")
num_noise = re.compile(">numnoise\*(.+)\*")
num_border = re.compile(">numborder\*(.+)\*")
num_clusterpts = re.compile(">numclusterpts\*(.+)\*")
num_pts = re.compile(">numpts\*(.+)\*")
num_procs = re.compile(">numprocs\*(.+)\*")
num_clusters = re.compile(">numcluster\*(.+)\*")
num_noise = re.compile(">numnoise\*(.+)\*")

# algos
serial_algo = re.compile(">serialalgo\*(.+)\*")
nd_algo = re.compile(">ndalgo\*(.+)\*")
algo = re.compile(">algo\*(.+)\*")

# params
param_datafile = re.compile(">datafile\*(.+)\*")
param_eps = re.compile(">eps\*(.+)\*")
param_minpts = re.compile(">minpts\*(.+)\*")
param_rho = re.compile(">rho\*(.+)\*")
param_dim = re.compile(">dim\*(.+)\*")
param_rounds = re.compile(">rounds\*(.+)\*")
# = re.compile(">\*(.+)\*")

# additional timers
timer1 = re.compile(">timer1\*(.+)\*")
timer2 = re.compile(">timer2\*(.+)\*")
timer3 = re.compile(">timer3\*(.+)\*")
timer4 = re.compile(">timer4\*(.+)\*")
timer5 = re.compile(">timer5\*(.+)\*")
timer6 = re.compile(">timer6\*(.+)\*")
timer7 = re.compile(">timer7\*(.+)\*")
timer8 = re.compile(">timer8\*(.+)\*")
timer9 = re.compile(">timer9\*(.+)\*")
timer10 = re.compile(">timer10\*(.+)\*")

def ExtractString(t_target, t_regex):
  found = t_regex.findall(t_target)
  if len(found) == 0:
    return ""
  else:
    return found[0]

def ExtractNumber(t_target, t_regex):
  found = t_regex.findall(t_target)
  if len(found) == 0:
    return -1
  else:
    if (float(found[0]) == 0):
      return 1e-10
    return float(found[0])
