"""
example for how to use double ended connect to connect the minima in an existing database

we will use as an example system the Lennard-Jones cluster with a small number of atoms.
Since we don't already have a database, for this example we'll build a small one using
basinhopping"
"""
import numpy as np
import sys
from ase import Atoms
from ase.calculators.lj import LennardJones
from amp import Amp
#from ase.calculators.lammpslib import LAMMPSlib
from amp.descriptor.cutoffs import Cosine,dict2cutoff
import ase.io


from amp.descriptor.gaussian import Gaussian

def readinputs(filename):
    f=open(filename, 'r')
    parameters = {}
    lines=f.readlines()
    for line in lines:
      if line.startswith('#'):
         continue
      fields = line.split('=')
      parameters[fields[0].strip()]=fields[1].replace("\n","").strip()
    return parameters

def type_convertion(para, value):
    if type(para) is int:
       return int(value)
    if type(para) is float:
       return float(value)
    if type(para) is str:
       return value
    if type(para) is list:
       return value.split()
    if type(para) is bool:
       return ast.literal_eval(value)


def log(thetas, g2s, gr, etas, filename):
   g_2 = open(filename,'w')
   output = "#     {:12s}".format('r')
   for i in range(len(etas)):
      output += "{:12s}".format(etas[i])
   g_2.write("{:s} {:12s}\n".format(output, "gr"))
   for i in range(len(thetas)):
      output = "{:12.6f}".format(thetas[i])
      for j in range(len(etas)):
         output += "{:14.12f}".format(g2s[etas[j]][i])
      g_2.write("{:s} {:12.8f}\n".format(output, gr[i]))
   return
   
def print_expectation(g2s, thetas):
    output = ""
    for k in g2s:
       output += "{:12.8f}".format(np.dot(g2s[k], thetas))
    print output

default_paras = dict(
    Rc = 3.0,
    etas = [10.0, 20.0,  40.0,  80.0,  160.0],
    Rs =   [0.0,   0.5,   1.0,   2.0,    3.0],
    gamma = 1,
    eta = 0.005,
    )
args = sys.argv
gr_file=open(args[1], 'r')
auto = False
if len(args)==6:
  auto = True
  eta = float(args[2])
  Rmin = float(args[3])
  Rmax = float(args[4])
  dRs = float(args[5])

paras=readinputs('input.ini')

for para in default_paras:
   try:
     paras[para] = type_convertion(default_paras[para], paras[para])
   except:
     paras[para] = default_paras[para]

Rc = paras['Rc']
cutoff={'name': 'Cosine', 'kwargs': {'Rc': Rc}}
cutoff_fxn = dict2cutoff(cutoff)

radium = []
lines = gr_file.readlines()
gr = []
Rij = 0
dr = 0.02
n_dot = int(Rc/dr)
for i in range(n_dot):
   gr.append(float(lines[i].split()[1]))
   radium.append(Rij)
   Rij += dr
gr=np.array(gr)
radium = np.array(radium)

g2s = {}
g2_orig={}
g2s_gr = {}
keys=[]
etas = [float(field) for field in paras['etas']]
Rs = [float(field) for field in paras['Rs']]

if auto:
  etas = []
  nRs = int((Rmax-Rmin)/dRs)
  Rs = []
  for i in range(nRs):
     etas.append(eta)
     Rs.append(Rmin+i*dRs)
  print Rs
  print etas


for eta_rs in zip(etas, Rs):
    Rij=0
    key = str(eta_rs[0])+'_'+str(eta_rs[1])
    keys.append(key)
    g2s[key] = []
    for i in range(n_dot):
       Rij += 0.02
       args_cutoff_fxn = dict(Rij=Rij)
       #if cutoff['name'] == 'Polynomial':
       #    args_cutoff_fxn['gamma'] = cutoff['kwargs']['gamma']
       term =  np.exp(-eta_rs[0] * ((Rij-eta_rs[1]) ** 2.) / (Rc ** 2.)) * \
                 cutoff_fxn(**args_cutoff_fxn)
       g2s[key].append(term)
    g2_orig[key] = np.array(g2s[key])
    g2s[key] = g2_orig[key] / sum(g2_orig[key])

    g2s_gr[key] = g2s[key] * gr
    g2s_gr[key] /= sum(g2s_gr[key])

output_start = args[1].split('.')[0]+'_'
log(radium, g2s, gr, keys, output_start+'g2s.dat')
log(radium, g2s_gr, gr, keys, output_start+'g2s_gr.dat')

print_expectation(g2s_gr, radium)
