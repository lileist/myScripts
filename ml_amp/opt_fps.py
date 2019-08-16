"""
example for how to use double ended connect to connect the minima in an existing database

we will use as an example system the Lennard-Jones cluster with a small number of atoms.
Since we don't already have a database, for this example we'll build a small one using
basinhopping"
"""
import numpy as np
import sys
from amp import Amp
from amp.descriptor.cutoffs import Cosine,dict2cutoff
import ase.io

from amp.descriptor.gaussian import Gaussian

args = sys.argv
#read in gr function
Rc =  4.0
gr = []
Rmin = 0.0
Rij = 0.0
dr = 0.02
radium = []
gr_file=open(args[1], 'r')
lines = gr_file.readlines()

n_dotall = int((Rc)/dr)

for i in range(n_dotall):
   fields = lines[i].split()
   if float(fields[0]) >= Rmin and float(fields[1])<= Rc:
      gr.append(float(fields[1]))
      radium.append(Rij)
   Rij += dr
gr=np.array(gr)
radium = np.array(radium)

unity = np.full((len(gr),),1.0)
#print unity
#Define cutoff function
cutoff={'name': 'Cosine', 'kwargs': {'Rc': Rc}}
cutoff_fxn = dict2cutoff(cutoff)

step = 0

def get_loss(parametervector):
    global step
    numb_fps = len(parametervector) / 2
    #fps = calculate_fps(parametervector.reshape((numb_fps, 2)))
    fps = calculate_fps(parametervector)
    completeness=None
    g2s_gr = {}
    fp_gr_sum = 0
    for key in fps.keys():
       g2s_gr[key] = fps[key] * gr
       #normalize independantly
       g2s_gr[key] /= np.sum(g2s_gr[key])
       if completeness is None:
          completeness = g2s_gr[key]  ** 2
       else:
          completeness += g2s_gr[key] **2
    keys = fps.keys()
    inner_product = 0
    for i in range(len(keys)-1):
       for j in range(i+1, len(keys)):
         inner_product += np.dot(g2s_gr[keys[i]], g2s_gr[keys[j]])
    if step == 0:
       print inner_product
    #diff = np.sum(unity-completeness)/float(len(gr))
    diff = inner_product
    print "{:6d} {:14.6f} {:14.6f} {:14.6f}".format(step, diff, min(parametervector), max(parametervector))
    #print parametervector
    step += 1
    #return np.sum(unity-completeness)
    return diff

#def calculate_fps(eta_rss):
def calculate_fps(etas):
    global Rmin
    g2s = {}
    g2_orig={}
    keys=[]
    #print cutoff_fxn
    #fp_sum = 0
    #fix Rs
    Rs   = [2.5, 2.8, 3.0,  3.2,  3.4, 3.6, 3.8]
    n_start =None
    n_end =None
    for eta_rs in zip(etas, Rs):
       Rij=-0.01
       key = str(eta_rs[0])+'_'+str(eta_rs[1])
       keys.append(key)
       g2s[key] = []
       for i in range(n_dotall):
          Rij += 0.02
          args_cutoff_fxn = dict(Rij=Rij)
          #if cutoff['name'] == 'Polynomial':
          #    args_cutoff_fxn['gamma'] = cutoff['kwargs']['gamma']
          term =  np.exp(-eta_rs[0] * ((Rij-eta_rs[1]) ** 2.) / (Rc ** 2.)) * \
                    cutoff_fxn(**args_cutoff_fxn)
          g2s[key].append(term)
          if Rij >= Rmin and n_start is None:
             n_start = i
       g2s[key] = np.array(g2s[key])
       #normalize fps independently
       #g2s[key] = g2_orig[key] / sum(g2_orig[key])
       g2s[key] = g2s[key][n_start:-1:1]
       g2s[key] = g2s[key] / np.sum(g2s[key])
       #normalize fps overall
       #fp_sum += sum(g2_orig[key])
    #for key in g2s.keys():
    #   g2s[key] = g2_orig[key] / fp_sum
    return g2s

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
          output += "{:12.8f}".format(g2s[etas[j]][i])
       g_2.write("{:s} {:12.8f}\n".format(output, gr[i]))
    return
   
def print_expectation(g2s, thetas):
    output = ""
    for k in g2s:
       output += "{:12.8f}".format(np.dot(g2s[k], thetas))
    print output


def main():
  #Initialize parameters
  #Rs   = [0., 0.5, 1.0,  1.5,  0.0, 0., 0.5]
  Rs   = [2.5, 2.8, 3.0,  3.2,  3.4, 3.6, 3.8]
  etas = [5, 10.,  20., 40., 100., 200., 200.0]
  parametervector = []
  for i in range(len(Rs)):
     parametervector.append(etas[i])
     #parametervector.append(Rs[i])
 
  #x0 = np.array(parametervector)
  x0 = parametervector
 
 
  optimizer = 'basinhopping'
 
  #from scipy.optimize import minimize as optimizer
  #optimizer_kwargs = {
  #                    'method' : 'BFGS',
  #                    'options': {'gtol': 1e-15, }
  #                   }
  if optimizer == 'BFGS':
     from scipy.optimize import minimize as optimizer
     optimizer_kwargs = {
                         'method' : 'BFGS',
                         'options': {'gtol': 1e-15, }
                        }
     #optimizer_kwargs = {'method':'BFGS', 'gtol': 1e-15, }
  elif optimizer == 'basinhopping':
     from scipy.optimize import basinhopping as optimizer
     optimizer(get_loss, x0, 
               niter=500,
               T = 0.1,
               stepsize = 0.01
               )
  elif optimizer == 'L-BFGS-B':
     from scipy.optimize import minimize as optimizer
     optimizer_kwargs = {
                         'method': 'L-BFGS-B',
                         'options': {'ftol': 1e-05,
                                     'gtol': 1e-08,
                                     'maxfun': 1000000,
                                     'maxiter': 1000000}
                        }
     import scipy
     from distutils.version import StrictVersion
     if StrictVersion(scipy.__version__) >= StrictVersion('0.17.0'):
         optimizer_kwargs['options']['maxls'] = 2000
  #optimizer(get_loss, x0, **optimizer_kwargs)
  #try:
  #except:
  #   print "Optimizer error"
  #output_start = args[1].split('.')[0]+'_'
  #log(radium, g2s, gr, keys, output_start+'g2s.dat')
  #log(radium, g2s_gr, gr, keys, output_start+'g2s_gr.dat')

  #print_expectation(g2s_gr, radium)

if __name__ == '__main__':
    main()

