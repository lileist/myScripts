"""
This code is used to calculate G4 function and then correlated with gr function
Input example:
Rc = 3.0
zetas = 10.0  20.0  40.0  80.0  160.0
gamma = 1
eta = 0.005
"""
import numpy as np
import sys, math
from ase import Atoms
from ase.calculators.lj import LennardJones
from amp import Amp
#from ase.calculators.lammpslib import LAMMPSlib
from amp.descriptor.cutoffs import Cosine,dict2cutoff
import ase.io
import sys
from amp.descriptor.gaussian import Gaussian
#from gaussian import Gaussian
 

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
"""
def calculate_G4(Rc, cutoff_fxn, theta_s,zeta, eta, gamma):
    fps = []
    ymin = -3.0
    xmin = -3.0
    R0 = np.array([0.0,0.0])
    for i in range(150):
        for j in range(150):
            nb_1 = np.array([xmin + i * 0.02, ymin + j * 0.02])
            nb_2 = np.array([xmin + j * 0.02, ymin + i * 0.02])
            R01_vector = nb_1 - R0
            R01 = np.linalg.norm(R01_vector)
            R02_vector = nb_2 - R0
            R02 = np.linalg.norm(R02_vector)
            R12_vector = nb_2 - nb_1
            R12 = np.linalg.norm(R12_vector)
            cos_theta_012 = np.dot(R01_vector, R02_vector) / R01 / R02

            theta = np.arccos(cos_theta_012)
#            term = (1. + gamma * cos_theta_012) ** zeta
            term = (1. + gamma * np.cos( theta - theta_s)) ** zeta
            term *= np.exp(-eta * (R01 ** 2. + R02 ** 2. + R12 ** 2.) /
                           (Rc ** 2.))
            _Rij = dict(Rij=R01)
            _Rik = dict(Rij=R02)
            _Rjk = dict(Rij=R12)
            term *= cutoff_fxn(**_Rij)
            term *= cutoff_fxn(**_Rik)
            term *= cutoff_fxn(**_Rjk)
            fps.append(term)
    return np.array(fps)
"""
def cal_cutoff(Rc, r):
    if r <= Rc:
      return 0.5*(1.0+np.cos(np.pi * r / Rc))
    if r> Rc:
      return 0

def calculate_G4(Rc, cutoff_fxn, theta_s, zeta, eta, gamma):
    fps = []
    thetas =[]
    ymin = -3.0
    xmin = -3.0
    R0 = np.array([0.0,0.0])
    r = Rc/3.0
    #r = 2.66
    nb_1 = np.array([-r, 0.0])
    theta = 0
    for i in range(150):
        R12 = r * np.sin(theta/2)
        if gamma==-1:
          term = (1. + gamma * np.cos(theta - (theta_s-np.pi))) ** zeta
        else:
          term = (1. + gamma * np.cos(theta - theta_s)) ** zeta
        term *= np.exp(-eta * (r ** 2. + r ** 2. + R12 ** 2.) / (Rc ** 2.))
        term *= cal_cutoff(Rc, r)
        term *= cal_cutoff(Rc, r)
        term *= cal_cutoff(Rc, R12)
        term *= 2**(1-zeta)
        fps.append(term)
        thetas.append(theta)
        theta += 0.02
    """
    for i in range(50):
        x = -r + i * 0.02
        print x
        nb_2 = np.array([ x, np.sqrt( r**2- x**2 )])
        R01_vector = nb_1 - R0
        R01 = np.linalg.norm(R01_vector)
        R02_vector = nb_2 - R0
        R02 = np.linalg.norm(R02_vector)
        R12_vector = nb_2 - nb_1
        R12 = np.linalg.norm(R12_vector)
        cos_theta_012 = np.dot(R01_vector, R02_vector) / R01 / R02
        term = (1. + gamma * cos_theta_012) ** zeta
        term *= np.exp(-eta * (r ** 2. + r ** 2. + R12 ** 2.) /
                       (Rc ** 2.))
        _Rij = dict(Rij=r)
        _Rik = dict(Rij=r)
        _Rjk = dict(Rij=R12)
        term *= cutoff_fxn(**_Rij)
        term *= cutoff_fxn(**_Rik)
        term *= cutoff_fxn(**_Rjk)
    """
#    if theta_s == 0.:
#       return np.array(fps)/(2.0*np.sum(np.array(fps))), np.array(thetas)
       
    #return np.array(fps)/np.sum(np.array(fps)), np.array(thetas)
    return np.array(fps), np.array(thetas)

def log(thetas, g2s, gr, etas, filename):
   Rij=0
   Rs=[]
   g_2 = open(filename,'w')
   output = "#     {:12s}".format('r')
   for i in range(len(etas)):
      output += "{:16s}".format(etas[i])
   g_2.write("{:s} {:16s}\n".format(output, "gr"))
   for i in range(len(thetas)):
      output = "{:16.6f}".format(thetas[i])
      for j in range(len(etas)):
         output += "{:16.8f}".format(g2s[etas[j]][i])
      g_2.write("{:s} {:16.8f}\n".format(output, gr[i]))
   return np.array(Rs)  

def print_expectation(g2s, thetas):
    output = ""
    for k in g2s:
       output += "{:12.8f}".format(np.dot(g4s_gr[k], thetas))
    print(output)

default_paras = dict(
    Rc = 3.0,
    zetas = [10.0, 20.0,  40.0,  80.0,  160.0],
    theta_s = [0., 1., 2.],
    gammas = [1,-1],
    eta = 0.005,
    to_normalize_theta_ss = [0.0, 0.15],
    ref_theta_s = 0.45,
    template = True
    )

args = sys.argv
paras=readinputs('input.ini')

for para in default_paras:
   try:
     paras[para] = type_convertion(default_paras[para], paras[para])
   except:
     paras[para] = default_paras[para]

Rc = paras['Rc']
cutoff={'name': 'Cosine', 'kwargs': {'Rc': Rc}}
cutoff_fxn = dict2cutoff(cutoff)

gr_file=open(args[1], 'r')
lines = gr_file.readlines()
gr = []
for i in range(150):
   gr.append(float(lines[i].split()[1]))
gr=np.array(gr)

g4s = {}
g4_orig={}
g4s_gr = {}
theta_ss = [float(theta_s) for theta_s in paras['theta_s']]
zetas = [float(zeta) for zeta in paras['zetas']]
gammas = [float(gamma) for gamma in paras['gammas']]
eta = paras['eta']
keys = []

ref_theta_s = paras['ref_theta_s']
to_normalize_theta_ss = [float(theta_s) for theta_s in paras['to_normalize_theta_ss']]

for zeta_theta in zip(zetas, theta_ss, gammas):
    zeta = zeta_theta[0]
    theta_s = zeta_theta[1]
    gamma = zeta_theta[2]
    key = str(zeta)+'_' + str(theta_s) +'_'+ str(gamma)
    keys.append(key)
    g4s[key],thetas = calculate_G4(Rc, cutoff_fxn, theta_s, zeta, eta, gamma)
    if theta_s == ref_theta_s:
       normal_max = max(g4s[key])

fp_sum = 0

def normalize(fp):
    max_xval = np.argmax(fp)
    if fp[0] > fp[max_xval]/100.0:
       return fp/(np.sum(fp)+np.sum(fp[2*max_xval:]))
    else:
       return fp/np.sum(fp)

for key in g4s.keys():
    
    #if float(key.split('_')[1]) in to_normalize_theta_ss:
    #  max_o = max(g4s[key])
    #  coeff = max_o / normal_max
    #  g4s[key] /= coeff
    if paras['template']:
       g4s[key] = normalize(g4s[key])
    #fp_sum += sum(g4s[key])
    g4s_gr[key] = g4s[key] * gr
    g4s_gr[key] = normalize(g4s_gr[key])
#for key in g4s.keys():
#   g4s[key] /= fp_sum
output_start  = args[1].split('.')[0]+'_'
log(thetas, g4s, gr, keys, output_start+'g4s.dat')
log(thetas, g4s_gr, gr, keys, output_start+'g4s_gr.dat')
print_expectation(g4s_gr, thetas)

