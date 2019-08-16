"""
Opt fingerprints with a given template
"""
import numpy as np
import sys
from amp import Amp
from amp.descriptor.cutoffs import Cosine,dict2cutoff
import ase.io

from amp.descriptor.gaussian import Gaussian

class MyBounds(object):
     def __init__(self, xmax=[2000.0, 8.0], xmin=[0.,0.0] ):
         self.xmax = np.array(xmax)
         self.xmin = np.array(xmin)
     def __call__(self, **kwargs):
         x = kwargs["x_new"]
         tmax = bool(np.all(x <= self.xmax))
         tmin = bool(np.all(x >= self.xmin))
         return tmax and tmin

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

#default_paras = dict(
#    Rc = 3.0,
#    Rmin = 0.0,
#    etas = [10.0, 20.0,  40.0,  80.0,  160.0],
#    Rs =   [0.0,   0.5,   1.0,   2.0,    3.0],
#    gamma = 1,
#    eta = 0.005,
#    )

args = sys.argv
#paras=readinputs('input.ini')
#for para in default_paras:
#   try:
#     paras[para] = type_convertion(default_paras[para], paras[para])
#   except:
#     paras[para] = default_paras[para]

#Rc =  paras['Rc']
#Rmin = paras['Rmin']
Rc = 6.0
Rmin = 1.0
Rmax = 8.0
Rij = 0.0
dr = 0.02
#read in gr function
radium = []
gr = []
gr_file=open(args[1], 'r')
lines = gr_file.readlines()

n_dotall = int((Rmax)/dr)

for i in range(n_dotall):
   fields = lines[i].split()
   if float(fields[0]) >= Rmin and float(fields[1])<= Rmax:
      gr.append(float(fields[1]))
      radium.append(Rij)
   Rij += dr
gr=np.array(gr)
radium = np.array(radium)

#unity = np.full((len(gr),),1.0)
#print unity
#Define cutoff function
cutoff={'name': 'Cosine', 'kwargs': {'Rc': Rc}}
cutoff_fxn = dict2cutoff(cutoff)

step = 0
#Calculate fp reference

def calculate_fp(eta_rs):
    Rij = -0.01
    g2s = []
    n_start = None
    for i in range(n_dotall):
       Rij += 0.02
       args_cutoff_fxn = dict(Rij=Rij)
       #if cutoff['name'] == 'Polynomial':
       #    args_cutoff_fxn['gamma'] = cutoff['kwargs']['gamma']
       term =  np.exp(-eta_rs[0] * ((Rij-eta_rs[1]) ** 2.) / (Rc ** 2.)) * \
                 cutoff_fxn(**args_cutoff_fxn)
       g2s.append(term)
       if Rij >= Rmin and n_start is None:
          n_start = i
    return np.array(g2s)[n_start:-1:1]

def calculate_prime(eta_rs):
    Rij = -0.01
    g2s = []
    n_start = None
    df_deta = []
    df_dRs = []
    for i in range(n_dotall):
       Rij += 0.02
       args_cutoff_fxn = dict(Rij=Rij)
       #if cutoff['name'] == 'Polynomial':
       #    args_cutoff_fxn['gamma'] = cutoff['kwargs']['gamma']
       term =  np.exp(-eta_rs[0] * ((Rij-eta_rs[1]) ** 2.) / (Rc ** 2.)) * \
                 cutoff_fxn(**args_cutoff_fxn)
       df_deta.append(term * (-( Rij - eta_rs[1]) ** 2. / (Rc ** 2.)))
       df_dRs.append(term * 2 * eta_rs[0] * (Rij-eta_rs[1]) / (Rc ** 2.))
       if Rij >= Rmin and n_start is None:
          n_start = i
    return  np.array(df_deta)[n_start:-1:1], np.array(df_dRs)[n_start:-1:1]

#Calculate Refrence fps
#ref_etas = [float(field) for field in paras['etas']]
#ref_Rs = [float(field) for field in paras['Rs']]
#fp_refs = {}
eta = float(args[2])
Rs  = float(args[3])
#for eta_Rs in zip(ref_etas, ref_Rs]:
fp_ref = calculate_fp([eta, Rs])
fp_ref /= np.sum(fp_ref)

def log(g2s, gr, filename='selected_fp.dat'):
    global radium
    g_2 = open(filename,'w')
    for i in range(len(g2s)):
       g_2.write("{:12.8f} {:12.8f} {:12.8f}\n".format(radium[i], gr[i], g2s[i]))
    g_2.close()
    return
   
def get_loss(parametervector):
    global step
    global Rs
    #global fp_ref
    #numb_fps = len(parametervector) / 2
    #numb_fps = 1
    #fp = calculate_fp(parametervector.reshape(numb_fps,2))
    fp = calculate_fp(parametervector)

    fp_gr = fp * gr
    #df_deta, df_dRs = calculate_prime(parametervector)
    #print fp_gr
    #normalize
    fp_gr /= np.sum(fp_gr)
    log(fp_gr, gr, filename=str(Rs)+'fp.dat')
    diff = np.sum(np.absolute(fp_gr - fp_ref))
    #print "{:6d} {:14.6f} {:14.6f} {:14.6f} ".format(step, diff, min(parametervector), 
    #      max(parametervector))
    #      -np.dot(df_deta, gr), -np.dot(df_dRs, gr))
    #print parametervector
    step += 1
    #return diff, np.array([np.sum(np.absolute(df_deta*gr)), np.sum(np.absolute(df_dRs*gr))])
    return diff

    

def calculate_fps(eta_rss):
    global Rmin
    g2s = {}
    g2_orig={}
    keys=[]
    #print cutoff_fxn
    #fp_sum = 0
    n_start =None
    n_end =None
    for eta_rs in eta_rss:
       Rij=-0.01
       key = str(eta_rs[0])+'_'+str(eta_rs[1])
       keys.append(key)
       g2s[key], n_start = calculate_fp(eta_rs)
       #normalize fps independently
       #g2s[key] = g2_orig[key] / sum(g2_orig[key])
       g2s[key] = g2s[key][n_start:-1:1]
       g2s[key] = g2s[key] / np.sum(g2s[key])
       #normalize fps overall
       #fp_sum += sum(g2_orig[key])
    #for key in g2s.keys():
    #   g2s[key] = g2_orig[key] / fp_sum
    return g2s

def print_expectation(g2s, thetas):
    output = ""
    for k in g2s:
       output += "{:12.8f}".format(np.dot(g2s[k], thetas))
    print output

log_min = open('bh_min'+str(Rs)+'.dat','w')
log_min.write("xValues    fValue     nAccepted\n")
def print_fun(x, f, accepted):
    global log_bh
    vals = ""
    for val in x:
       vals += "{:12.8f}".format(val)
    vals+="{:12.8f}".format(f) 
    vals+="{:6d}".format(accepted) 
    log_min.write("%s\n" % (vals))


def main():
  #Initialize parameters
  #Rs   = [0., 0.5, 1.0,  1.5,  0.0, 0., 0.5]
  global eta
  global Rs
  #etas   = [eta]
  #Rss = [Rs]
  etas = [100.]
  Rss = [0.]
  #g2_ref = calculate_fp([eta, Rs])
  
  parametervector = []
  for i in range(len(Rss)):
     parametervector.append(etas[i])
     parametervector.append(Rss[i])
 
  #x0 = np.array(parametervector)
  x0 = parametervector
 
 
  optimizer = 'basinhopping'
  #optimizer = 'L-BFGS-B'
 
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
     minimizer_kwargs = {"method": "L-BFGS-B",
                         'options': {
                                     'maxiter': 10,
                                     #'maxfun':1
                                     }
                        }
     #mybounds = MyBounds()
     results = optimizer(get_loss, x0, 
                         #minimizer_kwargs=minimizer_kwargs,
                         niter=500,
                         T = 0.01,
                         stepsize = 1.0,
                         #niter_success=1000,
                         disp=True,
                         callback=print_fun,
      #                   accept_test=mybounds
                         )
     print results
     sys.exit()
  elif optimizer == 'L-BFGS-B':
     from scipy.optimize import minimize as optimizer
     optimizer_kwargs = {
                         'method': 'L-BFGS-B',
                         'options': {'ftol': 1e-05,
                                     'gtol': 1e-18,
                                     'maxfun': 1000000,
                                     'maxiter':1500000,
                                     'disp': True}
                        }
     import scipy
     from distutils.version import StrictVersion
     if StrictVersion(scipy.__version__) >= StrictVersion('0.17.0'):
         optimizer_kwargs['options']['maxls'] = 2000
  optimizer(get_loss, x0, jac=False, **optimizer_kwargs)
  #try:
  #except:
  #   print "Optimizer error"
  #output_start = args[1].split('.')[0]+'_'
  #log(radium, g2s, gr, keys, output_start+'g2s.dat')
  #log(radium, g2s_gr, gr, keys, output_start+'g2s_gr.dat')

  #print_expectation(g2s_gr, radium)

if __name__ == '__main__':
    main()

