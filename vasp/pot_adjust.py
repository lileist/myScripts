#!/usr/bin/env python
"""
This code is used to find NELECT to match electrode potential to specific value relative to SHE (4.43 eV).
Method: calculate Ef (Fermi Energy relative to vaccum level), then U (electrode potential = (-Ef-SHE)/e)
Steepest desent algorithm used to minimize dU (U-U_target)
An example of inputfile (any line starts with '#' will be ignored):
    cp2k_inp = suppl.inp
    job_submit_script = qsub.hf
    job_submit_cmd    = sbatch
"""
import sys
import os
import errno
import ase
from ase import Atoms
import numpy
import subprocess

sys.stdout.flush()

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def vaccum_pot(idir, start_index, end_index):
#  spin = int(arg[2])
  pot_input = open('LOCPOT', 'r')
  while True:
      line = pot_input.readline()
      if not line:
         break
      for i in range(6):
          line = pot_input.readline()
      ions = [int(field) for field in line.split()]
      print "ions:", ions
      n_ions = sum(ions)
      print "total number of ions:", n_ions
      line = pot_input.readline()
      for i in range(n_ions):
          pot_input.readline()
      pot_input.readline()
      vlocal = []
      vlocal.append(0.0)
#      for i in range(spin):
      fields = [int(field) for field in pot_input.readline().split()]
      ngx = fields[0]
      ngy = fields[1]
      ngz = fields[2]
      print "ngx ngy ngz:",fields 
      nplwv = ngx*ngy*ngz
      #first element is not used
      for j in xrange(int(nplwv/5)+1):
          vlocal.extend([float(field) for field in pot_input.readline().split()])
      print len(vlocal)
      if len(vlocal)>=nplwv+1:
         break
  if idir == 'x':
     nout = ngx
  elif idir == 'y':
     nout = ngy
  else:
     nout = ngz
  scale = 1./float(nplwv/nout)
  vav = []
  vav.append(0.0)
  for i in xrange(nout):
     vav.append(0.0)
  if idir == 'x':
    for i in xrange(1, ngx+1):
        for j in xrange(1, ngy+1):
            for k in xrange(1, ngz+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[i] = vav[i]+ vlocal[ipl] * scale
  elif idir == 'y':
    for j in xrange(1, ngy+1):
        for k in xrange(1, ngz+1):
            for i in xrange(1, ngx+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[j] = vav[j]+ vlocal[ipl] * scale
  elif idir == 'z':
    for k in xrange(1, ngz+1):
        for j in xrange(1, ngy+1):
            for i in xrange(1, ngx+1):
               ipl = i + ( (j-1) + (k-1) * ngy) * ngx
               vav[k] = vav[k]+ vlocal[ipl] * scale
  output = open('loc_pot.dat','w')
  for i in xrange(1, len(vav)):
      output.write("%d %15.6f\n"%(i, vav[i]))
  output.close()
  return numpy.mean(vav[start_index:end_index+1])

def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
           raise

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

def read_old_u_log(filename, target_pot):
    f = open(filename, 'r')
    lines = f.readlines()
    nelect = []
    pot = []
    for line in lines:
        if line.startswith('step'):
           continue
        fields = line.split()
        nelect.append(float(fields[1]))
        pot.append(float(fields[2]))
    fit_paras = numpy.polyfit(nelect, pot, 1)
    k = fit_paras[0]
    b = fit_paras[1]
    return nelect, pot, (float(target_pot) - b)/k

def main():
    arg = sys.argv
    paras = readinputs(arg[1])
    inp = open(paras['incar'], 'r')
    lines = inp.readlines()
    #distances=paras['distance'].split()
    n_step = 0
    if 'nelect' in paras:
       nelect_new = float(paras['nelect'])
    else:
       nelect_new = None
    if 'guess' in paras:
       guess = int(paras['guess'])
    else:
       guess = None
    if 'linearFitting' in paras:
       linearFitting = int(paras['linearFitting'])
    else:
       linearFitting = None
    step_size = float(paras['step_size'])
    max_step = int(paras['max_step'])
    vaspsol = int(paras['vaspsol'])
    pot_she = float(paras['pot_SHE'])
    u_target = float(paras['applied_voltage'])
    convergence = float(paras['convergence'])
    u_history = int(paras['u_history'])
    du_limit = float(paras['du_limit'])
    k=0
    nelect_list = []
    du_list = []

    #used to restart calculation with history log file
    if u_history == 1:
       nelect_list, du_list, nelect_new = read_old_u_log('u_log.dat', u_target) 

    #structure optimization with given parameters in INCAR
    if 'geo_opt' in paras:
       print "##### vasp running"
       os.system(paras['run_vasp'])
       proc = subprocess.Popen("grep 'reached required accuracy' OUTCAR|tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
       vasp_done = proc.communicate()
       print "vasp run infor",vasp_done[0]
       proc_2 = subprocess.Popen("grep 'energy  without entropy' OUTCAR |tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
       vasp_output = proc_2.communicate()
       energy = vasp_output[0].split()[6]
       if "stopping structural energy minimisation" not in vasp_done[0]:
          print "job with nelect=", nelect_new, "is not converged"
          sys.exit()
       os.system("cp CONTCAR POSCAR")

    log_u = open('u_log.dat','w')
    log_u.write("{:5s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}{:22s}\n".format("step", "nelect_new", "vacc_level","fermi","u_new", "du_new", "k", "energy"))

    for n_step in range(max_step):

        if n_step == 1:
           # fermi level too high. incease e, fermi shift to lower level
           if du_new < 0: 
              nelect_new = nelect_new - step_size
           else:
              nelect_new = nelect_new + step_size

        #prepare INCAR
        output = open('INCAR', 'w')
        for line in lines:
            if 'NELECT' in line and nelect_new is not None: 
               output.write("%s %15.6f\n"%('NELECT = ', nelect_new))
               continue
            if 'ISTART'  in line:
               output.write("%s %d\n"%('ISTART = ', 1))
               continue
            if 'ICHARG'  in line:
               output.write("%s %d\n"%('ICHARG = ', 1))
               continue
            if 'NSW'  in line:
               output.write("%s %d\n"%('NSW = ', 0))
               continue
            if 'IBRION'  in line:
               output.write("%s %d\n"%('IBRION = ', -1))
               continue
            if 'NELMIN'  in line:
               output.write("%s %d\n"%('NELMIN = ', 2))
               continue
            output.write("%s"%(line))
        if vaspsol == 1:
           output.write("%s\n"%(' LSOL = .TRUE.'))
           output.write("%s\n"%(' EB_K = 78.4'))
           #output.write("%s\n"%(' TAU = 0.00'))
           output.write("%s\n"%(' LAMBDA_D_K = 3.0'))
           output.write("%s\n"%(' NC_K    =  0.0047300'))
        output.close()

        #run vasp
        print "##### vasp running"
        os.system(paras['run_vasp'])
        proc = subprocess.Popen("grep 'aborting loop because EDIFF is reached' OUTCAR|tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        vasp_done = proc.communicate()
        print "vasp run infor",vasp_done[0]
        proc_2 = subprocess.Popen("grep 'energy  without entropy' OUTCAR |tail -n 1", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        vasp_output = proc_2.communicate()
        energy = vasp_output[0].split()[6]
        if "EDIFF is reached" not in vasp_done[0]:
           print "job with nelect=", nelect_new, "is not converged"
           sys.exit()

        #read nelect from OUTCAR
        if nelect_new is None:
           proc_1 = subprocess.Popen("grep NELECT OUTCAR", shell=True, stdout=subprocess.PIPE)
           out = proc_1.communicate()
           nelect_new = float(out[0].split()[2])
           os.system("cp CONTCAR POSCAR")
           print "POSCAR updated"

        #calculate voltage
        doscar = open('DOSCAR', 'r')
        for i in range(6):
            line = doscar.readline()
        fermi = float(line.split()[3])
        doscar.close()
        vaccum_level = vaccum_pot(paras['idir'], int(paras['start_index']), int(paras['end_index']))
        print "vaccum level:", vaccum_level
        u_new = vaccum_level - fermi - pot_she
        #when u < u_target, take abs. since need du --> 0. if no abs, du --> negative infinite
        du_new = u_new - u_target
        
        nelect_list.append(nelect_new)
        du_list.append(du_new)

        log_u.write("%d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %s\n"%(n_step, nelect_new, vaccum_level, fermi, u_new, du_new, k, energy))
        log_u.flush()

        du_new = abs(du_new)
        if du_new < convergence:
           break
        #use small step size if close to target
        if du_new < du_limit and n_step > 0:
           step_size = float(paras['small_step_size'])
           print 'step size', step_size
        if n_step ==0:
           du_min = du_new
           nelect_min = nelect_new
           du_old = du_new
           nelect_old = nelect_new
           os.system("cp CONTCAR POSCAR")
           continue

        if du_new < du_min:
           du_min = du_new
           nelect_min = nelect_new

        k = (du_new - du_old)/(nelect_new - nelect_old)
        #nelect_new shows a linear relationship with nelect_new. guess 'nelect_new' by constructing a straight line
        if linearFitting:
           if n_step % linearFitting == 0 and len(nelect_list)>2:
              fit_paras = numpy.polyfit(nelect_list[-3:], du_list[-3:], 1)
              k = fit_paras[0]
              b = fit_paras[1]
              du_old = du_new
              nelect_old = nelect_new
              nelect_new = - b / k
              continue
        #guess target nelect based on most recent two dots
        if guess:
           if n_step == guess:
              b = du_new - k * nelect_new
              du_old = du_new
              nelect_old = nelect_new
              nelect_new = - b / k
              continue
        du_old = du_new
        nelect_old = nelect_new
        nelect_new = nelect_min - step_size * k
if __name__ == '__main__':
    main()
    
