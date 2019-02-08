"""
Evaluate AMP force field by comparing to the PES obtained from LJ
"""
import numpy as np
import logging
from pele.optimize import lbfgs_scipy, fire, lbfgs_cpp, mylbfgs
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
import sys, os
import pele.potentials.ase_ff as ase_ff
from pele.potentials.lj import LJ
from ase import Atoms
from ase.calculators.lj import LennardJones
from amp import Amp
from ase.io import read
from pele.storage import Database
from expectra.atoms_operator import match


class compareStructure(object):
#     def __int__(self, min1, min2):
#        self.min1 = min1
#        self.min2 = min2

     def __call__(self, min1, min2):
        p= ase.io.read('POSCAR')
        p1 = p.copy()
        p1.set_positions(min1.coords.reshape((7, -1)))
        p2 = p.copy()
        p2.set_positions(min2.coords.reshape((7, -1)))
        if match(p1, p2, 0.2, 1.3, True, False):
           return True
        return False

def get_rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def fetch_pes(db):
    min_energies = []
    ts_energies = []
    ts_id = []
    for m1 in db.minima(): 
#       print m1.id()
       min_energies.append(m1.energy)
       for m2 in db.minima():
           if m1==m2 or m1.id()>m2.id():
             continue
           try:
            ts_energies.append(db.getTransitionState(m1, m2).energy)
            ts_id.append([m1.id(),m2.id()])
           except:
            pass
    return np.array(min_energies), np.array(ts_energies), np.array(ts_id)

def calculate_spointE(pot, db, nsteps=1000):
    min_energies = []
    ts_energies = []
    ts_id = []
    for m1 in db.minima():
       #pot = ase_ff(atoms.set_positions(m1.coords.reshape((7, -1))))
       res1 = mylbfgs(m1.coords, pot, nsteps=nsteps)
       E1 = res1.energy
       min_energies.append(E1)
       for m2 in db.minima():
           if m1==m2 or m1.id()>m2.id():
             continue
           try:
            ts_coords = db.getTransitionState(m1, m2).coords
            res1 = mylbfgs(ts_coords, pot, nsteps=0)
            E1 = res1.energy
            ts_energies.append(E1)
            ts_id.append([m1.id(), m2.id()])
           except:
            pass
    return np.array(min_energies), np.array(ts_energies), np.array(ts_id)

def print_ts(tsid_lj,ts_lj, tsid_amp, ts_amp):
   for i in range(max(len(tsid_lj), len(tsid_amp))):
      if len(tsid_amp) > i and len(tsid_lj) > i:
         print "     ",tsid_lj[i], ts_lj[i], tsid_amp[i], ts_amp[i]
         continue
      if len(tsid_amp) <= i:
         print "      ",tsid_lj[i], ts_lj[i], None, None
         continue
      print "     ",None, None, tsid_amp[i], ts_amp[i]

def print_min(min_lj, min_amp):
   for i in range(max(len(min_lj), len(min_amp))):
      if len(min_amp) > i and len(min_lj) > i:
         print "     ",min_lj[i], min_amp[i]
         continue
      if len(min_amp) <= i:
         print "      ",min_lj[i], None
         continue
      print "     ",None, min_amp[i]

args = sys.argv

try:
  os.system("rm -rf amp-* opt.traj")
except:
  pass

wdir = os.getcwd()

amp_dir = wdir+"/../"+args[1]

os.chdir(amp_dir)
print "AMP directory:", os.getcwd()
db_amp = Database(db='amp_disconnGraph.db',accuracy=1e-5)
calc = Amp.load('amp.amp')
os.chdir(wdir)

db_lj = Database(db='lj.db', accuracy=1e-5)
#db = Database(db='amp_disconnGraph.db',accuracy=1e-3, compareMinima= compareStructure())

min_amp, ts_amp, tsid_amp = fetch_pes(db_amp)
min_lj, ts_lj, tsid_lj = fetch_pes(db_lj)

print "# of min (lj, amp):",len(min_lj), len(min_amp)
print "# of TS  (lj, amp):", len(ts_lj), len(ts_amp)

print_min(min_lj, min_amp)

print "Connected minsID (lj, amp) :"
print_ts(tsid_lj, ts_lj, tsid_amp, ts_amp)

flag = 0
if len(min_amp)!=len(min_lj) or len(ts_amp)!=len(ts_lj):
   print "#Non-equal number of minimum."
   print "==Applying AMP pot on lj particles"
   p1 = read('pos.xyz',index=0,format='xyz')
   p1.set_calculator(calc)
   pot = ase_ff(p1)
   min_amp_1, ts_amp_1, tsid_amp_1 = calculate_spointE(pot, db_lj)
   flag = 1
   print_min(min_lj, min_amp_1)
   print_ts(tsid_lj, ts_lj, tsid_amp_1, ts_amp_1)

   min_lj = np.append(min_lj, ts_lj)
   min_amp_1 = np.append(min_amp_1, ts_amp_1)
   mae_1 = ((min_lj - min_amp_1)**2).mean(axis=0)
   rmse_1 = get_rmse(min_lj,  min_amp_1)

   print "==Applying lj pot on AMP particles"
   pot = LJ()
   min_lj_1, ts_lj_1, tsid_lj_1 = calculate_spointE(pot, db_amp,nsteps=0)
   print_min(min_lj_1, min_amp)
   print_ts(tsid_lj_1, ts_lj_1, tsid_amp, ts_amp)

   min_lj_1 = np.append(min_lj_1, ts_lj_1)
   min_amp = np.append(min_amp, ts_amp)
   mae = ((min_lj_1 - min_amp)**2).mean(axis=0)
   rmse = get_rmse(min_lj_1,  min_amp)
else:
   mae_1=0
   rmse_1=0
   min_lj = np.append(min_lj, ts_lj)
   min_amp = np.append(min_amp, ts_amp)
   mae = ((min_lj - min_amp)**2).mean(axis=0)
   rmse = get_rmse(min_lj,  min_amp)
   
print "mae, rmse, non-equal (lj, amp):"
print '{:12.8f} {:12.8f} {:12.8f} {:12.8f} {:d}'.format(mae_1, rmse_1, mae, rmse, flag)
