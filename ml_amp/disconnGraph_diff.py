"""
Evaluate AMP force field by comparing to the PES obtained from LJ
"""
import numpy as np
import logging
from pele.optimize import lbfgs_scipy, fire, lbfgs_cpp, mylbfgs
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
import sys
#import pele.potentials.ase_ff as ase_ff
from ase import Atoms
from ase.calculators.lj import LennardJones
from amp import Amp
from pele.transition_states import orthogopt
import ase.io
from pele.storage import Database
#from expectra.atoms_operator import match


class compareStructure(object):
#     def __int__(self, min1, min2):
#        self.min1 = min1
#        self.min2 = min2
     def __call__(self, min1, min2):
        """
        p= ase.io.read('POSCAR')
        p1 = p.copy()
        p1.set_positions(min1.coords.reshape((7, -1)))
        p2 = p.copy()
        p2.set_positions(min2.coords.reshape((7, -1)))
        if match(p1, p2, 0.2, 1.3, True, False):
           return True
        """
        return False

def get_rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def fetch_pes(db):
    min_energies = []
    min_id =[]
    ts_energies = []
    ts_id = []
    barrier =[]
    reactionE = []
    for m1 in db.minima(): 
#       print m1.id()
       min_energies.append(m1.energy)
       min_id.append(m1.id)
       for m2 in db.minima():
           if m1==m2 or m1.id()>m2.id():
             continue
           try:
            saddleE = db.getTransitionState(m1, m2).energy
            ts_energies.append(saddleE)
            barrier.append(saddleE - m1.energy)
            reactionE.append(m2.energy - m1.energy)
            ts_id.append([m1.id(),m2.id()])
           except:
            pass
    return np.array(min_energies), np.array(ts_energies), np.array(ts_id), np.array(min_id), np.array(barrier),np.array(reactionE)

def calculate_spointE(atoms, db):
    min_energies = []
    ts_energies = []
    ts_id = []
    """
    for m1 in db.minima():
       pot = ase_ff(atoms.set_positions(m1.coords))
       res1 = mylbfgs(m1.coords, pot)
       E1 = res1.energy
       min_energies.append(E1)

    for ts in db.transition_states:
       pot = ase_ff(atoms.set_positions(ts.coords))
       res1 = mylbfgs(m1.coords, pot, nsteps=0)
       E1 = res1.energy
       ts_energies.append(E1)
       ts_id.append([ts.min1, ts.min2])
    """

    return np.array(min_energies), np.array(ts_energies), np.array(ts_id)

args = sys.argv

db_amp = Database(db=args[1],accuracy=1e-10)
db_lj = Database(db=args[2], accuracy=1e-10)
#db = Database(db='amp_disconnGraph.db',accuracy=1e-3, compareMinima= compareStructure())

min_amp, ts_amp, tsid_amp, minamp_id, barrier_amp, reactionE_amp = fetch_pes(db_amp)
min_lj, ts_lj, tsid_lj, minlj_id, barrier_lj, reactionE_lj = fetch_pes(db_lj)

print "# of min (lj, amp):",len(min_lj), len(min_amp)
print "# of TS  (lj, amp):", len(ts_lj), len(ts_amp)

print "Connected minsID (lj, amp) :"
for i in range(max(len(tsid_lj), len(tsid_amp))):
   if len(tsid_amp) > i and len(tsid_lj) > i:
      print "   ",tsid_lj[i], ts_lj[i], tsid_amp[i], ts_amp[i]
      continue
   if len(tsid_amp) <= i:
      print "    ",tsid_lj[i], ts_lj[i], None, None
      continue
   print "   ",None, None, tsid_amp[i], ts_amp[i]

if len(min_amp)!=len(min_lj):
   print "Non-equal number of minimum. Applying AMP pot on db_lj"
   calc = Amp.load('amp.amp')
   p1 = read('pos.xyz',index=0,format='xyz')
   p1.set_calculator(calc)
   min_amp, ts_amp, tsid_amp = calculate_spointE(p1, db_lj)

   for i in range(max(len(tsid_lj), len(tsid_amp))):
      if len(tsid_amp) > i and len(tsid_lj) > i:
         print "   ",tsid_lj[i], ts_lj[i], tsid_amp[i], ts_amp[i]
         continue
      if len(tsid_amp) <= i:
         print "    ",tsid_lj[i], ts_lj[i], None, None
         continue
      print "   ",None, None, tsid_amp[i], ts_amp[i]

min_lj = np.append(min_lj, ts_lj)
min_amp = np.append(min_amp, ts_amp)

relativeE_amp = np.append(barrier_amp, reactionE_amp)
relativeE_lj = np.append(barrier_lj, reactionE_lj)

mae = ((min_lj - min_amp)**2).mean(axis=0)
rmse = get_rmse(min_lj,  min_amp)
print "absolute E mae, rmse:"
print "    ",mae, rmse
mae = ((relativeE_amp - relativeE_lj)**2).mean(axis=0)
rmse = get_rmse(relativeE_amp,  relativeE_lj)
print "barrier and reaction energy mae, rmse:"
print "    ",mae, rmse
