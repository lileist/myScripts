"""Simple test of the Amp calculator, using Gaussian descriptors and neural
network model. Randomly generates data with the EMT potential in MD
simulations."""

import os,sys
from ase import Atoms, Atom, units
import ase.io
import numpy as np
from ase.optimize.lbfgs import LBFGS
from ase.optimize.fire import FIRE
#from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators.emt import EMT
from ase.md.langevin import Langevin
from ase.lattice.surface import fcc110
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet
from ase.constraints import FixAtoms
from expectra.atoms_operator import match
from amp import Amp
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.descriptor.cutoffs import Cosine,Polynomial
#from amp.model.neuralnetwork import NeuralNetwork
from amp.model.tflow import NeuralNetwork
from amp.model import LossFunction
from amp.utilities import TrainingConvergenceError



class NNML(object):
    def __int__(self, atoms, training_set=None, refine_steps=0,
                ml_module=None, lossfunction=None, 
                optimizer=None, optimizer_logfile=None, maxstep=0.1, dt=0.1, dtmax=0.2,):
       """
            Gs = {"H":  [{"type":"G2", "element":"Pd", "eta":20.},
                         {"type":"G2", "element":"Pd", "eta":40.},
                         {"type":"G2", "element":"Pd", "eta":80.},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":4.00, "zeta":1.0},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":10.00, "zeta":1.0},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":20.00, "zeta":1.0},
                        ],
                  "Pd": [{"type":"G2", "element":"H", "eta":20.},
                         {"type":"G2", "element":"H", "eta":40.},
                         {"type":"G2", "element":"H", "eta":80.},
                         {"type":"G2", "element":"Pd", "eta":0.05},
                         {"type":"G2", "element":"Pd", "eta":4.},
                         {"type":"G2", "element":"Pd", "eta":20.},
                         {"type":"G2", "element":"Pd", "eta":40.},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":2.00, "zeta":1.0},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":10.00, "zeta":1.0},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":20.00, "zeta":1.0},
                         {"type":"G4", "elements":["Pd", "Pd"],"eta":0.005, "gamma":40.00, "zeta":1.0},
                         {"type":"G4", "elements":["H", "Pd"],"eta":0.005, "gamma":1.00, "zeta":1.0},
                         {"type":"G4", "elements":["H", "Pd"],"eta":0.005, "gamma":4.00, "zeta":1.0},
                         {"type":"G4", "elements":["H", "Pd"],"eta":0.005, "gamma":20.00, "zeta":1.0},
                        ]}
           
      
           ml_module = Amp(descriptor=Gaussian(Gs=Gs,
                                               cutoff=Cosine(6.0)
                                              #cutoff=Polynomial(gamma=5, Rc=3.0)
                                               ),
                           cores=12,
                           model=NeuralNetwork(hiddenlayers=(50,), activation='sigmoid'))
       """
        self.atoms = atoms
        self.refine_steps = refine_steps
        #TODO: check if parameters are correctly set up
        self.training_set = training_set
        self.ml_module   = ml_module
        self.ml_module.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001,
                                                                     'force_rmse': 0.02},
                                                         force_coefficient=0.2)
        self.optimizer = optimizer
        self.optimizer_logfile = optimizer_logfile
        self.maxstep = maxstep
        self.dt = dt
        self.dtmax = dtmax
        
    def refining(self):
        """
        evaluate energy and force of atoms with pes_refiner
        """
        if self.refine_steps == 0:
           return self.atoms.get_forces()
        self.optimize(atoms=self.atoms, steps=self.refine_steps)
        return self.atoms.get_forces()

    def train(self, step=0):
        """
        training data with machine-learning module given by ml_module
        """
        if step!= 0:
           self.training_set.append(self.atoms)
        while True:
           try:
              #TODO: avoid to calculate fingerprints for thoes already done
              if os.path.exists('amp-untrained-parameters.amp')
                 #load nn model including lossfunction
                 self.ml_module = Amp.load('amp-untrained-parameters.amp')
              self.ml_module.train(self.training_set, overwirte=True)
           except TrainingConvergenceError:
              continue
           break
        return Amp.load('amp.amp')
 
    def run(self, fmax=0.05, steps=100000):
        atoms = self.atoms.copy()
        step = 0
        while step < steps:
            f=self.refine()
            if self.converged(f):
               return True
            #train macnine learning force field
            calc = train(step)
            atoms.set_calculator(calc)
            self.optimize(atoms, fmax, steps)
            step+=1
        return False

    def optimize(self, atoms=None, fmax=0.05,  steps=100000):
        #TODO: define a callable optimizer to avoid the re-construction of optimizer
        if self.optimizer.__name__ == "FIRE":
           opt = self.optimizer(atoms,
                                maxmove = self.maxstep,
                                dt = self.dt, dtmax = self.dtmax,
                                    logfile=self.optimizer_logfile)
        else:
           opt = self.optimizer(atoms,
                                logfile=self.optimizer_logfile,
                                maxstep=self.maxstep)
        print "Geometry optimization is running"
        opt.run(fmax=fmax, steps=steps)
        self.atoms.set_positions(atoms.get_positions())
     """
        #TODO: neb optimizer
        nim = 7  # number of images, including end points
        band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', ss=False)
    
     # to restart, uncomment the following lines which read the previous optimized images into the band
     #    for i in range(1,nim-1):
     #        filename = str(i)+'.CON'
     #        b = read(filename,format='vasp')
     #        band.path[i].set_positions(b.get_positions())
     #        band.path[i].set_cell(b.get_cell())
    
        opt = neb.qm_ssneb(band, maxmove =0.2, dt=0.1)
     #Uncomment to use fire optimization algorithm
     #    opt = neb.fire_ssneb(band, maxmove =0.05, dtmax = 0.03, dt=0.03)
        opt.minimize(forceConverged=0.02, maxIterations = 200)
     """

    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()
        if hasattr(self.atoms, 'get_curvature'):
            return ((forces**2).sum(axis=1).max() < self.fmax**2 and
                    self.atoms.get_curvature() < 0.0)
        return (forces**2).sum(axis=1).max() < self.fmax**2
