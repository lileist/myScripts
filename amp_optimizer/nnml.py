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


class NNML(object):
    def __int__(self):
    

"""
Gs = {"H": [{"type":"G2", "element":"Pd", "eta":20.},
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


calc = Amp(descriptor=Gaussian(Gs=Gs,
                               cutoff=Cosine(6.0)
                               #cutoff=Polynomial(gamma=5, Rc=3.0)
                               ),
           cores=12,
           model=NeuralNetwork(hiddenlayers=(50,), activation='sigmoid'))
"""
calc = Amp.load('amp-untrained-parameters.amp')
os.system('rm amp-untrained-parameters.amp')
calc.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001,
                                                    'force_rmse': 0.02},
                                        force_coefficient=0.2)

calc.train(images='train.traj')

