import os,sys
import time
import ase.io
import numpy as np
from amp import Amp
from ase.io import Trajectory
from math import sqrt
from amp.utilities import TrainingConvergenceError

def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
           raise

#TODO: optimizer wrapper that wraps all kinds of optimizer and select a defined optimizer in a case/switch manner
#      Need to modify optimize.py in ASE or TSASE to enable run(atoms,fmax, steps):
#      def run(self, atoms=None, fmax, steps):
#          if atoms is not None:
#             self.atoms = atoms
class optimizer(object):
      def __init__(self, atoms, optimizer=None,
                   maxstep=0.1, dt=0.1, dtmax=0.2,
                   trajectory = 'geoopt.traj',
                   logfile='geoopt.log'):
          self.atoms=atoms
          self.optimizer = optimizer
          self.maxstep = maxstep
          self.dt = dt
          self.dtmax = dtmax
          self.trajectory = trajectory
          self.logfile = logfile

      def get_optimizer(self):
          opt = getattr(slef, self.optimizer.__name__, lambda:"Invalid Optimizer")
          return opt

      def neb(self,):
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
          return "neb optimizer"

      def dimer(self,):
          return "neb optimizer"

      def fire(self, ):
          opt = self.optimizer(atoms=self.atoms,
                               maxmove = self.maxstep,
                               dt = self.dt, dtmax = self.dtmax,
                                   trajectory = self.trajectory,
                                   logfile=self.logfile)
          return opt
      def others(self,):
          opt = self.optimizer(atoms=self.atoms,
                               logfile=self.logfile,
                               trajectory = self.trajectory,
                               maxstep=self.maxstep)
          return opt

class NNML(object):
    """
    Optimizer that uses machine-learning methods to get a rough PES
    """

    def __init__(self, atoms, training_set=None, refine_steps=0, presteps=10,logfile=None,
                ml_module=None, lossfunction=None, regressor=None,
                optimizer=None, optimizer_logfile=None, maxstep=0.1, dt=0.1, dtmax=0.2,):
        """
           ml_module = Amp(descriptor=Gaussian(Gs=Gs,
                                               cutoff=Cosine(6.0)
                                              #cutoff=Polynomial(gamma=5, Rc=3.0)
                                               ),
                           cores=12,
                           model=NeuralNetwork(hiddenlayers=(50,), activation='sigmoid'))
        """
        self.atoms = atoms
        self.presteps = presteps
        self.refine_steps = refine_steps
        self.logfile = open(logfile, 'w')
        #TODO: check if parameters are correctly set up
        self.training_set = training_set
        self.ml_module   = ml_module
        if lossfunction is not None:
           self.ml_module.model.lossfunction = lossfunction
        if regressor is not None:
           self.ml_module.model.regressor = regressor
#        self.ml_module.model.lossfunction = LossFunction(convergence={'energy_rmse': 0.001,
#                                                                     'force_rmse': 0.02},
#                                                         force_coefficient=0.2)
        self.optimizer = optimizer
        self.optimizer_logfile = optimizer_logfile
        self.maxstep = maxstep
        self.dt = dt
        self.dtmax = dtmax
        self.cwd = os.getcwd()

    def refine(self):
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
        workdir = self.cwd+'/train'+str(step)
        make_dir(workdir)
        if step!= 0:
           self.training_traj.write(self.atoms)
           #self.training_set.append(self.atoms)
           self.training_set=Trajectory('training.traj','r')
        os.chdir(workdir)
        retries = 0
        while True:
           retries+=1
           print "  Training not converged, retry:", retries
           try:
              #TODO: avoid to calculate fingerprints for thoes already done
              if os.path.exists('amp-untrained-parameters.amp'):
                 #load nn model including lossfunction
                 self.ml_module = Amp.load('amp-untrained-parameters.amp')
              self.ml_module.train(self.training_set, overwrite=True)
           except TrainingConvergenceError:
              continue
           break
        os.chdir(self.cwd)
        return Amp.load(workdir+'/amp.amp')
 
    def run(self, fmax=0.05, steps=100000):
        step = 0
        self.fmax = fmax
        opt_traj = Trajectory('outer_opt.traj', 'w')

        if self.training_set is None:
           self.training_traj = Trajectory('training.traj', 'a', self.atoms)
           #self.training_set = self.optimize(atoms=self.atoms, steps=self.presteps, memo_interval=1)
           self.optimize(atoms=self.atoms, steps=self.presteps, memo_interval=1)
           #self.training_set = [atoms for atoms in Trajectory('geoopt.traj')]
           self.training_set = Trajectory('training.traj','r')
        #sys.exit()
        while step < steps:
            atoms = self.atoms.copy()
            f=self.refine()
            self.log(f, step)
            print "================"
            print "  Refined force:", f
            opt_traj.write(self.atoms)
            if self.converged(f):
               return True
            #train macnine learning force field
            calc = self.train(step)
            atoms.set_calculator(calc)
            print "  Predicted force:", atoms.get_forces()
            self.optimize(atoms, steps)
            print "================"
            step+=1
        return False

    def optimize(self, atoms=None, steps=100000, memo_interval=None):
        #TODO: define a callable optimizer to avoid the re-construction of optimizer
        #print "  Before opt:", atoms.get_positions()
        if self.optimizer.__name__ == "FIRE":
           opt = self.optimizer(atoms,
                                maxmove = self.maxstep,
                                dt = self.dt, dtmax = self.dtmax,
                                    trajectory = 'geoopt.traj',
                                    logfile=self.optimizer_logfile)
        else:
           opt = self.optimizer(atoms,
                                logfile=self.optimizer_logfile,
                                trajectory = 'geoopt.traj',
                                maxstep=self.maxstep)
        print "  Relax geometry with the machine-learning force field"
        if memo_interval is not None:
           #traj = []
           #def traj_memo(atoms=self.atoms):
           #    traj.append(atoms)
               #epot=atoms.get_potential_energy()
               #ekin=atoms.get_kinetic_energy()
           #opt.attach(traj_memo, interval=memo_interval)
           opt.attach(self.training_traj.write, interval=memo_interval)
        opt.run(fmax=self.fmax, steps=steps)
        #print "  After opt:", atoms.get_positions()
        self.atoms.set_positions(atoms.get_positions())
        #if memo_interval is not None:
        #   return traj


    def log(self, forces, step):
        fmax = sqrt((forces**2).sum(axis=1).max())
        e = self.atoms.get_potential_energy()
        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            if step == 0:
                self.logfile.write(
                    '%s  %4s %8s %15s %12s\n' %
                    (' ' * len(name), 'Step', 'Time', 'Energy', 'fmax'))
            self.logfile.write('%s:  %3d %02d:%02d:%02d %15.6f %12.4f\n' %
                               (name, step, T[3], T[4], T[5], e, fmax))
            self.logfile.flush()

    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()
        if hasattr(self.atoms, 'get_curvature'):
            return ((forces**2).sum(axis=1).max() < self.fmax**2 and
                    self.atoms.get_curvature() < 0.0)
        return (forces**2).sum(axis=1).max() < self.fmax**2
