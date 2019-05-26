"""
This neural-network machine-learning optimizer for NEB.
Initially contributed by Lei Li
Other contributors:
"""
import os,sys
import time, copy
import ase.io
import numpy as np
from amp import Amp
from ase.io import Trajectory
from math import sqrt
from amp.utilities import TrainingConvergenceError
from scipy.optimize import minimize
from ase.optimize.optimize import Optimizer
from ase.optimize import FIRE
from ase.calculators.lj import LennardJones
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError
from tsase import neb
from tsase.neb.minimizer_ssneb import minimizer_ssneb

def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
           raise

class pseudoCalculator(Calculator):
    implemented_properties = ['energy', 'forces']
    default_parameters = {}
    nolabel = True

    def __init__(self, energy=None, forces=None,**kwargs):
        Calculator.__init__(self, **kwargs)
        self.energy = energy
        self.forces = forces

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        self.results['energy'] = self.energy
        self.results['forces'] = self.forces

class NNML(minimizer_ssneb):
    """
    Optimizer that uses machine-learning methods to get a rough PES
    """

    def __init__(self, band, restart=None, logfile='-', trajectory=None,
                 ml_module=None, max_training_cycle=10, lossfunction=None, regressor=None,
                 optimizer=FIRE, optimizer_logfile='ml_opt.log', maxstep=0.1, dt=0.1, dtmax=0.2,
                 force_consistent=None
                ):
        """
           ml_module = Amp(descriptor=Gaussian(Gs=Gs,
                                               cutoff=Cosine(6.0)
                                              #cutoff=Polynomial(gamma=5, Rc=3.0)
                                               ),
                           cores=12,
                           model=NeuralNetwork(hiddenlayers=(50,), activation='sigmoid'))
        """
        minimizer_ssneb.__init__(self, band)

        #duplicated path to store amp calculator
        self.band_ml = copy.deepcopy(band)
        self.nimages = self.band.numImages
        #use a pre-defined approxPot to avoid unphysical structure
        #self.approxPot = LennardJones(epsilon=0.65, sigma=2.744)
        #self.approxPot_replica = atoms.copy()
        #self.approxPot_replica.set_calculator(self.approxPot)

        #self.logfile = open(logfile, 'w')
        #TODO: check if parameters are correctly set up
        self.ml_module   = ml_module
        if lossfunction is not None:
           self.ml_module.model.lossfunction = lossfunction
        if regressor is not None:
           self.ml_module.model.regressor = regressor

        #Optimizer used to relax geometry on ML PES
        self.optimizer = optimizer
        self.optimizer_logfile = optimizer_logfile
        self.progress_log = open('progress.log','w')
        self.maxstep = maxstep
        self.dt = dt
        self.dtmax = dtmax
        self.force_consistent = force_consistent
        self.training_set = []
        self.training_traj = Trajectory('training.traj','w')
        self.ml_e = None
        self.ml_log = open('ml_opt.log', 'w')
     
        self.train_endpoints = True
        self.calc_endpoints = True
        self.cwd = os.getcwd()
        self.function_calls =0
        self.force_calls = 0

    def relax_model(self, r0s):
        """
        Minimization on ml PES
        """
        if not os.path.exists(self.cwd + '/amp_neb'):
           make_dir(self.cwd + '/amp_neb')
        os.chdir(self.cwd + '/amp_neb')

        for i in range(self.nimages):
            # making a directory for each image, which is nessecary for vasp to read last step's WAVECAR
            # also, it is good to prevent overwriting files for parallelizaiton over images
            fdname = '0'+str(i)
            if not os.path.exists(fdname): os.mkdir(fdname)
        
        if self.calc_endpoints:
           os.chdir(self.cwd + '/amp_neb/00')
           self.band_ml.path[0].get_potential_energy()
           os.chdir(self.cwd + '/amp_neb/0'+str(self.nimages-1))
           self.band_ml.path[-1].get_potential_energy()
           self.calc_endpoints = False
           os.chdir(self.cwd + '/amp_neb')

        opt = neb.fire_ssneb(self.band_ml, maxmove =self.maxstep, dtmax = self.dtmax, dt=self.dt)
        self.progress_log.write("  Relax geometry with the machine-learning force field\n")
        for i in range(self.nimages-2):
           self.band_ml.path[i+1].set_positions(r0s[i])
        opt.minimize(forceConverged=0.10, maxIterations = 400)

        r1s = np.array([self.band_ml.path[i+1].get_positions() for i in range(self.nimages-2)])
       
        #if np.all((r1s-r0s) < 0.002 ):
        #   print "atoms not moved"
        os.chdir(self.cwd)
        return r1s

    def update(self, rs, es, fs):
        """
        training data with machine-learning module given by ml_module
        """
        print "Start to train"
        if self.train_endpoints:
           print "End points used for training"
           nimages = self.nimages
           n_offset = 0
           self.train_endpoints = False
        else:
           nimages = self.nimages - 2
           n_offset = 1

        for i in range(nimages):
           self.band.path[i+n_offset].set_positions(rs[i])
        #self.approxPot_replica.set_positions(r)
        #f = f - self.approxPot_replica.get_forces()
        #e = e - self.approxPot_replica.get_potential_energy()
          
           pseudoAtoms = self.band.path[i+n_offset].copy()
           pseudoAtoms.set_calculator(pseudoCalculator(energy= es[i], \
                                                       forces= fs[i]))
           pseudoAtoms.get_potential_energy()
           pseudoAtoms.get_forces()
           self.training_set.append(pseudoAtoms)
           self.training_traj.write(pseudoAtoms)

        #self.training_set=Trajectory('training.traj','r')
        #os.chdir(workdir)
        if os.path.exists('amp-fingerprint-primes.ampdb'):
           os.system('rm -rf amp-fingerprint-primes.ampdb')
        if os.path.exists('amp-fingerprints.ampdb'):
           os.system('rm -rf amp-fingerprints.ampdb')
        if os.path.exists('amp-fingerprints.ampdb'):
           os.system('rm -rf amp-neighborlists.ampdb')
        #if os.path.exists('amp.amp'):
        #   os.system('rm amp.amp')
        #if os.path.exists('amp-untrained-parameters.amp'):
           #load nn model including lossfunction
        #   os.system('rm amp-untrained-parameters.amp')
        try:
           self.progress_log.write("Train ml model\n")
           self.ml_module.train(images='training.traj', overwrite=True)
        except TrainingConvergenceError:
           os.system('mv amp-untrained-parameters.amp amp.amp')
           pass
        #load ml model
        try:
           ml_calc = Amp.load('amp.amp')
        except:
           ml_calc = Amp.load('amp-untrained-parameters.amp')
           pass
        
        for i in range(self.nimages):
            self.band_ml.path[i].set_calculator(ml_calc)

    def get_path_ref(self):
        #get positions for intermidate images
        start_time = time.time()
        print "band.forces"
        self.band.forces()
        end_time = time.time()
        print "band forces time:", end_time-start_time

        if self.train_endpoints:
           rs = [self.band.path[i].get_positions()\
                        for i in range(self.nimages)]
           es = [self.band.path[i].u \
                for i in range(self.nimages)]
           fs = [self.band.path[i].f \
                 for i in range(self.nimages)]
        else:
           rs = [self.band.path[i+1].get_positions()\
                        for i in range(self.nimages-2)]
           es = [self.band.path[i+1].u \
                for i in range(self.nimages-2)]
           fs = [self.band.path[i+1].f \
                 for i in range(self.nimages-2)]

        print "get fe time:", time.time()-end_time
        return np.array(rs), np.array(es), np.array(fs)

    def set_path_rs(self, rs):
        for i in range(self.nimages-2):
           self.band.path[i+1].set_positions(rs[i])
    
    def step(self):
        print "1st cycle:"   
        r0s, e0s, f0s = self.get_path_ref()
        #update ml model
        self.update(r0s, e0s, f0s)
     
        #relax atoms on ml-rough PES
        r1s = self.relax_model(r0s)

        self.set_path_rs(r1s)

        r1s, e1s, f1s = self.get_path_ref()

        self.function_calls += 1
        self.force_calls += 1

        count = 0
        self.progress_log.write("# New step started:\n")
        self.progress_log.flush()
        while np.all(e1s >= e0s):

            self.update(r1s, e1s, f1s)
            r1s = self.relax_model(r0s)

            self.set_path_rs(r1s)
            r1s, e1s, f1s = self.get_path_ref()

            self.function_calls += 1
            self.force_calls += 1
            self.progress_log.write("  Opted with ML: {:3d} {:16.8f} {:16.8f}\n".format(count, e0, e1))
            self.progress_log.flush()

            count += 1
            if count == 30:
                raise RuntimeError('A descent model could not be built')

