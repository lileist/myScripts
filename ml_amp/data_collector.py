"""
This code is designed to collect structure-energy-force information from geometry optimization based on
an existing trajectory file. (Contributed by Lei Li)
Example inputs (config.ini):
[main]
structurefile = 1000K.xyz
fileformat = xyz
structure_slice = 1200:5100:100
calc_type = vasp

[calculator]
xc      = PBE
prec    = Medium
istart  = 1
icharg  = 1
ispin   = 1
encut   = 300
ismear  = 0
sigma   = 0.1
nelm    = 200
nelmin  = 4
ediff   = 1e-5
algo    = Normal
lwave   = true
lcharg  = true
lreal   = Auto
lplane  = true
npar    = 4
nsim    = 6
"""
import ConfigParser
from ase.io import read
from ase.optimize.fire import FIRE

class calc_constructor:
    def __init__(self, calc_type, calc_paras):
       self.calc_type = calc_type
       self.calc_paras=calc_paras

    def get_calculator(self):
       calc_paras = self.calc_paras
       if self.calc_type == 'vasp':
          from ase.calculators.vasp import Vasp
          calc  = Vasp(xc     = calc_paras.get('calculator','xc'),       # for MD, coarse prec
                       prec   = calc_paras.get('calculator','prec'),
                       istart = calc_paras.getint('calculator','istart'),
                       icharg = calc_paras.getint('calculator','icharg'),
                       ispin  = calc_paras.getint('calculator','ispin'),
                       encut  = calc_paras.getfloat('calculator','encut'),
                       ismear = calc_paras.getint('calculator','ismear'),
                       sigma  = calc_paras.getfloat('calculator','sigma'),
                       nelm   = calc_paras.getint('calculator','nelm'),
                       nelmin = calc_paras.getint('calculator','nelmin'),
                       ediff  = calc_paras.getfloat('calculator','ediff'),
                       algo   = calc_paras.get('calculator','algo'),
                       lwave  = calc_paras.getboolean('calculator','lwave'),
                       lcharg = calc_paras.getboolean('calculator','lcharg'),
                       lreal  = calc_paras.get('calculator','lreal'),
                       lplane = calc_paras.getboolean('calculator','lplane'),
                       npar   = calc_paras.getint('calculator','npar'),
                       nsim   = calc_paras.getint('calculator','nsim')
                       )
       return calc

def opt_structure(atoms, interval, opt_calculator, optimizer):
    atoms.set_calculator(opt_calculator)
    opt= optimizer(atoms=atoms, maxmove = 1.0, dt = 0.2, dtmax = 1.0, logfile='geo_opt.dat', trajectory='geo_opt.traj')

    opt_traj = []
    e_log = []
    forces = []
    def log_traj(atoms=atoms):
        opt_traj.append(atoms.copy())
        epot=atoms.get_potential_energy(force_consistent=True)
        f = atoms.get_forces()
        e_log.append(epot)
        forces.append(f)
    opt.attach(log_traj, interval=1)
    opt.run(fmax=0.02, steps=300)
    return opt_traj, e_log, forces

def main():
    #config = ConfigParser.SafeConfigParser()
    config = ConfigParser.ConfigParser()
    config.read('config.ini')
    #print config
    main_paras = dict(config.items('main'))
    #print main_paras
    traj = read(main_paras['structurefile'],
                index=main_paras['structure_slice'], 
                format=main_paras['fileformat'])
    print("# of structures:", len(traj))
    calc = calc_constructor(calc_type=config.get('main','calc_type'), 
                            calc_paras=config)

    calc = calc.get_calculator()

    logfile = open('trajectory.xyz', 'w')
    numb_structure = 40
    cellsize = 15
    atoms_index = 0
    for atoms in traj:
        atoms.set_cell([[cellsize,0,0],[0,cellsize,0],[0,0,cellsize]],scale_atoms=False)
        atoms.set_pbc((True, True, True))
        atoms.center()
        opt_traj, e_log, fs = opt_structure(atoms=atoms, interval=1, opt_calculator=calc, optimizer=FIRE)
        if len(opt_traj) < numb_structure:
           interval = 1
        else:
           interval = int(len(opt_traj)/numb_structure)
        for i in range (len(opt_traj)):
            if i % interval == 0 or i == len(opt_traj)-1:
               if i == len(opt_traj)-1:
                  label = 'local minimum'
               else:
                  label = 'nonstable'
               logfile.write("%d\n"%(len(atoms)))
               logfile.write("%s  %s: %15.6f\n"%(label, 'Energy',e_log[i]))
               for atom in opt_traj[i]:
                   j = atom.index
                   logfile.write("%s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n"%(atom.symbol, atom.x,
                                 atom.y, atom.z, fs[i][j][0], fs[i][j][1], fs[i][j][2]))
            logfile.flush()
        print("{:d} iterations were took to optimize the {:d} structure:".format(i,atoms_index))
        atoms_index += 1
    logfile.close()

if __name__ == '__main__':
    main()

