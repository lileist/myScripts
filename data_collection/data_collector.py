import ConfigParser

config = ConfigParser()
config['vasp']={'xc':'PBE',
                'prec':'Medium',
                'istart': 1,
                'icharg': 2,
                'ispin': 1,
                'encut': 300.,
                'ismear': 0,
                'sigma': 0.1,
                'nelm': 200,
                'nelmin': 4,
                'ediff': 1e-5,
                'algo': 'Fast',
                'lwave': 'false',
                'lcharg': 'false',
                'lreal': 'Auto',
                'lplane': 'true',
                'npar': 4,
                'nsim': 6,
               }


class construct_calculator:
    def __init__(self, calc_type, calc_paras):
       self.calc_type = calc_type
       self.calc_paras=calc_paras

    def get_calculator(self):
       calc_paras = self.calc_paras
       if self.calc_type == 'vasp':
          from ase.calculators.vasp import Vasp
          calc  = Vasp(xc     = calc_paras.get('xc')       # for MD, coarse prec
                       prec   = calc_paras.get('prec'),
                       istart = calc_paras.getint('istart'),
                       icharg = calc_paras.getint('icharg'),
                       ispin  = calc_paras.getint('ispin'),
                       encut  = calc_paras.getfloat('encut'),
                       ismear = calc_paras.getint('ismear'),
                       sigma  = calc_paras.getfloat('sigma'),
                       nelm   = calc_paras.getint('nelm'),
                       nelmin = calc_paras.getint('nelmin'),
                       ediff  = calc_paras.getint('ediff'),
                       algo   = calc_paras.get('algo'),
                       lwave  = calc_paras.getboolean('lwave'),
                       lcharg = calc_paras.getboolean('lcharg'),
                       lreal  = calc_paras.get('lreal'),
                       lplane = calc_paras.getboolean('lplane'),
                       npar   = calc_paras.getint('npar'),
                       nsim   = calc_paras.getint('nsim')
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
        epot=atoms.get_potential_energy()
        f = atoms.get_forces()
        e_log.append(epot)
        forces.append(f)
    opt.attach(log_traj, interval=1)
    opt.run(fmax=0.02, steps=300)
    return opt_traj, e_log, forces

def main():
    arg = sys.argv
    paras = readinputs(arg[1])

    config = ConfigParser()
    config.read('config.ini')
  
    traj = read(config.get('main','structurefile'),index=":", format=config.get('main','fileformat'))

    calc = calc_constructor(calc_type=config.get('main','calc_type'), 
                            calc_paras=config['calculator'])

    calc = calc.get_calculator()

    logfile = open('trajectory.xyz', 'w')
    numb_structure = 40

    for atoms in traj:
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
    logfile.close()

if __name__ == '__main__':
    main()

