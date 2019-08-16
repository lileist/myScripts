import os,sys
from ase.io import read, Trajectory
import pandas as pd

args = sys.argv
start = int(args[1])
end  = int(args[2])
current = os.getcwd()
traj = Trajectory('akmc.traj','w')
for i in range(start, end+1):
   os.chdir(current+'/'+str(i))
   process = pd.read_table('processtable', delimiter = r'\s+', skiprows = [0])
   process.columns = ['processID', 'saddle_e', 'prefactor', 'productID', 'product_e', 'product_prefactor', 'barrier', 'rate', 'repeats']
   for j in process['processID']:
      for prename in ['reactant_','saddle_','product_']:
         traj.write(read('procdata/'+prename+str(j)+'.con',format='eon'))

  

