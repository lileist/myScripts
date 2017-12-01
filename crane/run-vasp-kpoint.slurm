#!/bin/sh
#SBATCH -J znte
#SBATCH --nodes=1 --ntasks-per-node=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=168:00:00
#SBATCH --partition=batch
#SBATCH --error=TestJob.stderr
#SBATCH --output=TestJob.stdout

module load compiler/intel/12
module load intel-mkl/12
module load openmpi/1.6

#ulimit -s unlimited

nodelist >  ./machines.LINUX
echo $SLURM_JOB_ID > ./jobid
#mpiexec /home/zeng/programs/vasp5.3.3/vasp.5.3/vasp-intel12 > out
#mpiexec /home/zeng/programs/vasp5.2.2-NEB/vasp.5.2/vasp-intel11-mkl10 > out
# To run vasp5.3.3 across multiple nodes
#mpiexec /bin/sh -c 'ulimit -s unlimited; /home/zeng/programs/vasp5.3.3/vasp.5.3/vasp-intel12 > out'

for i in 6 8 10 12 14 16 18 20 22 24 26 28 30
do
rm WAVECAR
#i=$(echo "$i"|bc -l) 
#a=$(echo "$i"|bc -l) 
## echo "$a $b"


cat >KPOINTS << CATY
K-Points
0
Gamma
$i $i 1
0 0 0
CATY
echo '==========poscar============'  >> out
cat POSCAR  >> out
echo '============================='  >> out

mpiexec /home/zeng/programs/vasp5.3.3/vasp.5.3/vasp-intel12 > out

echo -n "0 0 $i :  " >> volum_energy.txt
grep "energy without entropy" OUTCAR |tail -1|awk '{print $8}' >> volum_energy.txt
#grep "energy without entropy" OUTCAR |tail -1
echo $i >> contcar-save
cat CONTCAR > $i
echo ' '
grep "energy without entropy" OUTCAR |tail -1   >> out
grep "EENTRO =" OUTCAR | tail -1   >> out
done

