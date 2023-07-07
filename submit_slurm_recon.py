import os
import itertools
import time

reps = range(6,21)
spikes = [5,10,50,100,10000,100000]

for rep in reps:
  for spike in spikes:
    name = f'salm{spike}rep{rep}'
    command = f"sbatch --export=ALL,N={name} SLURM_reconstitutelibs.sh"
    print(command)
    os.system(command)
    time.sleep(5)
