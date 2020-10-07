import subprocess
import numpy as np

compiled_fname 	= "ehrenfest_parallel.py"
job_script     	= "sub.sh"
input_0 		= "input0.txt"	#"input0.txt"
input_i			= "input.txt"

p = np.loadtxt(input_0, dtype=int)
ncore = int(p[1])
#Input text files
for seed in range(1,ncore+1):
	np.savetxt(input_i, p.astype(int), fmt='%i')
	file = open(input_i, 'a')
	file.write(str(seed) +"\n")
	file.close()

	file = open("bashrun.sh", "w")
	file.write("#!/bin/bash\n")
	file.write("mkdir "+ str(seed) +"\n")
	
	file.write("cp "+ compiled_fname+ " " + str(seed) + "\n")
	file.write("cp "+ job_script	+ " " + str(seed) + "\n")
	file.write("mv "+ input_i		+ " " + str(seed) + "/"+ input_0+"\n")

	file.write("cd "+ str(seed) + "\n")
	#file.write("python "+ compiled_fname + "\n")		
	file.write("qsub "+ job_script + "\n")
	file.write("cd .. \n") 
	file.close()

	subprocess.call("(chmod +x bashrun.sh)", shell = True)
	subprocess.call("(./bashrun.sh)", shell = True)






