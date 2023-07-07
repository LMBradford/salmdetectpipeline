import glob
import os

# Find input files
FILES = glob.glob("*_R1.fq.gz") #Files must be in working directory

#Using os.system,
#open the zipped file
#Take every 4th line, starting from the 1st (name lines in fastq)
#Remove the first character (@ in fastq)
#Cut everything after the space 
#Print to new file, retaining sample name

for file in FILES:
	name = file.split("_")[0]
	cmd = f"zcat {file} | sed -n '1~4p' | cut -c2- | cut -d ' ' -f 1 > {name}_names.txt"
	print(cmd)
	os.system(cmd)
