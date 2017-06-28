import os

def executeBashCommand(bashCommand,doPrint=True):
  import subprocess
  print "executing command %s \n" % (bashCommand)
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  if doPrint:
    print p.communicate()[0]
    return ""
  else:
    return p.communicate()[0]

work_dir = "/mnt/storage/home/vsfishman/tmp/m5s/"
transcr_index="/mnt/storage/home/vsfishman/tmp/m5s/transcriptome_data/mm10known"
genome_index="/mnt/storage/home/vsfishman/tmp/m5s/Bowtie2Index/genome"
fastq_dir = work_dir+"fastq/"

file_names=range(2,19)
for i,v in enumerate(file_names):
		file_names[i]="%(v)03d"%{"v":v}


#FILLING FILE NAMES
allfiles=os.listdir(fastq_dir)
files={}
for i in allfiles:
	for j in file_names:
		if (i.find(j) == 6):
			files[fastq_dir+i.split(".")[0]] = j
			break

for i in files.keys():
	f=open("tophut.q","w")
	f.write("#!/bin/bash \n#PBS -V -r n \n#PBS -l select=1:ncpus=4:mem=100gb,cput=64:10:00,walltime=32:15:00\n#PBS -q vkop2q\n#PBS -N thut"+files[i]+"\n#PBS -j oe\ndate\n#!/bin/sh\ndate\n")
	f.write("cd "+work_dir+"\n")
	if (not os.path.isfile(i+".fastq.trimmed")) and (not os.path.isfile(i+".fastq.trimmed.gz")): #if we do not have any trimmed file
		f.write("java -jar /mnt/storage/home/vsfishman/Distr/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 "+i+".fastq "+i+".fastq.trimmed HEADCROP:3\n") #create it
		f.write("gzip "+i+".fastq "+"\n") #and gztip initial fastq
	
	
	f.write("tophat --prefilter-multihits -p 3 -o "+work_dir+"tophat_out/"+files[i]+"_trimmed --transcriptome-index "+transcr_index+" "+genome_index+" "+i+".fastq.trimmed*\n") #run tophat with trimmed or trimmed.gz file
	
	if (not os.path.isfile(i+".fastq.trimmed.gz")): #gzip trimmed file if not done before
		f.write("gzip "+i+".fastq.trimmed"+"\n")

	f.write("python2.7 infect_align_pareser.py "+work_dir+"tophat_out/"+files[i]+"_trimmed/unmapped.bam\n")
	f.write("cp "+work_dir+"tophat_out/"+files[i]+"_trimmed/align_summary.txt "+work_dir+files[i]+"_align_summary.txt")
	f.write("cp "+work_dir+"tophat_out/"+files[i]+"_trimmed/unmapped.bam.fasta.infect_bwtout.log "+work_dir+files[i]+"_unmapped.bam.fasta.infect_bwtout.log")
	f.close()
	executeBashCommand("qsub tophut.q")