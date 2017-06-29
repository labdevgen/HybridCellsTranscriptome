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

work_dir = "/mnt/storage/home/vsfishman/tmp/m5s_2017/"
transcr_index="/mnt/storage/home/vsfishman/tmp/m5s_2017/transcriptome_data/mm10known"
genome_index="/mnt/storage/home/vsfishman/HiC/bin/bowtie2/index/mm10"
fastq_dir = work_dir+"fastq/"

#FILLING FILE NAMES
allfiles=os.listdir(fastq_dir)
files=dict([(fastq_dir+i,"s"+i.split("_")[0]+"_2017") for i in allfiles if i.endswith(".fastq.gz")])

for i in files.keys():
	f=open("tophut.q","w")
	f.write("""#!/bin/bash
#PBS -V -r n 
#PBS -l select=1:ncpus=4:mem=100gb,cput=64:10:00,walltime=32:15:00
#PBS -q vkop2q
#PBS -N thut"""+files[i]+"""
#PBS -j oe
date"""+
"\nexport PATH=$HOME/Distr/bowtie2-2.2.1:$PATH\n" #change bowtie2 version to old 2.2.1
)
	
	f.write("cd "+work_dir+"\n")
	
	# if (not os.path.isfile(i+".fastq.trimmed")) and (not os.path.isfile(i+".fastq.trimmed.gz")): #if we do not have any trimmed file
		# f.write("java -jar /mnt/storage/home/vsfishman/Distr/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 4 "+i+".fastq "+i+".fastq.trimmed HEADCROP:3\n") #create it
		# f.write("gzip "+i+".fastq "+"\n") #and gztip initial fastq
	
	
	f.write("tophat --prefilter-multihits -p 12 -o "+work_dir+"tophat_out/"+files[i]+" --transcriptome-index "+transcr_index+" "+genome_index+" "+i+"\n") #run tophat with trimmed or trimmed.gz file
	
#	f.write("python2.7 infect_align_pareser.py "+work_dir+"tophat_out/"+files[i]+"_trimmed/unmapped.bam\n")
	f.write("cp "+work_dir+"tophat_out/"+files[i]+"/align_summary.txt "+work_dir+files[i]+"_align_summary.txt\n")
#	f.write("cp "+work_dir+"tophat_out/"+files[i]+"/unmapped.bam.fasta.infect_bwtout.log "+work_dir+files[i]+"_unmapped.bam.fasta.infect_bwtout.log")
	f.close()
#	executeBashCommand("qsub tophut.q")