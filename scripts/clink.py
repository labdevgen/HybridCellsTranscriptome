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
tophat_out_dir="/mnt/storage/home/vsfishman/tmp/m5s_2017/tophat_out/"
gff_file="/mnt/storage/home/vsfishman/tmp/m5s_2017/transcriptome_data/mm10known.gff"
fasta_file="/mnt/storage/home/vsfishman/HiC/bin/bowtie2/index/mm10.fa"

allfiles=os.listdir(tophat_out_dir)

for i in allfiles:
	f=open("clink.q","w")
	f.write("""#!/bin/bash
#PBS -V -r n 
#PBS -l select=1:ncpus=12:mem=24gb,cput=164:10:00,walltime=32:15:00
#PBS -q bl2x220g7q
#PBS -N clink"""+i+"""
#PBS -j oe
date"""+
"\nexport PATH=$HOME/Distr/bowtie2-2.2.1:$PATH\n"+ #change bowtie2 version to old 2.2.1
"export PATH=$HOME/Distr/samtools-0.1.19/:$PATH\n" #change samtools version to old 0.1.19
)
	f.write("cd "+work_dir+"\n")
	f.write("cufflinks -b "+fasta_file+" -G "+gff_file+" --multi-read-correct -p 12 -o clinks_out/"+i+" --no-update-check "+tophat_out_dir+i+"/accepted_hits.bam\n")
	f.close()
#	executeBashCommand("qsub clink.q")