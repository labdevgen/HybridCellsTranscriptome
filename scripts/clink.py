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
tophat_out_dir="/mnt/storage/home/vsfishman/tmp/m5s/tophat_out/"
gff_file="/mnt/storage/home/vsfishman/tmp/m5s/transcriptome_data/mm10known.gff"
fasta_file="/mnt/storage/home/vsfishman/tmp/m5s/Bowtie2Index/genome.fa"

file_names=range(1,2)
for i,v in enumerate(file_names):
		file_names[i]="%(v)03d"%{"v":v}


for i in file_names:
	f=open("clink.q","w")
	f.write("#!/bin/bash \n#PBS -V -r n \n#PBS -l select=1:ncpus=4:mem=60gb,cput=64:10:00,walltime=32:15:00\n#PBS -q vkop2q\n#PBS -N clink"+i+"\n#PBS -j oe\ndate\n#!/bin/sh\ndate\n")
	f.write("cd "+work_dir+"\n")
	f.write("cufflinks -b "+fasta_file+" -G "+gff_file+" --multi-read-correct -p 4 -o clinks_out/"+i+" --no-update-check "+tophat_out_dir+i+"_trimmed/accepted_hits.bam\n")
	f.close()
	executeBashCommand("qsub clink.q")