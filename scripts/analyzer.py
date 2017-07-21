import numpy as np
import numpy.ma as ma
import math
import matplotlib.pyplot as plt


################################functions###############################
def del_intersection_from_list(l1,l2):
  return np.delete(l1,np.nonzero(np.in1d(l1,l2)))

def getMask(exp,sample1,sample2,significance,fold_change=-1, FPKM_threshold=0):
#return a boolean array where all results of comparison of sample1 and sample2 are NOT "significance"
#and fold change is < fold_change
#and FPKM_threshold for either sample1 and sample2 <= FPKM_threshold
#masking with this array (i.e. eleminating elements of original array where masking array is true) will leave only valid samples
 
#  MAX_UNSIGNIFICANT_DIFFERENCE=math.log(100,2)
  r1=np.logical_or(np.logical_and(exp[:,4]==sample1,exp[:,5]==sample2),np.logical_and(exp[:,4]==sample2,exp[:,5]==sample1))
  r1=np.logical_and(r1,exp[:,13]==significance)
#  if significance == 'no':
#    r1=np.logical_and(r1,np.absolute(exp[:,9].astype(np.float))<=MAX_UNSIGNIFICANT_DIFFERENCE)
#  if significance == 'yes':
#    r1=np.logical_or(r1,np.absolute(exp[:,9].astype(np.float))>MAX_UNSIGNIFICANT_DIFFERENCE)
  if fold_change != -1:
    fold_change=math.log(fold_change,2)
    r1=np.logical_and(r1,np.absolute(exp[:,9].astype(np.float))>=fold_change)
  if FPKM_threshold > 0:
    frkpm_t_array=np.logical_or(np.absolute(exp[:,7].astype(np.float))>=FPKM_threshold,np.absolute(exp[:,8].astype(np.float))>=FPKM_threshold)
    r1=np.logical_and(r1,frkpm_t_array)
  return np.logical_not(r1)

#Column number	Column name	Example	Description
#1	Tested id	XLOC_000001	A unique identifier describing the transcipt, gene, primary transcript, or CDS being tested
#2	gene	Lypla1	The gene_name(s) or gene_id(s) being tested
#2	gene_id	Lypla1	The gene_name(s) or gene_id(s) being tested
#4	locus	chr1:4797771-4835363	Genomic coordinates for easy browsing to the genes or transcripts being tested.
#5	sample 1	Liver	Label (or number if no labels provided) of the first sample being tested
#6	sample 2	Brain	Label (or number if no labels provided) of the second sample being tested
#7	Test status	NOTEST	Can be one of OK (test successful), NOTEST (not enough alignments for testing), LOWDATA (too complex or shallowly sequenced), HIDATA (too many fragments in locus), or FAIL, when an ill-conditioned covariance matrix or other numerical exception prevents testing.
#8	FPKMx	8.01089	FPKM of the gene in sample x
#9	FPKMy	8.551545	FPKM of the gene in sample y
#10	log2(FPKMy/FPKMx)	0.06531	The (base 2) log of the fold change y/x	
#11	test stat	0.860902	The value of the test statistic used to compute significance of the observed change in FPKM
#12	p value	0.389292	The uncorrected p-value of the test statistic
#13	q value	0.985216	The FDR-adjusted p-value of the test statistic
#14	significant	no	Can be either "yes" or "no", depending on whether p is greater then the FDR after Benjamini-Hochberg correction for multiple-testing

################################end functions###############################

################################initialization of variables###############################
print "loading data"

FPKM_threshold_N=10

fibs=["tmf2","tmf5"]
esc=["tme14","tme17"]
tauGFP="tau-GFP4N"
m5s="m5S4N"

#fpkm_file_filds_N = {"m5S":1,""tmf2":2,"tmf5":3,"tme14":4,"tauGFP":5,"tme17":6,m5S-4N":7,"tau-GFP4N":8} #numbers of cols (0-basez) in fpkm_tracking file for each cell line, where it's fpkm is writen.
fpkm_file_filds_N = {"tmf2":2,"tmf5":3,"tme14":4,"tme17":6,"m5S4N":7,"tau-GFP4N":8} #numbers of cols (0-basez) in fpkm_tracking file for each cell line, where it's fpkm is writen.
for i in fpkm_file_filds_N.keys(): #numbers of cols (0-basez) in fpkm_tracking file for each cell line, where it's fpkm is writen.
    fpkm_file_filds_N[i] = (fpkm_file_filds_N[i]-1)*4+9
    
#loading "gene_exp.diff" and "genes.fpkm_tracking" files
dt = [('test_id', 'S18'), ('gene_id', 'S55'), ('gene', 'S55'), ('locus', 'S34'), ('sample_1', 'S10'), ('sample_2', 'S10'), ('status', 'S6'), ('value_1', '<f8'), ('value_2', '<f8'), ('log2fold_change', '<f8'), ('test_stat', '<f8'), ('p_value', '<f8'), ('q_value', '<f8'), ('significant', 'S3')]
exp = np.transpose(np.loadtxt("gene_exp.diff", dtype=dt,skiprows=1,unpack=True))

#dt = [('f0', 'S18'), ('f1', 'S1'), ('f2', 'S1'), ('f3', 'S18'), ('f4', 'S18'), ('f5', 'S100'), ('f6', 'S34'), ('f7', 'S1'), ('f8', 'S1'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'), ('f12', 'S2'), ('f13', '<f8'), ('f14', '<f8'), ('f15', '<f8'), ('f16', 'S6'), ('f17', '<f8'), ('f18', '<f8'), ('f19', '<f8'), ('f20', 'S6'), ('f21', '<f8'), ('f22', '<f8'), ('f23', '<f8'), ('f24', 'S2'), ('f25', '<f8'), ('f26', '<f8'), ('f27', '<f8'), ('f28', 'S2'), ('f29', '<f8'), ('f30', '<f8'), ('f31', '<f8'), ('f32', 'S2')]
dt = [('f0', 'S18'), ('f1', 'S1'), ('f2', 'S1'), ('f3', 'S18'), ('f4', 'S18'), ('f5', 'S104'), ('f6', 'S34'), ('f7', 'S1'), ('f8', 'S1'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'), ('f12', 'S2'), ('f13', '<f8'), ('f14', '<f8'), ('f15', '<f8'), ('f16', 'S6'), ('f17', '<f8'), ('f18', '<f8'), ('f19', '<f8'), ('f20', 'S6'), ('f21', '<f8'), ('f22', '<f8'), ('f23', '<f8'), ('f24', 'S2'), ('f25', '<f8'), ('f26', '<f8'), ('f27', '<f8'), ('f28', 'S2'), ('f29', '<f8'), ('f30', '<f8'), ('f31', '<f8'), ('f32', 'S2'), ('f33', '<f8'), ('f34', '<f8'), ('f35', '<f8'), ('f36', 'S6'), ('f37', '<f8'), ('f38', '<f8'), ('f39', '<f8'), ('f40', 'S6')]
fpkms = np.transpose(np.loadtxt("genes.fpkm_tracking", dtype=dt,skiprows=1,unpack=True))

global genes_dict #in this dictionary stored arrays of genes sorted by gene category and clone name
genes_dict={}
genes_dict["proper"]={}
genes_dict["not_reprogrammed"]={}
#genes_dict["intermediate_nsign"]={}
#genes_dict["intermediate_sign"] = {}
genes_dict["intermediate"] = {} #just a sum of both intermediate (Significant and NotSignigicant)
genes_dict["new_exp_profile"] = {}


################################end initialization of variables###############################

################################Initial analyzses of genes: number of genes, of statistically signif tests, of low expressed genes and etc###############################
print "Counting"
 
all_genes = np.unique(exp[:,0])
all_genes_N = len(all_genes)

too_low_expression_list = []
for i in fpkms:
  fpkm_values = []
  for j in fpkm_file_filds_N.keys():
    fpkm_values.append(float(i[fpkm_file_filds_N[j]]))
  if np.max(fpkm_values) < FPKM_threshold_N:
    too_low_expression_list.append(i[0])
 

filtered_genes=np.array(too_low_expression_list)

genes_with_notOK_test_values = np.unique(exp[np.where(exp[:,6]!="OK")][:,0])

filtered_genes=np.concatenate((genes_with_notOK_test_values,filtered_genes))

filtered_genes = np.unique(filtered_genes)

print "genes where statistical tests failed or expression is too low ",len(filtered_genes)


same_in_tauGFP_and_m5s = ma.masked_array(exp[:,0],getMask(exp,tauGFP,m5s,'no')).compressed()
same_in_tauGFP_and_m5s = del_intersection_from_list(same_in_tauGFP_and_m5s,filtered_genes) #filter statisticaly unsignificant genes from this list
print "Genes with same expression between m5s/tauGFP ",len(same_in_tauGFP_and_m5s)

different_tauGFP_m5s=ma.masked_array(exp[:,0],getMask(exp,tauGFP,m5s,'yes',FPKM_threshold=FPKM_threshold_N,fold_change=2)).compressed()

different_tauGFP_m5s_rejected=ma.masked_array(exp[:,0],getMask(exp,tauGFP,m5s,'yes')).compressed()
different_tauGFP_m5s_rejected=del_intersection_from_list(different_tauGFP_m5s_rejected,different_tauGFP_m5s)

different_tauGFP_m5s = del_intersection_from_list(different_tauGFP_m5s,filtered_genes) #filter statisticaly unsignificant genes from this list
different_tauGFP_m5s_rejected = del_intersection_from_list(different_tauGFP_m5s_rejected,filtered_genes) #filter statisticaly unsignificant genes from this list

print "Genes with too small difference between m5s and tauGFP ",len(different_tauGFP_m5s_rejected)
print "Genes different between m5s and tauGFP ",len(different_tauGFP_m5s)

print "out of ",all_genes_N

################################end Initial analyzses of genes: number of genes, of stat sign tests, of low expressed genes and etc###############################

def compare_with_parent(sample,same_parent,different_parent,genes_to_count):
  global genes_dict
  print "Sample ",sample
  
  same_as_in_same_parent = ma.masked_array(exp[:,0],getMask(exp,i,same_parent,'no')).compressed()
  same_as_in_different_parent = ma.masked_array(exp[:,0],getMask(exp,i,different_parent,'no')).compressed()
  
  #step one. Genes that are properly reprogrammed
  proper_genes = np.intersect1d(same_as_in_same_parent,genes_to_count) #condition 1 - keep only genes where parents have diferent
   
  intersec_idxs=np.nonzero(np.in1d(proper_genes,same_as_in_different_parent)) #condition 2 - remove genes where expression is not significantly different between hybrid and both parents
  intermediate_genes = np.array(proper_genes[intersec_idxs],copy=True)
  proper_genes = np.delete(proper_genes,intersec_idxs)
  print "Properly reprogrammed ",len(proper_genes)
  genes_dict["proper"][sample]=proper_genes
  
  #now let's calculate not reprogrammed genes
  not_reprogrammed = np.intersect1d(same_as_in_different_parent,genes_to_count) #condition 1 - remove genes where parents have same expression
  
  intersec_idxs=np.nonzero(np.in1d(not_reprogrammed,intermediate_genes)) #condition 2 - remove genes where expression is not significantly different between hybrid and both parents
  not_reprogrammed = np.delete(not_reprogrammed,intersec_idxs)

  print "Completetly not reprogrammed ",len(not_reprogrammed)
  genes_dict["not_reprogrammed"][sample]=not_reprogrammed
#  genes_dict["intermediate_nsign"][sample]=intermediate_genes
  genes_dict["intermediate"][sample]=intermediate_genes


################################stage 1##############################
#stage 1. How complete is reprogramming. Select all genes, that change there level properly (and were different in tauGFP and m5s)

print "Stage1. How complete is reprogramming."
for i in fibs:
  compare_with_parent(i,m5s,tauGFP,different_tauGFP_m5s)
for i in esc:
  compare_with_parent(i,tauGFP,m5s,different_tauGFP_m5s)

################################end stage 1##############################

################################stage 2. Incomplete reprogramming##############################
#stage 2. Genes, that have intermediate level of reprogramming: level is changed and it is between tauGFP and m5s 

print "Stage2. Incomplete reprogramming"
for i in fibs+esc:
  print i
  diff_from_m5s = ma.masked_array(exp[:,0],getMask(exp,i,m5s,'yes',FPKM_threshold=FPKM_threshold_N,fold_change=-1)).compressed()
  diff_from_m5s_rejected = ma.masked_array(exp[:,0],getMask(exp,i,m5s,'yes')).compressed()
  
  diff_from_tauGFP = ma.masked_array(exp[:,0],getMask(exp,i,tauGFP,'yes',FPKM_threshold=FPKM_threshold_N,fold_change=-1)).compressed()
  diff_from_tauGFP_rejected = ma.masked_array(exp[:,0],getMask(exp,i,tauGFP,'yes')).compressed()
  
  differ_from_both_parents = np.intersect1d (diff_from_m5s,diff_from_tauGFP)
  differ_from_both_parents = np.intersect1d(differ_from_both_parents,different_tauGFP_m5s)
  
  differ_from_both_parents_rejected = np.intersect1d (diff_from_m5s_rejected,diff_from_tauGFP_rejected)
  differ_from_both_parents_rejected = np.intersect1d (differ_from_both_parents_rejected,different_tauGFP_m5s)
  differ_from_both_parents_rejected = del_intersection_from_list(differ_from_both_parents_rejected,differ_from_both_parents)

  differ_from_both_parents_between=[]
  differ_from_both_parents_not_between=[]
  for j in differ_from_both_parents:
    fpkms_idx = np.where(fpkms[:,0]==j)
    m5s_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[m5s]][0][0])
    tauGFP_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[tauGFP]][0][0])
    sample_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[i]][0][0])
    if m5s_fpkm>sample_fpkm>tauGFP_fpkm or m5s_fpkm<sample_fpkm<tauGFP_fpkm:
      differ_from_both_parents_between.append(j)
    else:
      
      if abs(sample_fpkm-m5s_fpkm) > abs(sample_fpkm-tauGFP_fpkm): 
	tmp=tauGFP_fpkm 
      else: 
	tmp=m5s_fpkm
      
      if (sample_fpkm== tmp == 0):
	differ_from_both_parents_rejected=np.append(differ_from_both_parents_rejected,j)
      elif (sample_fpkm, tmp) == 0:
	differ_from_both_parents_not_between.append(j)
      elif max(tmp,sample_fpkm)/min(tmp,sample_fpkm)> math.pow(2,0.5):
	differ_from_both_parents_not_between.append(j)
      else:
	differ_from_both_parents_rejected=np.append(differ_from_both_parents_rejected,j)
  
  print "Rejected because both clones have too low expression ",len(differ_from_both_parents_rejected)      
  print "Total number of genes expressed different from both parents ",len(differ_from_both_parents)
  print "From this, between them: ",len(differ_from_both_parents_between)," not between ",len(differ_from_both_parents_not_between)
#  genes_dict["intermediate_sign"][i]=differ_from_both_parents_between
  genes_dict["intermediate"][i]=np.append(genes_dict["intermediate"][i],differ_from_both_parents_between)
  genes_dict["new_exp_profile"][i]=differ_from_both_parents_not_between

################################end stage 2. Incomplete reprogramming##############################

################################Image: Bar-plot for ration between  "proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile" genes#########################
 
print "Preparing image"
#Bar-plot for ration between  "proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile" genes

plt_index = np.arange(len(esc+fibs))
plt_width = 0.35
plt_colors=[]
#dict keys are: "proper","not_reprogrammed","intermediate","intermediate_nsign","intermediate_sign","new_exp_profile"

colors = {"proper":"b","not_reprogrammed":"g","intermediate_nsign":"c","intermediate_sign":"r","intermediate":"r","new_exp_profile":"y"}

genes_dict_percantage = {}
for i in genes_dict.keys():
  genes_dict_percantage[i]={}
  for j in genes_dict[i].keys():
    genes_dict_percantage[i][j]=0.0

for i in esc+fibs:
    total_number_of_genes = 0.0
    for j in sorted(genes_dict.keys()):
      total_number_of_genes += len(genes_dict[j][i])
    for j in sorted(genes_dict.keys()):
      print i,j,len(genes_dict[j][i]),(len(genes_dict[j][i])/total_number_of_genes)*100
      genes_dict_percantage[j][i] = (len(genes_dict[j][i])/total_number_of_genes)*100

fig = plt.figure()
ax = plt.subplot(111)

plots = []
b=np.zeros(len(fibs+esc))
#for j in ["proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile"]:
for j in ["proper","not_reprogrammed","intermediate","new_exp_profile"]:
  values = [genes_dict_percantage[j][i] for i in sorted(fibs+esc)]
  p = ax.bar(plt_index, values, plt_width, color=colors[j],bottom=b)
  b += values
  plots.append(p)
  p.legend = j
 
plt.ylabel('% of genes')
plt.xticks(plt_index+plt_width/2., sorted(fibs+esc))
plt.yticks(np.arange(0,101,10))

# Shink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])

# Put a legend below current axis
#ax.legend([p[0] for p in plots], ["Reprogrammed","Not Reprogrammed","Intermediate","Novel expression pattern"],loc='upper center', \
ax.legend([p[0] for p in plots], ["Subgroup 1","Subgroup 2","Subgroup 3","Subgroup 4"],loc='upper center', \
	bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=2)
 
plt.savefig("Genes_classification.pdf",dpi=600)
plt.clf()

#plot a picture with ratio of analyzed for each clone genes from total number of genes in transcriptome

################################End Image: Bar-plot for ration between  "proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile" genes#########################

################################Image:Genes statistics (total number, low expressed, etc)#########################

fig = plt.figure()
ax = plt.subplot(111)

analyzed_genes = 0
for j in fibs+esc:
  for i in genes_dict.keys():
    analyzed_genes += float(len(genes_dict[i][j]))

analyzed_genes /= len(fibs+esc)
additionaly_rejected_genes = len(different_tauGFP_m5s)-analyzed_genes
				
values = [all_genes_N,round(len(filtered_genes)+additionaly_rejected_genes),len(same_in_tauGFP_and_m5s),len(different_tauGFP_m5s_rejected),round(len(different_tauGFP_m5s)-additionaly_rejected_genes)]

p=ax.bar(np.arange(len(values)),values,plt_width)

plt.ylabel('Number of genes')

labels = ["Total\ngenes","Not enough\ndata","Same expr.\nin parental lines","Minor expr. difference\nin parental lines","Differently expr.\nin parental lines"]
plt.xticks(np.arange(len(labels))+plt_width/2., labels)
ax.tick_params(axis='x', labelsize=8)

plt.savefig("Genes_analyzed.pdf",dpi=600)
plt.clf()

################################end image:Genes statistics (total number, low expressed, etc)#########################

################################Stage 3.Analizys of fold change in different categories of genes, plot a picture for fold change#########################
print "Stage 3.Analizys of fold change in different categories of genes"

#during this step, also best examples of genes  category "not_reprogrammed" will be found
#best means with maximum difference between parents

number_of_examples_to_search=10
print "Best examples of genes  the category not_reprogrammed:"
category="not_reprogrammed"
examples={}
for j in  [esc,fibs]:
    shared_genes=np.intersect1d(genes_dict["not_reprogrammed"][j[0]], genes_dict["not_reprogrammed"][j[1]])
    j=j[0] #converting  j from list to as scalar with a name equal to first element of the list
    examples[j]={}    
    for k in shared_genes:
      fpkms_idx = np.where(fpkms[:,0]==k)
      m5s_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[m5s]][0][0])
      tauGFP_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[tauGFP]][0][0])
      if m5s_fpkm  != 0 and tauGFP_fpkm != 0:
	fold_change = max(tauGFP_fpkm,m5s_fpkm)/min(tauGFP_fpkm,m5s_fpkm)
	fold_change = round (fold_change)
	if fold_change in examples[j].keys():
	  examples[j][fold_change].append(k)
	else:
	  examples[j][fold_change]=[k]
    s=""
    t=number_of_examples_to_search
    for i in sorted(examples[j].keys(),reverse=True):
      if t==0: break
      for k in examples[j][i]:
	s+='"'+str(k)+'" '
	t-=1
	if t ==0: break
    print j,s
  
  
genes_fold_change = {}
for i in genes_dict.keys():
  genes_fold_change[i]={}
  for j in genes_dict[i].keys():
    genes_fold_change[i][j]={}


print "Outputing genes from cat.4, subcat. 'novel' with difference > 10 times:"
for i in genes_dict.keys():
  for j in genes_dict[i]:
    for k in genes_dict[i][j]:

    #determine frkpm for the gene and parent
      fpkms_idx = np.where(fpkms[:,0]==k)
      m5s_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[m5s]][0][0])
      tauGFP_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[tauGFP]][0][0])
      sample_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[j]][0][0])

      if (i == "new_exp_profile" or i=="intermediate"): 
	if abs(sample_fpkm-m5s_fpkm) < abs(sample_fpkm-tauGFP_fpkm): #this will calculate the fold cahange to closest parent
	  parent_fpkm=m5s_fpkm
	else:
	  parent_fpkm=tauGFP_fpkm
      else:
	if j in fibs: #this will calculate the fold cahange to dominant parent
	  parent_fpkm = m5s_fpkm
	else:
	  parent_fpkm = tauGFP_fpkm
    
      #count
      if parent_fpkm != 0 and sample_fpkm != 0:
	fold_change = max(parent_fpkm,sample_fpkm)/min(parent_fpkm,sample_fpkm)
	fold_change = round (fold_change)
	if (fold_change > 10):
	  fold_change = 10.0
	  if i == "new_exp_profile": 
	    print i,j,k,str(fpkms[fpkms_idx,6][0][0]).split(":")[0]
	    
	  if (i == "proper" and fold_change >2): 
	    print "!!!!",i,j,k,str(fpkms[fpkms_idx,6][0][0]).split(":")[0]

	if fold_change in genes_fold_change[i][j].keys():
	  genes_fold_change[i][j][fold_change] += 1.0
	else:
	  genes_fold_change[i][j][fold_change] = 1.0

#transform frequences into %
for i in genes_dict.keys():
  for j in genes_dict[i].keys():
    total_genes = sum (genes_fold_change[i][j].values())
    for k in genes_fold_change[i][j].keys():
      genes_fold_change[i][j][k] = (genes_fold_change[i][j][k]/total_genes)*100

print "Saving picture"
#plot figures      
for i in fibs+esc:
  fig = plt.figure()
  ax = plt.subplot(111)
  plots=[]
#  for j in ["proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile"]:
  for j in ["proper","not_reprogrammed","intermediate","new_exp_profile"]:
    x_values = sorted(genes_fold_change[j][i].keys())
    y_values = [genes_fold_change[j][i][v] for v in sorted(genes_fold_change[j][i].keys())]
    plots.append(ax.plot(x_values,y_values,color = colors[j],ls='-',marker='o'))
  
  plt.ylabel('% of genes')
  plt.xticks(np.arange(0,11))
  plt.yticks(np.arange(0,101,5))
  labels = [str(t)+"-"+str(t+1) for t in xrange(0,11)]
  labels[0] = ""
  labels[1] = "<2"
  labels[-1] = ">10"
  ax.set_xticklabels(labels)
 
  # Shink current axis's height by 10% on the bottom
  box = ax.get_position()
  ax.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width, box.height * 0.7])
  ax.set_xlabel('Fold change between expression in parent and hybrid')
#  ax.xaxis.set_label_coords(x=0,y=0)  
  # Put a legend below current axis

  ax.legend([p[0] for p in plots], ["Subgroup 1","Subgroup 2","Subgroup 3","Subgroup 4"],loc='upper center', \
#  ax.legend([p[0] for p in plots], ["Reprogrammed","Not Reprogrammed","Intermediate state(NS)","Intermediate state(S)","Novel expression pattern"],loc='upper center', \    
	bbox_to_anchor=(0.5, -1*box.height*0.2),fancybox=True, shadow=True, ncol=2)
  plt.savefig(i+"_fold_change.pdf",dpi=600)
  plt.clf()
  
print "Stage 4. Overlap of not reprogrammed genes:"
print "For ESC-like clones: "
overlap=np.intersect1d(genes_dict["not_reprogrammed"][esc[0]], genes_dict["not_reprogrammed"][esc[1]])
print len(overlap)
np.savetxt("ESC_notreprogrammed_genes.txt",overlap,fmt="%s")

print "For Fibroblasts-like clones: "
overlap=np.intersect1d(genes_dict["not_reprogrammed"][fibs[0]], genes_dict["not_reprogrammed"][fibs[1]])
print len(overlap)
np.savetxt("Fib_notreprogrammed_genes.txt",overlap,fmt="%s")

print "For all clones: "
overlap1=np.intersect1d(genes_dict["not_reprogrammed"][fibs[0]], genes_dict["not_reprogrammed"][fibs[1]])
overlap2=np.intersect1d(genes_dict["not_reprogrammed"][esc[0]], genes_dict["not_reprogrammed"][esc[1]])
overlap3= np.intersect1d(overlap1,overlap2)
print len(overlap3)
np.savetxt("All_notreprogrammed_genes.txt",overlap3,fmt="%s")

print "Stage 5. How genes from different categories are distributed on chromosomes"
chromosomes = {} #for each clone it will be a dictionary with chromosome_Number-->Number_of_genes_in_each_category_in_this_chromosome
chromosomes_total_gene_number = {} #dictionary with chromosome_Number-->Total_Number_of_genes_n_this_chromosome
for i in genes_dict.keys():
  chromosomes[i]={}
  for j in genes_dict[i].keys():
    chromosomes[i][j]={}

for i in fpkms[:,6]:
  chrm = str(i).split(":")[0].split('chr')[-1] #field 6 is locus in format "chr7:45567794-45589710"
  if chrm == 'X':
	chrm = "20"
  if chrm in chromosomes_total_gene_number.keys():
    chromosomes_total_gene_number[chrm] += 1
  else:
    chromosomes_total_gene_number[chrm] = 1.0 

all_chrms=[]
for i in genes_dict.keys():
  for j in genes_dict[i].keys():
    for k in genes_dict[i][j]:
      fpkms_idx = np.where(fpkms[:,0]==k)
      chrm = str(fpkms[fpkms_idx,6]).split(":")[0].split('chr')[-1] #field 6 is locus in format "chr7:45567794-45589710"
      if chrm == 'X':
	chrm = "20"
  
      if chrm in chromosomes[i][j].keys():
	chromosomes[i][j][chrm] += 1.0/chromosomes_total_gene_number[chrm]
	all_chrms.append(chrm)
      else:
	chromosomes[i][j][chrm] = 1.0/chromosomes_total_gene_number[chrm]

all_chrms.sort(key=lambda item: (len(item), item))
all_chrms=np.unique(np.array(all_chrms))

for i in genes_dict.keys():
  for j in genes_dict[i].keys():
    for k in all_chrms:
      if not k in chromosomes[i][j]:
	chromosomes[i][j][k] = 0

print "Saving picture"
#plot figures      
fig = plt.figure()
ax = plt.subplot(111)
markers = ["o","v","x","+"]
legends={"proper":"Subgroup 1","not_reprogrammed":"Subgroup 2","intermediate":"Subgroup 3","new_exp_profile":"Subgroup 4"}
for ind,i in enumerate(esc+fibs):
  plots=[]
  Legend_text=[]
#  for j in ["proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile"]:
  for j in ["proper","not_reprogrammed","intermediate","new_exp_profile"]:
    x_values = range(len(all_chrms))
    y_values = [chromosomes[j][i][v] for v in all_chrms]
    plots.append(ax.plot(x_values,y_values,color = colors[j],ls='-',marker=markers[ind]))
    Legend_text.append(i+" "+legends[j])
  
plt.ylabel("% of genes")
plt.xticks(np.arange(0,len(all_chrms)))
ax.set_xticklabels([i for i in all_chrms],rotation=90)
# Shink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.4,
               box.width, box.height * 0.6])
ax.set_xlabel('Chromosome')
#  ax.xaxis.set_label_coords(x=0,y=0)  
# Put a legend below current axis
ax.legend([p[0] for p in plots], Legend_text,loc='upper center', \
#  ax.legend([p[0] for p in plots], ["Reprogrammed","Not Reprogrammed","Intermediate state(NS)","Intermediate state(S)","Novel expression pattern"],loc='upper center', \    
	bbox_to_anchor=(0.5, -1*box.height*0.35),fancybox=True, shadow=True, ncol=3)
plt.savefig("Distribution_by_chromosomes.pdf",dpi=600)
plt.clf()
  

print "Stage 6. What happens with genes that have same expreesion between m5s and tauGFP"
genes_dict["changed_only_in_hybrids"]={}
for i in fibs+esc:
  print i
  diff_from_m5s = ma.masked_array(exp[:,0],getMask(exp,i,m5s,'yes',FPKM_threshold=FPKM_threshold_N,fold_change=-1)).compressed()
  diff_from_tauGFP = ma.masked_array(exp[:,0],getMask(exp,i,tauGFP,'yes',FPKM_threshold=FPKM_threshold_N,fold_change=-1)).compressed()

  differ_from_both_parents = np.intersect1d(diff_from_m5s,diff_from_tauGFP)
  differ_from_both_parents = np.intersect1d(differ_from_both_parents,same_in_tauGFP_and_m5s)
  print "Genes that change there there status: ",len(differ_from_both_parents)," out of ",len(same_in_tauGFP_and_m5s)," genes that are same in tauGFP and m5s"
  genes_dict["changed_only_in_hybrids"][i]=differ_from_both_parents

#Plot with ratio of genes changed_only_in_hybrids

plt_index = np.arange(len(esc+fibs))
ax = plt.subplot(111)
#for j in ["proper","not_reprogrammed","intermediate_nsign","intermediate_sign","new_exp_profile"]:
values = np.array([100.0/len(same_in_tauGFP_and_m5s)*len(genes_dict["changed_only_in_hybrids"][i]) for i in sorted(fibs+esc)])

p1 = ax.bar(plt_index, 100.0-values, plt_width, color=colors.values()[1])
p2 = ax.bar(plt_index, values, plt_width, color=colors.values()[0],bottom=100-values)

plt.ylabel('% of genes')
plt.xticks(plt_index+plt_width/2., sorted(fibs+esc))
plt.yticks(np.arange(0,101,10))

# Shink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])

# Put a legend below current axis
ax.legend([p1[0],p2[0]], ["Expression pattern is not changed","Expression pattern is changed"],loc='upper center', \
	bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=1)
plt.savefig("Genes_classification_2.pdf",dpi=600)
plt.clf()


  
#Plot with fold change of genes changed_only_in_hybrids
genes_fold_change["changed_only_in_hybrids"]={}
for i in fibs+esc: genes_fold_change["changed_only_in_hybrids"][i] = {}
for j in genes_dict["changed_only_in_hybrids"]:
    #determine parent
    if j in fibs:
      parent = m5s
    else:
      parent = tauGFP
    for k in genes_dict["changed_only_in_hybrids"][j]:
      #determine frkpm for the gene and parent
      fpkms_idx = np.where(fpkms[:,0]==k)
      parent_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[parent]][0][0])
      sample_fpkm = float(fpkms[fpkms_idx,fpkm_file_filds_N[j]][0][0])
      #count
      if parent_fpkm != 0 and sample_fpkm != 0:
	fold_change = max(parent_fpkm,sample_fpkm)/min(parent_fpkm,sample_fpkm)
	fold_change = round (fold_change)
	if (fold_change > 10):
	  fold_change = 10.0
	  print j,k,str(fpkms[fpkms_idx,6][0][0]).split(":")[0]
	if fold_change in genes_fold_change["changed_only_in_hybrids"][j].keys():
	  genes_fold_change["changed_only_in_hybrids"][j][fold_change] += 1.0
	else:
	  genes_fold_change["changed_only_in_hybrids"][j][fold_change] = 1.0
	  
#transform frequences into %
for j in genes_dict["changed_only_in_hybrids"]:
    total_genes = sum (genes_fold_change["changed_only_in_hybrids"][j].values())
    for k in genes_fold_change["changed_only_in_hybrids"][j].keys():
      genes_fold_change["changed_only_in_hybrids"][j][k] = (genes_fold_change["changed_only_in_hybrids"][j][k]/total_genes)*100
      
#plot figure
ax = plt.subplot(111)
plots=[]
for ind,i in enumerate(sorted(fibs+esc)):
    x_values = sorted(genes_fold_change["changed_only_in_hybrids"][i].keys())
    y_values = [genes_fold_change["changed_only_in_hybrids"][i][v] for v in sorted(genes_fold_change["changed_only_in_hybrids"][i].keys())]
    plots.append(ax.plot(x_values,y_values,color = sorted(colors.values())[ind],ls='-',marker='o'))
  
plt.ylabel('% of genes')
plt.xticks(np.arange(0,11))
plt.yticks(np.arange(0,101,5))
labels = [str(t)+"-"+str(t+1) for t in xrange(0,11)]
labels[0] = ""
labels[1] = "<2"
labels[-1] = ">10"
ax.set_xticklabels(labels)
 
# Shink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.3,
               box.width, box.height * 0.7])
ax.set_xlabel('Fold change between expression in parent and hybrid')
# Put a legend below current axis
ax.legend([p[0] for p in plots], sorted(fibs+esc),loc='upper center', \
  bbox_to_anchor=(0.5, -1*box.height*0.2),fancybox=True, shadow=True, ncol=2)
plt.savefig("genes_different_only_in_hybrids_fold_change.pdf",dpi=600)
plt.clf()

##########################################Save genes to file#######################
f=open("Genes.txt","w")
for ind,i in enumerate(fpkms[:,0]):
  s=str(i)+"\t"
  s += fpkms[ind,6] #chromosome
  for j in esc+fibs:
    category="-"
    for k in genes_dict.keys():
      if (i in genes_dict[k][j]):
	category = k
	break
    s+= "\t"+str(j)+"\t"+str(category)
  f.write(s+"\n")
f.close()
