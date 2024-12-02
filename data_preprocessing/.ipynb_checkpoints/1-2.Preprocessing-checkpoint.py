#!/usr/bin/env python
# coding: utf-8

# In[1]:


DBPATH = "/home/yc2553/projects/TRE_directionality/TRE_Directionality_DB/"
DBNAME = "2024_TRE_Directionality.db"
DBRESDIR = "RMSResources/"
import sys
import logging
logger = logging.getLogger("rmspool")
logger.setLevel("INFO")
logger.addHandler(logging.StreamHandler(sys.stdout))
from rmsp import rmscore
rms = rmscore.ResourceManagementSystem(f"{DBPATH}{DBNAME}", f"{DBPATH}{DBRESDIR}")
from rmsp import rmsbuilder
rmspool = rmsbuilder.RMSProcessWrapPool(rms, 8)
rmsb = rmsbuilder.RMSUnrunTasksBuilder(rmspool)
from rmsp import rmsutils
rms.set_scriptID("93b9b3d5f2934e639886ea387b0c4f50")



# In[ ]:


PROJECT_DIR_r = "/home/yc2553/projects/TRE_directionality/Resources/"
PROJECT_DIR_o = "/home/yc2553/projects/TRE_directionality/output/"
PROJECT_DIR_d = "/home/yc2553/projects/TRE_directionality/PROcap/"


# In[2]:


from commonhelper import load_func
RMSTemplate_preprocessing_with_UMI_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_preprocessing_with_UMI_20240501")
RMSTemplate_merge_reseqs_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_merge_reseqs_20240501")
RMSTemplate_merge_replicates_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_merge_replicates_20240501")
RMSTemplate_peak_calling_pints_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_peak_calling_pints_20240501")
RMSTemplate_generate_PROcap_gbratio_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_generate_PROcap_gbratio_20240501")
RMSTemplate_generate_replicates_correlation_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_generate_replicates_correlation_20240501")
RMSTemplate_generate_PROcap_RNA_length_20240501 = load_func("/fs/cbsuhy01/storage/kl945/github_repositories/Test3/TestDB/RMSLibrary/PROcapAnalysis/20240501/scripts/procapanalysis.py", "RMSTemplate_generate_PROcap_RNA_length_20240501")


# # PRO-cap preprocessing

# In[3]:


import os
import aldentools.methods as apm
from commonhelper import convert_to_bool as ctb
import pandas as pd
procap_info = pd.read_excel(f"{PROJECT_DIR_o}supp_tabels/suppTable1.PROcap_metainfo.xlsx", engine='openpyxl')


# In[ ]:


seq_lib_names = ["C1a", "C1b", "CTCF_U1", "CTCF_U2", "CTCF_T1", "CTCF_T2"]


# In[5]:


gencode_annotations = {"Human":f"{PROJECT_DIR_r}Gencode/gencode.v37.annotation.gff3.gz"}
chrom_size_files = {"Human":f"{PROJECT_DIR_r}Genome/hg38.chrom.sizes.filtered"}
star_nospikein_ref_dirs = {"Human":f"{PROJECT_DIR_r}STAR_genome_index/hg38_rDNA_20240515"}
star_spikein_ref_dirs = {"Human":f"{PROJECT_DIR_r}STAR_genome_index/hg38_rDNA_dm6_20240519"}


# ## Preprocessing

# In[6]:


df = procap_info[procap_info["Sequenced Library Name"].isin(seq_lib_names)]
files = []
for _, r in df.iterrows():
	rna_strand = "reverse"
	
	seq_lib_name = r["Sequenced Library Name"]
	r1_name = r["Sequencing Read 1"]
	r2_name = r["Sequencing Read 2"]
	
	if r["3' adapter"] == "VRA3 + 6N" and r["5' adapter"] == "VRA5 + 6N":
		adapter_seq1 = "TGGAATTCTCGGGTGCCAAGG"
		adapter_seq2 = "GATCGTCGGACTGTAGAACTCTGAAC"
		umi_loc='per_read'
		umi_len=6
	elif r["3' adapter"] == "VRA3 + 6N" and r["5' adapter"] == "VRA5":
		adapter_seq1 = "TGGAATTCTCGGGTGCCAAGG"
		adapter_seq2 = "GATCGTCGGACTGTAGAACTCTGAAC"
		umi_loc='read1'
		umi_len=6
	else:
		print("ERROR:", r["3' adapter"], r["5' adapter"], seq_lib_name)
		continue
	try:
		if ctb(r["S2 spike-in"].strip()):
			star_genome_dir = star_spikein_ref_dirs[r["Species"]]
		else:
			star_genome_dir = star_nospikein_ref_dirs[r["Species"]]
	except:
		print("ERROR:", r["S2 spike-in"], seq_lib_name)
		continue
	try:
		rms.file_from_path(f"{PROJECT_DIR_d}RawData/{r1_name}")
		rms.file_from_path(f"{PROJECT_DIR_d}RawData/{r2_name}")
	except:
		print("ERROR:", "Missing files", seq_lib_name)
		continue
	RMSTemplate_preprocessing_with_UMI_20240501(
		rmsb,
		**{"conda_env": "rms_procap_analysis_20240501",
						"project_dir": f"{PROJECT_DIR_d}",
						"prefix": f"{seq_lib_name}",
						"genome_dir": star_genome_dir,
						"chrom_size_file": chrom_size_files["Human"],
						"r1": f"{PROJECT_DIR_d}RawData/{r1_name}",
						"r2": f"{PROJECT_DIR_d}RawData/{r2_name}",
						"rna_strand": rna_strand,
						"umi_len": umi_len,
						"umi_loc": umi_loc,
						"adapter_seq_r1": adapter_seq1,
						"adapter_seq_r2": adapter_seq2,
						"thread": 32,
					   "remove_spliced_reads":True
		  }
	)


# In[7]:


sample_to_namedreplicates = {
	"C1": ["C1a", "C1b"],
	"CTCF_U": ["CTCF_U1", "CTCF_U2"],
	"CTCF_T": ["CTCF_T1", "CTCF_T2"],
}


# In[8]:


for sample, namedreplicates in sample_to_namedreplicates.items():
	RMSTemplate_merge_replicates_20240501(
		rmsb,
		**{"conda_env": "rms_procap_analysis_20240501",
			"project_dir": f"{PROJECT_DIR_d}",
			"prefixes": namedreplicates,
			"chrom_size_file": f"{PROJECT_DIR_r}Genome/Human/hg38.chrom.sizes.filtered",
		  }
	)


# In[6]:


sample_to_namedsample = {
	"C1": "brm_C1a_and_C1b_erm",
	"CTCF_U": "brm_CTCF_U1_and_CTCF_U2_erm",
	"CTCF_T": "brm_CTCF_T1_and_CTCF_T2_erm",
}


# In[7]:


for sample, namedsample in sample_to_namedsample.items():
	RMSTemplate_peak_calling_pints_20240501(
		rmsb,
		**{"conda_env": "rms_procap_analysis_20240501",
			"project_dir": f"{PROJECT_DIR_d}",
			"prefix": namedsample,
			"exp_type":"PROcap", 
			"thread": 32,
		  }

	)


# ## QC

# ### Gene body ratio

# In[11]:


for seq_lib_name in seq_lib_names:
	RMSTemplate_generate_PROcap_gbratio_20240501(
		rmsb,
		**{
			"conda_env": "rms_procap_analysis_20240501",
			"project_dir": f"{PROJECT_DIR_d}",
			"gbratio_annotation_file":f"{PROJECT_DIR_d}QC/Annotations/gencode.v37.annotation_proteincoding_nonoverlapping_long_transcripts.gtf.gz",
			"prefix": seq_lib_name,
		  }
	)


# ### Median RNA length

# In[12]:


for seq_lib_name in seq_lib_names:
	RMSTemplate_generate_PROcap_RNA_length_20240501(
		rmsb,
		conda_env='rms_procap_analysis_20240501',
		project_dir = f"{PROJECT_DIR_d}",
		prefix = seq_lib_name,
	)


# ### Replicates correlation

# In[14]:


for s, replicates in sample_to_namedreplicates.items():
	if len(replicates) <= 1:
		continue
	RMSTemplate_generate_replicates_correlation_20240501(
		rmsb,
		conda_env='rms_procap_analysis_20240501',
		project_dir = f"{PROJECT_DIR_d}",
		prefixes = replicates,
		bin_size=50,
		chrom_size_file=f"{PROJECT_DIR_r}Genome/Human/hg38.chrom.sizes.filtered.no.chrY"
	)
	


# In[15]:


rmspool.nthread=16
u = rmsb.execute_builder()


# In[23]:


rmspool.print_stat()

