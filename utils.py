def get_distance_to_center(inputfiles, outputfile):
	"""
	Display elements in heatmap in the following order 
	- Divergent: sort by the distance between two prominent TSSs
	- Unidirectional: sort by the distance between the prominent TSS and the center of overlapping DNase peaks
	"""

	from biodata.delimited import DelimitedReader
	import json

	distances = {}
	for inputfile in inputfiles:
		with DelimitedReader(inputfile) as dr:
			for cols in dr:
				chrom, start, end, start2, end2, tss1, tss2, sign = cols[:8]
				center = (int(end)+int(start))//2
				distances["_".join([chrom, start, end])] = abs(int(tss1) - int(center))
	with open(outputfile, "w") as f:
		json.dump(distances, f)
		
		

def generate_feature_metaplot(df, palette, hue_order, ax, test=True, errorbar="ci"):
	"""
	Generate a metaplot for a given feature. 
	"""

	import pandas as pd
	import seaborn as sns
	
	if test:
		frames = []
		for n in range(len(hue_order)):
			frames.append(df[df["Label"]==hue_order[n]].head(10))
		df = pd.concat(frames)

	# lineplot drops NAs from the DataFrame before plotting
	sns.lineplot(data=df, x="Position", y="Feature", hue="Label", hue_order=hue_order, palette=palette, errorbar=errorbar, ax=ax)



def generate_feature_heatmap(df, vlims, cmap, cbar, cbar_ax, cbar_kws, ax, sort_file=None, test=True, yticklabels=False):
	"""
	Generate a heatmap for a given feature. 
	"""

	import json
	import seaborn as sns

	if sort_file:
		with open(sort_file, "r") as f:
			sort_dict = json.load(f)
		df["sort_value"] = df.index.map(sort_dict)
		df = df.sort_values("sort_value", ascending=False)
		df = df.drop(columns="sort_value")

	if test:
		df = df.head(10)
					
	# Cells with missing values are automatically masked (will not be shown).
	sns.heatmap(df, vmin=vlims[0], vmax=vlims[-1], cmap=cmap, cbar=cbar, cbar_ax=cbar_ax, cbar_kws=cbar_kws, xticklabels=False, yticklabels=yticklabels, ax=ax)



def bin_values(df, bin_size=10):
    """
    Use bins to merge 5' end signals within bins as they're very sparse
	Bins are tiled from the anchor point outward using a (-5, +4) scheme for bin_size=10,
    such that the anchor bin (position 0) spans -half_bin to +(half_bin-1) bp,
    and each flanking position ±n is 10n bp from the anchor.
    """
    
    import pandas as pd
    import numpy as np
    
    results = []
    center = len(df.columns) // 2
    half_bin = bin_size // 2
    n_flank_bins = (center - half_bin) // bin_size
    
    for index, row in df.iterrows():
        new_row = []
        # left flank (outermost to innermost)
        for n in range(n_flank_bins, 0, -1):
            start = center - half_bin - n * bin_size
            new_row.append(row.iloc[start:start+bin_size].sum())
        
        # center bin: spans center-half_bin to center+(half_bin-1)
        new_row.append(row.iloc[center-half_bin:center+half_bin].sum())
        
        # right flank
        for n in range(n_flank_bins):
            start = center + half_bin + n * bin_size
            new_row.append(row.iloc[start:start+bin_size].sum())
        
        results.append(new_row)
    
    df_bins = pd.DataFrame(results, index=df.index)
    
    return df_bins



#--------------------------------------------------------------------------------------------------------
# Motif logo
#--------------------------------------------------------------------------------------------------------

def ppm2im(ppm):
	"""
	Goal:
	Convert position probability matrix to information matrix (bits).
	"""

	import numpy as np
	import pandas as pd
	
	exp = 4*(-0.25)*np.log2(0.25)
	results = []
	for index, row in ppm.iterrows():
		obs = sum([-row[c]*np.log2(row[c]) if row[c] > 0 else 0 for c in ppm.columns])
		results.append([(exp-obs)*row[c] for c in ppm.columns])
	df = pd.DataFrame(results, index=ppm.index, columns=ppm.columns)
	return df



def read_meme(filename):
	import numpy as np
	
	pwms = {}
	motif_name = None
	pwm_lines = []
	
	def save_pwm(name, lines):
		if lines:
			# only keep lines with exactly 4 floats
			clean_lines = [line for line in lines if len(line) == 4]
			if clean_lines:
				pwms[name] = np.array(clean_lines, dtype=float)
	
	with open(filename) as f:
		for line in f:
			line = line.strip()
			if line.startswith("MOTIF"):
				# save previous motif
				save_pwm(motif_name, pwm_lines)
				motif_name = line.split()[1]
				pwm_lines = []
			elif line.startswith("letter-probability matrix"):
				pwm_lines = []
			else:
				parts = line.split()
				# only numeric lines with 4 values
				if len(parts) == 4:
					try:
						pwm_lines.append([float(x) for x in parts])
					except ValueError:
						pass
		# save last motif
		save_pwm(motif_name, pwm_lines)
	
	return pwms



# Refer to codes from https://github.com/kundajelab/ProCapNet/blob/main/src/figure_notebooks/other_motif_utils.py and https://github.com/kundajelab/ProCapNet/blob/main/src/figure_notebooks/Fig2_modisco_results.ipynb
def compute_per_position_ic(ppm, pseudocount=0.001):
	import numpy as np
	
	background = np.array([0.25] * 4)
	alphabet_len = len(background)
	
	ppm_with_pseudocount = (ppm+pseudocount)/(1 + pseudocount*alphabet_len)
	ppm_logodds = np.log(ppm_with_pseudocount) * ppm / np.log(2)
	background_logodds = np.log(background) * background / np.log(2)
	ic = (ppm_logodds - background_logodds[None,:])
	return np.sum(ic, axis=1)



def trim_motif_by_thresh(pfm, trim_threshold=0.3, pad=2):
	import numpy as np
	
	trim_thresh = np.max(pfm) * trim_threshold
	pass_inds = np.where(pfm >= trim_thresh)[0]
	
	start = max(np.min(pass_inds) - pad, 0)
	end = min(np.max(pass_inds) + pad + 1, len(pfm) + 1)
	
	return start, end



def plot_motif_on_ax(array, ax, title=None, visible=False, fontsize=12):
	import pandas as pd
	import logomaker
	
	assert len(array.shape) == 2 and array.shape[-1] == 4, array.shape
	df = pd.DataFrame(array, columns=['A', 'C', 'G', 'T'])
	df.index.name = 'pos'
	
	crp_logo = logomaker.Logo(df, ax=ax, baseline_width=0)
	
	ax.set_ylim(min(df.sum(axis=1).min(), 0), df.sum(axis=1).max())
	ax.set_xticks([])
	if not visible:
		crp_logo.style_spines(visible=False)
		ax.set_yticks([])
	
	if title:
		ax.set_title(title, fontsize=fontsize)
	
	return crp_logo



def plot_PROcap(array, ax, pcolor, ncolor, width=2):	
	import pandas as pd
	
	orientations = {"fwd": pcolor, "rev": ncolor}
	length = array.shape[1]
	results = [[i, array[0][i], "fwd"] for i in range(length)]
	results += [[i, -array[1][i], "rev"] for i in range(length)]		
		
	df = pd.DataFrame(results, columns=["position", "reads", "orientation"])
	
	for orientation, color in orientations.items():
		df_orientation = df[df["orientation"] == orientation]
		ax.bar(df_orientation["position"], df_orientation["reads"], color=color, align="edge", width=width)



def plot_contrib_scores(array, ax):
	import pandas as pd
	import logomaker
	
	df = pd.DataFrame(array[:, 250:750].T, columns=['A', 'C', 'G', 'T'])
	color_dict = {'A': '#62a479', 'C': '#7184be', 'G': '#f9be6b', 'T': '#f36a69'}
	crp_logo = logomaker.Logo(df, ax=ax, color_scheme=color_dict, baseline_width=0)



#--------------------------------------------------------------------------------------------------------
# Motif enrichment (HOMER)
#--------------------------------------------------------------------------------------------------------

def generate_homer_input(es, outputfile):
	"""
	Generate input files conforming to HOMER's format
	"""
	from biodata.delimited import DelimitedWriter
	
	with DelimitedWriter(outputfile) as dw:
		for e in es:
			chrom, start, end = e.split("_")
			dw.write([chrom, start, end, "_".join([chrom, start, end]), ".", "+"])



def run_homer(homer_dir, target, outdir, bg=None, motif_file=None, denovo=True):
	"""
	Run HOMER for motif enrichment
	"""
	import subprocess
	from os.path import exists
	
	# HOMER Motif Database (http://homer.ucsd.edu/homer/motif/motifDatabase.html): This database is maintained as part of HOMER and is mostly based on the analysis of public ChIP-Seq data sets. These motifs are often referred to in the HOMER software as 'known' motifs since their degeneracy thresholds have been optimized by HOMER, unlike motifs found in JASPAR or other public databases.
	# http://homer.ucsd.edu/homer/ngs/peakMotifs.html: If custom background regions are provided ("-bg <peak/BED file>"), HOMER will automatically ensure that these regions do NOT overlap with the target regions (using mergePeaks). Custom regions will still be normalized for GC-content.
	# http://homer.ucsd.edu/homer/introduction/basics.html: The greatest advantage to using known motifs is found when you have a limited set of target sequences. The less data that is available or the weaker the true signal, it is difficult for de novo motif finding to accurately define a signal that is significant. Known motifs have the advantage of many less degrees of freedom and in may cases find the correct motifs when the enrichment falls below the 1e-10 thresholds for reliability when considering de novo results.

	if exists(outdir):
		subprocess.run(f"rm -r {outdir}", shell=True)
		
	commands = ["findMotifsGenome.pl",
				 target,
				 "hg38",
				 outdir,
				 "-size given"
				]
	if bg:
		commands.extend(["-bg", bg])
	
	if not denovo:	
		commands.append("-nomotif")
	
		
	if motif_file:
		commands.append(f"-mknown {motif_file}")
		
	subprocess.run(" ".join(commands), shell=True)



