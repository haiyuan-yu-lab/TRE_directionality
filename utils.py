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
	"""

	import pandas as pd
	import numpy as np
		
	results = []
	ns = (len(df.columns) - 1) // bin_size
	for index, row in df.iterrows():
	    new_row = []
	    for n in range(ns):
	        new_row.append(row[n*bin_size:(n+1)*bin_size].sum())
	    new_row.append(row.iloc[-1])
	    results.append(new_row)
	df2 = pd.DataFrame(results, index=df.index)

	return df2



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



# Get codes from Taskiran_et_al_Supplemental_Code/Fly/utils.py
def plot_a(ax, base, left_edge, height, color):

	import numpy as np
	import matplotlib
	
	a_polygon_coords = [
		np.array([
			[0.0, 0.0],
			[0.5, 1.0],
			[0.5, 0.8],
			[0.2, 0.0],
		]),
		np.array([
			[1.0, 0.0],
			[0.5, 1.0],
			[0.5, 0.8],
			[0.8, 0.0],
		]),
		np.array([
			[0.225, 0.45],
			[0.775, 0.45],
			[0.85, 0.3],
			[0.15, 0.3],
		])
	]

	for polygon_coords in a_polygon_coords:
		ax.add_patch(matplotlib.patches.Polygon((np.array([1, height])[None, :] * polygon_coords + np.array([left_edge, base])[None, :]), facecolor=color, edgecolor=color))



# Get codes from Taskiran_et_al_Supplemental_Code/Fly/utils.py
def plot_c(ax, base, left_edge, height, color):

	import numpy as np
	import matplotlib

	ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height, facecolor=color, edgecolor=color))
	ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height, facecolor='white', edgecolor='white'))
	ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height, facecolor='white', edgecolor='white', fill=True))



# Get codes from Taskiran_et_al_Supplemental_Code/Fly/utils.py
def plot_g(ax, base, left_edge, height, color):
	
	import numpy as np
	import matplotlib
	
	ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height, facecolor=color, edgecolor=color))
	ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height, facecolor='white', edgecolor='white'))
	ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height, facecolor='white', edgecolor='white', fill=True))
	ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.825, base + 0.085 * height], width=0.174, height=0.415 * height, facecolor=color, edgecolor=color, fill=True))
	ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.625, base + 0.35 * height], width=0.374, height=0.15 * height, facecolor=color, edgecolor=color, fill=True))



# Get codes from Taskiran_et_al_Supplemental_Code/Fly/utils.py
def plot_t(ax, base, left_edge, height, color):

	import numpy as np
	import matplotlib
	
	ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.4, base], width=0.2, height=height, facecolor=color, edgecolor=color, fill=True))
	ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge, base + 0.8 * height], width=1.0, height=0.2 * height, facecolor=color, edgecolor=color, fill=True))



# Get codes from Taskiran_et_al_Supplemental_Code/Fly/utils.py
default_colors = {0: 'green', 1: 'blue', 2: 'orange', 3: 'red'}
default_plot_funcs = {0: plot_a, 1: plot_c, 2: plot_g, 3: plot_t}
def plot_weights_given_ax(ax, array,
                          height_padding_factor,
                          length_padding,
                          subticks_frequency,
                          highlight,
                          colors=default_colors,
                          plot_funcs=default_plot_funcs):

	import numpy as np
	import matplotlib
	
	if len(array.shape) == 3:
		array = np.squeeze(array)
	assert len(array.shape) == 2, array.shape
	if array.shape[0] == 4 and array.shape[1] != 4:
		array = array.transpose(1, 0)
	assert array.shape[1] == 4
	max_pos_height = 0.0
	min_neg_height = 0.0
	heights_at_positions = []
	depths_at_positions = []
	for i in range(array.shape[0]):
		acgt_vals = sorted(enumerate(array[i, :]), key=lambda x: abs(x[1]))
		positive_height_so_far = 0.0
		negative_height_so_far = 0.0
		for letter in acgt_vals:
			plot_func = plot_funcs[letter[0]]
			color = colors[letter[0]]
			if letter[1] > 0:
				height_so_far = positive_height_so_far
				positive_height_so_far += letter[1]
			else:
				height_so_far = negative_height_so_far
				negative_height_so_far += letter[1]
			plot_func(ax=ax, base=height_so_far, left_edge=i, height=letter[1], color=color)
		max_pos_height = max(max_pos_height, positive_height_so_far)
		min_neg_height = min(min_neg_height, negative_height_so_far)
		heights_at_positions.append(positive_height_so_far)
		depths_at_positions.append(negative_height_so_far)
	
	for color in highlight:
		for start_pos, end_pos in highlight[color]:
			assert start_pos >= 0.0 and end_pos <= array.shape[0]
			min_depth = np.min(depths_at_positions[start_pos:end_pos])
			max_height = np.max(heights_at_positions[start_pos:end_pos])
			ax.add_patch(
				matplotlib.patches.Rectangle(xy=[start_pos, min_depth],
											 width=end_pos - start_pos,
											 height=max_height - min_depth,
											 edgecolor=color, fill=False))
	
	ax.set_xlim(-length_padding, array.shape[0] + length_padding)
	ax.xaxis.set_ticks(np.arange(0.0, array.shape[0] + 1, subticks_frequency))
	height_padding = max(abs(min_neg_height) * (height_padding_factor),
						 abs(max_pos_height) * (height_padding_factor))
	ax.set_ylim(min_neg_height - height_padding, max_pos_height + height_padding)
	return ax



# Modify codes from Taskiran_et_al_Supplemental_Code/Fly/utils.py
def plot_weights(array, ax,
                 height_padding_factor=0.2,
                 length_padding=1.0,
                 subticks_frequency=20,
                 colors=default_colors,
                 plot_funcs=default_plot_funcs,
                 highlight={}):

	import matplotlib
	
	y = plot_weights_given_ax(ax=ax, array=array,
                              height_padding_factor=height_padding_factor,
                              length_padding=length_padding,
                              subticks_frequency=subticks_frequency,
                              colors=colors,
                              plot_funcs=plot_funcs,
                              highlight=highlight)
	return ax



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



