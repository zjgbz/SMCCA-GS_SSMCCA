import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

argvs = sys.argv
cca_seed = int(argvs[1])
feature_factor = int(argvs[2])

size_factor = 100
nrow = 4 * size_factor
ncol_rel = 5
noise_sd = 5
cv_num = 15
cca_fashion = "sup-mcca"
# matrix_seed_list = list(range(2, 21))
matrix_seed_list = [1]
noise_mode = "noise-reinf"

ncol = 8 * size_factor * feature_factor
rep = 1 * size_factor * feature_factor

group_dict = {"group1":np.arange(0, 1 * rep), "group2":np.arange(1 * rep, 2 * rep), "group3":np.arange(2 * rep, 3 * rep),
			"group4":np.arange(3 * rep, 4 * rep), "group5":np.arange(4 * rep, 5 * rep), "group_noise":np.arange(5 * rep, 8 * rep)}

for matrix_seed in matrix_seed_list:
	load_dir = os.path.join(".", f"matrix-seed{matrix_seed}_cca-seed{cca_seed}_size{size_factor}_dim{feature_factor}_{noise_mode}")
	save_dir = load_dir
	
	w_name_list = ["w1", "w2"]
	for w_name in w_name_list:
		w_filename_prefix = f"{cca_fashion}_cca-seed{cca_seed}_matrix-seed{matrix_seed}_nrow{nrow}_ncol-rel{ncol_rel}_rep{rep}_noise-sd{noise_sd}_{noise_mode}_ncol{ncol}_{w_name}"
		w_filename = f"{w_filename_prefix}.tsv"
		w_dir_filename = os.path.join(load_dir, w_filename)
		w = pd.read_csv(w_dir_filename, sep = "\t", header = None, index_col = None)
		
		axs_x_num = 4
		axs_y_num = 4
		fig, axs = plt.subplots(axs_x_num, axs_y_num, sharex=True, sharey=True, gridspec_kw = {'wspace':0., 'hspace':0.1}, figsize=(16,16))
		my_suptitle = fig.suptitle(f"{cca_fashion} {nrow} x {ncol}, {w_name}",x=0.5,y=0.92, fontsize=20)
		raw_ymax_list = []
		raw_ymin_list = []
		for axs_x in range(axs_x_num):
			for axs_y in range(axs_y_num):
				cv_i = axs_x + axs_y * axs_x_num
				if cv_i != axs_x_num * axs_y_num - 1:
					for group_i in group_dict:
						group_i_idx = group_dict[group_i]
						w_group_i = w.iloc[group_i_idx, cv_i].values.reshape(1, -1)[0]
						axs[axs_y, axs_x].scatter(group_i_idx, w_group_i, label = group_i)
						# axs[axs_y, axs_x].set_ylim([-0.0625, 0.0625])
						axs[axs_y, axs_x].text(0.85, 0.9, f'cv {cv_i + 1}', fontsize = 20, horizontalalignment='center', verticalalignment='center', transform=axs[axs_y, axs_x].transAxes)
						raw_ymax_i = np.amax(w_group_i)
						raw_ymin_i = np.amin(w_group_i)
						raw_ymax_list.append(raw_ymax_i)
						raw_ymin_list.append(raw_ymin_i)
		raw_ymax = max(raw_ymax_list)
		raw_ymin = min(raw_ymin_list)
		raw_ylim = max(abs(raw_ymax), abs(raw_ymin))
		custom_ylim = (-abs(raw_ymax) * 6 / 5, abs(raw_ymax) * 6 / 5)
		plt.setp(axs, ylim=custom_ylim)
		fig_filename = f"{w_filename_prefix}_manha_cv{cv_num}.png"
		fig_dir_filename = os.path.join(save_dir, fig_filename)
		plt.savefig(fig_dir_filename, bbox_inches='tight', dpi = 600, bbox_extra_artists=[my_suptitle])
		plt.close()