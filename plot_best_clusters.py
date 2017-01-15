### Author: Edward Huang

import file_operations
import math
import sys
import matplotlib
import numpy as np

### This script plots a scatterplot for four networks. Network without GO,
### network with GO for lambda=0.9, WGCNA clusters, and preliminary clusters
### from Sheng's prediction networks.

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

def get_auc(pts):
    '''
    Gets the AUC curve for a list of points. For each threshold of p-values (x-
    axis), y-axis is the number of points that has a lower p-value.
    '''
    auc_pts = []
    # Get the sorted p-values in our list of points.
    p_values = sorted([pt[0] for pt in pts])
    num_pts = len(pts)
    for i in range(len(p_values)):
        p_value = p_values[i]
        num_better = num_pts - i
        auc_pts += [(p_value, num_better)]
    return auc_pts

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: %s data_type run_num go/go_auc/dbgap/gwas' % sys.argv[0]
        exit()
    data_type, run_num, plot_type = sys.argv[1:]
    assert (data_type == 'mouse' or data_type.isdigit()) and run_num.isdigit()
    assert plot_type in ['go', 'dbgap', 'go_auc', 'gwas']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]
    print data_type

    colors, mode_marker_list = ['blue', 'red', 'black'], ['<', 'v', 'o']

    highest_p, highest_y = 0, 0
    for mode_index, mode in enumerate(['wgcna', 'wlogv_no_go', 'wlogv_go']):
        pts = []

        if mode == 'wgcna':
            fname = './wgcna/results/%s_results/clus_info.tsv' % data_type
        elif mode == 'wlogv_go':
            fname = ('./results/%s_results/wlogv/clus_info_go/clus_info_go_'
                '%s.tsv' % (data_type, run_num))
            cheating_pts = []
        elif mode == 'wlogv_no_go':
            fname = ('./results/%s_results/wlogv/clus_info_no_go/clus_info_'
                'no_go_%s.tsv' % (data_type, run_num))

        f = open(fname, 'r')
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.split()
            # Skip the clusters with fewer than 10 genes.
            if int(line[4]) < 30:
                continue

            if plot_type == 'dbgap':
                enrichment_p = line[15]
            elif plot_type == 'gwas':
                enrichment_p = line[17]
            else:
                enrichment_p = line[5]
            top_enrichment_p = -math.log(float(enrichment_p), 10)
            point = (top_enrichment_p, float(line[3]))

            # Add the point to the cheating_pts if the genes labeled by
            # the best enriched GO term have significantly worse means
            # than those not labeled.
            # TODO: Change these indices to find the mean columns.
            if mode == 'wlogv_go' and float(line[19]) < 1e-5 and float(line[20]
                ) < float(line[23]):
                cheating_pts += [point]
            else:
                pts += [point]
        f.close()

        # Plot AUC points.
        if 'auc' in plot_type:
            pts = get_auc(pts)

        # Update the largest y-axis value.
        highest_p = max(highest_p, max(pt[0] for pt in pts))
        highest_y = max(highest_y, max(pt[1] for pt in pts))
        num_high_points = len([p for p in pts if p[0] >= 10])

        mode_label = '%s, %d' % (mode, num_high_points)
        if 'auc' in plot_type:
            plt.plot(*zip(*pts), color=colors[mode_index], label=mode_label)
        else:
            plt.scatter(*zip(*pts), color=colors[mode_index], label=mode_label,
                marker=mode_marker_list[mode_index])

        # Plot the cheating clusters.
        if mode == 'wlogv_go' and 'auc' not in plot_type:
            if len(cheating_pts) == 0:
                continue
            num_high_points = len([p for p in cheating_pts if p[0] >= 10])
            mode_label = mode + ', %d' % num_high_points
            plt.scatter(*zip(*cheating_pts), color=colors[mode_index],
                label='cheating ' + mode_label, marker='x')
        # Horizontal line of 0.8 * median of WGCNA cluster ratios.
        elif mode == 'wgcna' and 'auc' not in plot_type:
            plt.axhline(np.median([pt[1] for pt in pts]) * 0.8)

    # Construct the plot details.
    plt.axvline(10)
    if 'auc' not in plot_type:
        plt.title('Cluster in-densities vs. Best enrichment p-values')
        plt.ylabel('in-density/(in-density + out-density)')
    else:
        plt.title('Number of clusters vs. Best enrichment p-values')
        plt.ylabel('Number of clusters')
    plt.xlabel('negative log of lowest GO enrichment p-value')
    
    if 'auc' in plot_type:
        plt.legend(loc='upper right')
    else:
        plt.legend(loc='lower right')
    plt.ylim(0, highest_y * 1.1)
    plt.xlim(0, highest_p * 1.2)
    plt.show()

    # Save extra plots in the compiled results folder.
    if plot_type in ['go', 'go_auc']:
        pylab.savefig('./results/plots_and_tables/%s_%s.png' % (data_type,
            plot_type))

    folder = './results/%s_results/comparison_plots' % data_type
    pylab.savefig('%s/%s_comparison_plot_%s.png' % (folder, plot_type, run_num))
    plt.close()