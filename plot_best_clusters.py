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

def get_auc(point_list):
    '''
    Gets the AUC curve for a list of points. For each threshold of p-values (x-
    axis), y-axis is the number of points that has a lower p-value.
    '''
    auc_point_list = []
    # Get the sorted p-values in our list of points.
    sorted_point_list = sorted(point_list, key=lambda x: x[0])
    num_point_list = len(point_list) - len([pt for pt in point_list if pt[1] < indensity_threshold])

    for i, (p_value, density) in enumerate(sorted_point_list):
        if density < indensity_threshold:
            continue
        num_better = num_point_list - i
        auc_point_list += [(p_value, num_better)]
    return auc_point_list

def process_header_file(line):
    '''
    Get the indices of the relevant columns.
    '''
    global num_gene_idx, dbgap_enrichment_idx, gwas_enrichment_idx
    global go_enrichment_idx, ratio_idx, t_test_idx, labeled_mean_idx
    global unlabeled_mean_idx
    num_gene_idx = line.index('Number of genes')
    dbgap_enrichment_idx = line.index('Top DBGAP p-value')
    gwas_enrichment_idx = line.index('Top GWAS p-value')
    go_enrichment_idx = line.index('1st GO p-value')
    ratio_idx = line.index('In/(In+Out)')
    if 't-test p-value' in line:
        t_test_idx = line.index('t-test p-value')
        labeled_mean_idx = line.index('Labeled mean')
        unlabeled_mean_idx = line.index('Unlabeled mean')

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
        point_list = []

        if mode == 'wgcna':
            fname = ('./results/%s_results/wgcna/clus_info_go/clus_info_go_'
                '%s.tsv' % (data_type, run_num))
        elif mode == 'wlogv_go':
            fname = ('./results/%s_results/wlogv/clus_info_go/clus_info_go_'
                '%s.tsv' % (data_type, run_num))
            cheating_point_list = []
        elif mode == 'wlogv_no_go':
            fname = ('./results/%s_results/wlogv/clus_info_no_go/clus_info_'
                'no_go_%s.tsv' % (data_type, run_num))

        f = open(fname, 'r')
        for i, line in enumerate(f):
            line = line.strip().split('\t')
            if i == 0:
                process_header_file(line)
                continue
            # Skip the clusters with fewer than 30 genes.
            if int(line[num_gene_idx]) < 30:
                continue

            if plot_type == 'dbgap':
                enrichment_p = line[dbgap_enrichment_idx]
            elif plot_type == 'gwas':
                enrichment_p = line[gwas_enrichment_idx]
            elif 'go' in plot_type:
                enrichment_p = line[go_enrichment_idx]
            top_enrichment_p = -math.log(float(enrichment_p), 10)
            point = (top_enrichment_p, float(line[ratio_idx]))

            # Add the point to the cheating_point_list if the genes labeled by
            # the best enriched GO term have significantly worse means
            # than those not labeled.
            # TODO: Change these indices to find the mean columns.
            if mode == 'wlogv_go' and float(line[t_test_idx]) < 1e-5 and float(
                line[labeled_mean_idx]) < float(line[unlabeled_mean_idx]):
                cheating_point_list += [point]
            else:
                point_list += [point]
        f.close()

        if mode == 'wgcna':
            global indensity_threshold
            indensity_threshold = np.median([pt[1] for pt in point_list]) * 0.8

        # Plot AUC points.
        if 'auc' in plot_type:
            point_list = get_auc(point_list)

        # Update the largest y-axis value.
        highest_p = max(highest_p, max(pt[0] for pt in point_list))
        highest_y = max(highest_y, max(pt[1] for pt in point_list))
        num_high_points = len([p for p in point_list if p[0] >= 10])

        mode_label = '%s, %d' % (mode, num_high_points)
        if 'auc' in plot_type:
            plt.plot(*zip(*point_list), color=colors[mode_index],
                label=mode_label)
        else:
            plt.scatter(*zip(*point_list), color=colors[mode_index],
                label=mode_label, marker=mode_marker_list[mode_index])

        # Plot the cheating clusters.
        if mode == 'wlogv_go' and 'auc' not in plot_type:
            if len(cheating_point_list) == 0:
                continue
            num_high_points = len([p for p in cheating_point_list if p[0] >= 10])
            mode_label = mode + ', %d' % num_high_points
            plt.scatter(*zip(*cheating_point_list), color=colors[mode_index],
                label='cheating ' + mode_label, marker='x')
        # Horizontal line of 0.8 * median of WGCNA cluster ratios.
        if 'auc' not in plot_type:
            plt.axhline(indensity_threshold)

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