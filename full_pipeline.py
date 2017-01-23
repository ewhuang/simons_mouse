### Author: Edward Huang

from multiprocessing import Pool
import os
import subprocess
import sys
import time

### Run time: About 20 minutes.

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s mouse/tcga_idx wlogv/wgcna run_num' % sys.argv[0]
        exit()
    data_type, clus_method, run_num = sys.argv[1:]
    assert (data_type == 'mouse' or data_type.isdigit()) and run_num.isdigit()
    assert clus_method in ['wlogv', 'schaeffer', 'oclode', 'wgcna']

    str_format_lst = [data_type, clus_method, run_num]

    if clus_method != 'wgcna':
        command = 'create_clustering_input.py'
        print command
        subprocess.call(['python', command, data_type, run_num])

        # Multiple processes for simulated annealing.
        pool = Pool(processes=2)
        for network_type in ['go', 'no_go']:
            command = 'simulated_annealing.py'
            print command
            arg = ['python', command] + str_format_lst + [network_type]
            pool.apply_async(subprocess.call, (arg,))
        pool.close()
        pool.join()
        # Remove the log file.
        os.remove('log')

    # Evaluate clusters for all clusters.
    command = 'evaluate_clustering.py'
    print command
    subprocess.call(['python', command] + str_format_lst)

    # Don't pool this for loop, since we are pooling in the script already.
    for label_type in ['go', 'dbgap', 'gwas']:
        command = 'compute_label_enrichments.py'
        print command
        subprocess.call(['python', command] + str_format_lst + [label_type])

    if clus_method != 'wgcna':
        command = 'cheating_evaluation.py'
        print command
        subprocess.call(['python', command] + str_format_lst)

    # Creating cluster summarization tables.
    command = 'cluster_info_summary.py'
    print command
    subprocess.call(['python', command] + str_format_lst)

    # Don't bother pooling, since plotting is fast enough.
    if clus_method != 'wgcna':
        for plot_type in ['go', 'go_auc', 'dbgap', 'gwas']:
            command = 'plot_best_clusters.py'
            subprocess.call(['python', command, data_type, run_num, plot_type])

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("Full pipeline took %s seconds..." % (time.time() - start_time))