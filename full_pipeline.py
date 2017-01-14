### Author: Edward Huang

import subprocess
import sys
import time

### Runs the full pipeline, from creating the co-expression files to plotting
### the results.
### Run time: About half an hour.

def main():
    if len(sys.argv) != 4:
        print ('Usage:python %s mouse/tcga_cancer_index objective_function '
                'run_num' % sys.argv[0])
        exit()
    data_type, objective_function, run_num = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    assert objective_function in ['wlogv', 'schaeffer', 'oclode']
    assert run_num.isdigit()

    # Create clustering input. Also creates files necessary for ProsNet.
    command = 'python create_clustering_input.py %s %s' % (data_type, run_num)
    print command
    subprocess.call(command, shell=True)

    clus_format_str = (data_type, objective_function, run_num)
    
    # Perform simulated annealing on coexpression network.
    for network_type in ['go', 'no_go']:
        command = 'python simulated_annealing.py %s %s %s ' % clus_format_str
        command += network_type
        print command
        subprocess.call(command, shell=True)

    # Evaluate clusters.
    command = 'python evaluate_clustering.py %s %s %s' % clus_format_str
    print command
    subprocess.call(command, shell=True)

    # Computing label enrichments.
    for label_type in ['go', 'dbgap']:
        command = 'python compute_label_enrichments.py %s %s %s ' % (
            clus_format_str)
        command += label_type
        print command
        subprocess.call(command, shell=True)

    # Determine if wlogv_go cheats.
    command = 'python cheating_evaluation.py %s %s %s' % clus_format_str
    print command
    subprocess.call(command, shell=True)

    # Creating cluster summarization tables.
    command = 'python cluster_info_summary.py %s %s %s' % clus_format_str
    print command
    subprocess.call(command, shell=True)

    # Plotting.
    for plot_type in ['go', 'go_auc', 'dbgap']:
        command = 'python plot_best_clusters.py %s %s %s' % (data_type, run_num,
            plot_type)
        print command
        subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("Full pipeline took %s seconds..." % (time.time() - start_time))