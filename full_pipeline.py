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
    assert data_type in ['mouse', 'prosnet_mouse'] or data_type.isdigit()
    assert objective_function in ['wlogv', 'schaeffer', 'oclode']
    assert run_num.isdigit()

    # Create clustering input.
    command = 'python create_clustering_input.py %s %s' % (data_type, run_num)
    subprocess.call(command, shell=True)

    #TODO. Runs prosnet on the clustering input, and then reconstructs the
    # co-expression network. Currently not running prosnet.
    if data_type == 'prosnet_mouse':
        print 'Currently not re-running prosnet.'
        # command = 'python low_dimensional_nodes_prosnet.py mouse %s' % run_num
        # subprocess.call(command, shell=True)

        command = 'python create_prosnet_clustering_input.py %s %s' % (
            data_type, run_num)
        subprocess.call(command, shell=True)

    clus_format_str = (data_type, objective_function, run_num)
    # Perform simulated annealing.
    command = 'python simulated_annealing.py %s %s %s no_go' % clus_format_str
    subprocess.call(command, shell=True)

    command = 'python simulated_annealing.py %s %s %s go 0' % clus_format_str
    subprocess.call(command, shell=True)

    # Evaluate clusters.
    command = 'python evaluate_clustering.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)

    # Computing GO enrichments.
    command = 'python compute_go_enrichment.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)

    # Computing DBGAP enrichments.
    command = 'python compute_dbgap_enrichment.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)

    # Creating cluster summarization tables.
    command = 'python cluster_info_summary.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)

    # Plotting.
    command = 'python plot_best_clusters.py %s %s' % (data_type, run_num)
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))