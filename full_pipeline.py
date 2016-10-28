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
    subprocess.call(command, shell=True)

    # Runs prosnet on the clustering input.
    command = 'python low_dimensional_nodes_prosnet.py %s %s' % (data_type,
        run_num)
    subprocess.call(command, shell=True)

    # Reconstruct network from ProsNet vectors.
    command = 'python create_prosnet_clustering_input.py %s %s' % (data_type,
        run_num)
    subprocess.call(command, shell=True)

    clus_format_str = (data_type, objective_function, run_num)
    prosnet_format = ('prosnet_' + data_type, objective_function, run_num)
    
    # Perform simulated annealing on coexpression network.
    command = 'python simulated_annealing.py %s %s %s no_go' % clus_format_str
    subprocess.call(command, shell=True)
    command = 'python simulated_annealing.py %s %s %s go 0' % clus_format_str
    subprocess.call(command, shell=True)

    # Perform simulated annealing on prosnet network.
    command = 'python simulated_annealing.py %s %s %s no_go' % prosnet_format
    subprocess.call(command, shell=True)
    command = 'python simulated_annealing.py %s %s %s go' % prosnet_format
    subprocess.call(command, shell=True)

    # Evaluate clusters.
    command = 'python evaluate_clustering.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)
    command = 'python evaluate_clustering.py %s %s %s' % prosnet_format
    subprocess.call(command, shell=True)

    # Computing GO enrichments.
    command = 'python compute_go_enrichment.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)
    command = 'python compute_go_enrichment.py %s %s %s' % prosnet_format
    subprocess.call(command, shell=True)

    # Computing DBGAP enrichments.
    command = 'python compute_dbgap_enrichment.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)
    command = 'python compute_dbgap_enrichment.py %s %s %s' % prosnet_format
    subprocess.call(command, shell=True)

    # Creating cluster summarization tables.
    command = 'python cluster_info_summary.py %s %s %s' % clus_format_str
    subprocess.call(command, shell=True)
    command = 'python cluster_info_summary.py %s %s %s' % prosnet_format
    subprocess.call(command, shell=True)

    # Plotting.
    command = 'python plot_best_clusters.py %s %s' % (data_type, run_num)
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))