### Author: Edward Huang

import file_operations
import subprocess
import sys
import time

### Run time: About half an hour.

def get_tcga_disease_list():
    disease_list = []
    f = open('./data/tcga_data/tcga_diseases.txt', 'r')
    for line in f:
        disease_list += [line.strip()]
    f.close()
    return disease_list

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s mouse/tcga_cancer_index objective_function run_num' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    objective_function = sys.argv[2]
    assert objective_function in ['wlogv', 'schaeffer', 'oclode', 'prosnet']
    run_num = sys.argv[3]
    assert run_num.isdigit()

    # Create clustering input.
    if objective_function != 'prosnet':
        command = 'python create_clustering_input.py %s %s' % (data_type, run_num)
        subprocess.call(command, shell=True)

        # Perform simulated annealing.
        command = 'python simulated_annealing.py %s %s %s no_go' % (data_type,
            objective_function, run_num)
        subprocess.call(command, shell=True)

        command = 'python simulated_annealing.py %s %s %s go 0' % (data_type,
            objective_function, run_num)
        subprocess.call(command, shell=True)
    # PROSNET has its own clustering system.
    else:
        command = 'python prosnet_kmeans.py %s %s' % (data_type, run_num)
        subprocess.call(command, shell=True)        

    # Evaluate clusters.
    command = 'python evaluate_clustering.py %s %s %s' % (data_type,
        objective_function, run_num)
    subprocess.call(command, shell=True)

    # Computing GO enrichments.
    command = 'python compute_go_enrichment.py %s %s %s' % (data_type,
        objective_function, run_num)
    subprocess.call(command, shell=True)

    # Creating cluster summarization tables.
    command = 'python cluster_info_summary.py %s %s %s' % (data_type,
        objective_function, run_num)
    subprocess.call(command, shell=True)

    # Plotting.
    command = 'python plot_best_clusters.py %s %s' % (data_type, run_num)
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))