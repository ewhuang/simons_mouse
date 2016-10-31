### Author: Edward Huang

import os
import subprocess
import sys
import time

### This script runs Sheng's PROSNET on a given network.
#!/bin/sh

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s data_type run_num' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    run_num = sys.argv[2]
    assert run_num.isdigit()

    results_folder = './Sheng/prosnet/result/prosnet_%s' % data_type
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    for method in ['go', 'no_go']:
        command = ('./Sheng/prosnet/model/embed -node ./data/%s_data/'
                    'dca_networks_%s/dca_genes_%s_%s.txt -link ./data/'
                    '%s_data/dca_networks_%s/dca_edges_%s_%s.txt -output '
                    '%s/dca_genes_%s_%s.vec -binary 0 -size 100 -negative 5'
                    '-samples 1 -iters 100 -threads 24 -model 2 -depth 10 '
                    '-restart 0.8 -edge_type_num 3 -rwr_ppi 1 -rwr_seq 1 '
                    '-train_mode 1' % (data_type, method, method, run_num,
                        data_type, method, method, run_num, results_folder,
                        method, run_num))
        # command = './Sheng/prosnet/model/embed -node ./data/%s_data/' % (
        #     data_type)
        # command += 'dca_networks_%s/dca_genes_%s_%s.txt ' % (method, method,
        #     run_num)
        # command += '-link ./data/%s_data/dca_networks_%s/dca_edges_%s_%s.txt ' % (
        #     data_type, method, method, run_num)
        # command += '-output ./Sheng/prosnet/result/dca_genes_%s_%s.vec ' % (
        #     method, run_num)
        # command += '-binary 0 -size 100 -negative 5 -samples 1 -iters 100 '
        # command += '-threads 24 -model 2 -depth 10 -restart 0.8 -edge_type_num '
        # command += '3 -rwr_ppi 1 -rwr_seq 1 -train_mode 1'
        subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))