#include "cs-grn.h"
#include "ReadInputs.h"
#include "Cluster.h"
#include <thread>
using namespace std;

double abs(double x) {
    if (x < 0) return -x;
    return x;
}

// main
int main (int argc, char * const argv[]) {
    if (argc < 7) {
        printf("usage: %s <numclusters> <numtrials> <couplingConstant> <orthfile> <numSpc> <spc1nw> <spc2nw> ... -s <startclustering> -t <maxtemp>\n", argv[0]);
        exit(1);
    }
    clock_t begin = clock();
    int argbase = 1;
    // read numclusters, numTrials, couplingConstant, orthfile, numSpc
    int numclusters = atoi(argv[argbase++]);
    numclusters +=1; // add one more cluster
    
    int numTrials = atoi(argv[argbase++]);
    double couplingConstant = atof(argv[argbase++]);
    char *orthfile = argv[argbase++];
    int numSpc = atoi(argv[argbase++]);
    
    // create network object
    SpeciesNetwork **sns = new SpeciesNetwork*[numSpc];
    // read network files
    for (int i=0; i<numSpc; i++) {
        char *spc = argv[argbase++];
        sns[i] = ReadNetworkFile(spc);
    }
    // read parameters: startclusfn '-s' and maxtemp '-t'
    char startclusfn[1024]; startclusfn[0] = 0;
    double maxtemp = -1;
    while (argbase < argc) {
        if (!strcmp(argv[argbase], "-s")) {
            argbase++;
            strcpy(startclusfn, argv[argbase++]);
            continue;
        }
        if (!strcmp(argv[argbase], "-t")) {
            argbase++;
            maxtemp = atof(argv[argbase++]);
            continue;
        }
    }
    // parse orth file
    Orthology *orth = ReadOrthologyFile(orthfile, sns, numSpc);
    // create cluster object
    Cluster **c = new Cluster *[numTrials];
    
//    thread t_ary[NUM_THREADS]; // multi-threading
    
    for (int i=0; i<numTrials; i++) {
        fprintf(stderr,"launch from main\n");
        // initialize clustering parameters
        c[i] = new Cluster(sns, numSpc, orth, numclusters, couplingConstant);
        // use ground truth cluster label as start point
        if (startclusfn[0]!=0) {
            fprintf(stderr,"Using seed clustering file %s\n",startclusfn);
            c[i]->Preset(startclusfn);
        }
        // set maxtemp
        if (maxtemp >= 0) c[i]->SetMaxTemp(maxtemp);
        // simulated annealing to find parameters with best cost
        
//        t_ary[i] = thread(&Cluster::LearnGroundState,c[i]); // multi-threading

        c[i]->LearnGroundState();
    }
    // multi-threading
//    for (int i=0;i<NUM_THREADS;i++){
//        if (t_ary[i].joinable())
//            t_ary[i].join();
//    }
    
    double bestCost = 0.0;
    int bestTrial = -1;
    // go through all trails, and record the best cost
    for (int i=0; i<numTrials; i++) {
        // if (c[i]->CurrentCost < bestCost) { // the best cost in a trail is CurrentCost
        if (c[i]->BestCost < bestCost) { // the best cost in a trail is BestCost
            bestTrial = i;
            // bestCost = c[i]->CurrentCost;
            bestCost = c[i]->BestCost;
        }
    }
    
    printf("Best clustering, with cost %g is:\n", bestCost);
    c[bestTrial]->Print();
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    fprintf(stderr,"Total time: %f s\n",elapsed_secs);
    return 0;
}
