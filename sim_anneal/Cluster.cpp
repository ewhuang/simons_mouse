#include "Cluster.h"
#include "time.h"

// mutex m;

Cluster::Cluster(SpeciesNetwork **specnws, int numS, Orthology *orthology, int numC, double cC) {
    sns = specnws;
    numSpecies = numS;
    orth = orthology;
    numClusters = numC;
    couplingConstant = cC;
    maxtemp = MAXTEMP;
    
    UndoLogSize = 0;
    snsNumNodes = new int[numSpecies];
    
    for (int i=0; i<numSpecies; i++) {
        snsNumNodes[i] = sns[i]->numNodes;
    }
    // total number of edges for each species
    totalNumEdges = new int[numSpecies];
    // degree for each node in each species
    degree = new int*[numSpecies];
    for (int i=0; i<numSpecies; i++) {
        degree[i] = new int[snsNumNodes[i]];
    }
    // count total number of edges and degree
    TotalNumEdges();
    
    
    snsClusterAssn = new int *[numSpecies];
    bestClusterAssn = new int *[numSpecies];
    srand(time(NULL));
    // srand(1); // fix random seed
    
    for (int i=0; i<numSpecies; i++) {
        snsClusterAssn[i] = new int[snsNumNodes[i]];
        bestClusterAssn[i] = new int[snsNumNodes[i]];
        for (int j=0; j<snsNumNodes[i]; j++) {
            snsClusterAssn[i][j] = rand() % numClusters;
        }
    }
    CurrentCostClus = Cost(snsClusterAssn);
    CurrentCost = CurrentCostClus + couplingConstant * OrthCost(snsClusterAssn);
    
    // initialize BestCost
    BestCost = 0.0;
}

void Cluster::Print() {
    for (int i=0; i<numSpecies; i++) {
        //        for (hash_map<const char*, int, hash<const char*>, eqstr>::iterator iter = sns[i]->nodeName2Id.begin(); iter != sns[i]->nodeName2Id.end(); ++iter) {
        for (std::unordered_map<string, int>::iterator iter = sns[i]->nodeName2Id.begin(); iter != sns[i]->nodeName2Id.end(); ++iter){
            //            const char *name = (const char *)(iter->first);
            string name = iter->first;
            int id = sns[i]->nodeName2Id[name];
            // printf("Species %d\tGene %s\tCluster %d\n",i,name,snsClusterAssn[i][id]);
            printf("Species %d\tGene %s\tCluster %d\n",i,name.c_str(),bestClusterAssn[i][id]);
        }
    }
}

void Cluster::PrintBrief() {
    for (int i=0; i<numSpecies; i++) {
        printf("%d",snsClusterAssn[i][0]);
        for (int j=1; j<snsNumNodes[i]; j++) {
            printf("\t%d",snsClusterAssn[i][j]);
        }
        printf("\n");
    }
}

void Cluster::SetMaxTemp(double mt) {
    maxtemp = mt;
}

void Cluster::Preset(char *filename) {
    // m.lock();
    FILE *F = fopen(filename, "r");
    if (F==NULL) return;
    
    char line[linelen];
    fgets(line, linelen, F); // skip the single header line
    while (fgets(line, linelen, F)) {
        char spc[linelen];
        char id[linelen];
        char clusterid[linelen];
        char d1[linelen], d2[linelen], d3[linelen];
        // format of line is "Species <spc> Gene <gene> Cluster <cluster>"
        sscanf(line, "%s %s %s %s %s %s",(char *)&d1, (char *)&spc, (char *)&d2, (char *)&id, (char *)&d3, (char *)&clusterid);
        
        // convert spc id to numerics
        if (orth->spcName2Id.find(spc) == orth->spcName2Id.end()) {
            printf("Error: seed-clustering file %s has species name that wasnt seen in network files\n", filename);
            exit(1);
        }
        int spcid = orth->spcName2Id[spc];
        
        // convert gene id to numerics
        if (sns[spcid]->nodeName2Id.find(id) == sns[spcid]->nodeName2Id.end()) {
            continue; // this may happen because clusters file includes information about genes that are not there in
            // ... the network file (gene node with degree 0). It doesnt matter where (which cluster) we place such genes.
        }
        int geneid = sns[spcid]->nodeName2Id[id];
        
        // record information
        snsClusterAssn[spcid][geneid] = atoi(clusterid);
        fprintf(stderr,"Initial clustering has species %d gene %s assigned to cluster %d\n", spcid, id, snsClusterAssn[spcid][geneid]);
    }
    fclose(F);
    CurrentCostClus = Cost(snsClusterAssn);
    CurrentCost = CurrentCostClus + couplingConstant * OrthCost(snsClusterAssn);
    // m.unlock();
}

void Cluster::TotalNumEdges(){
    for (int spc = 0; spc < numSpecies; spc++) {
        // compute degree of each node
        // int *degree = new int[snsNumNodes[spc]];
        // initialize degree value
        for (int i=0; i<snsNumNodes[spc]; i++) degree[spc][i] = 0;
        
        for (int i=0; i<snsNumNodes[spc]; i++) {
            for (int j=i+1; j<snsNumNodes[spc]; j++) {
                // if there is an edge between i and j
                if (sns[spc]->IsEdge(i,j)) {
                    degree[spc][i]++;
                    degree[spc][j]++;
                }
            }
        }
        totalNumEdges[spc] = 0;
        for (int i=0; i<snsNumNodes[spc]; i++) {
            totalNumEdges[spc] += degree[spc][i];
        }
        // total number of edges, degree of both nodes added one each time
        totalNumEdges[spc] /= 2;
    }
}

void Cluster::LearnGroundState() {
//    m.lock();
    fprintf(stderr,"Called LearnGroundState\n");
    double tempmax = maxtemp;
    double tempmin = tempmax/1000;
    
    int totalNumNodes = 0;
    for (int i=0; i<numSpecies; i++) totalNumNodes += snsNumNodes[i];
    
    int nochangeiter = 0; // how many iterations have we seen no change in score
    for(double temp=tempmax; temp >= tempmin; temp *= 0.9) { // modify this
       if (nochangeiter > 500) break;
        // nochangeiter = 0;
        // a "sweep" of the network does about 50 x as many changes as there are nodes overall
        for (int i=0; i<50*totalNumNodes; i++) {
            if (nochangeiter > 500) break;
            fprintf(stderr, "C = %g\n", CurrentCost); // prev log
            if (CurrentCost < BestCost){
                BestCost = CurrentCost;
                for (int i=0; i<numSpecies; i++){
                    for (int j=0; j<snsNumNodes[i];j++){
                        bestClusterAssn[i][j] = snsClusterAssn[i][j];
                    }
                }
                // fprintf(stderr, "best cost: %f\n", BestCost);
            }
            // printf("%g\n",CurrentCost);
            // PrintBrief();
            
            double OldCost = CurrentCost;
            double OldCostClus = CurrentCostClus;
            // int **oldClusterAssn = Copy(snsClusterAssn);
            
            // propose new assignment, by randomly perturbing the current assignment
            // Perturb(snsClusterAssn, int(totalNumNodes/100));
            // large number of perturb gives worse resutls
            Perturb(snsClusterAssn, 1);
            double deltaCostCluster = DeltaCost(snsClusterAssn);
            
            // score new assignment
            // change this to dealtaCost
            double NewCostClus = OldCostClus + deltaCostCluster;
            double NewCost = NewCostClus + couplingConstant * OrthCost(snsClusterAssn);
            
            double deltaCost = NewCost - OldCost;
            
//            fprintf(stderr, "Delta Cost: %g\n", deltaCostCluster);
            
            // bad move
            if (deltaCost > 0) {
                // reject change with pr = 1 - exp(-delta/T)
                if (double(rand())/RAND_MAX >= exp(-deltaCost/temp)) {
                    // if new cost is more, i.e., worse,
                    // delta is more positive, thus exp(-delta/T) is closer to 0,
                    // thus reject probability is larger
                    
                    // fprintf(stderr, "REJECTED MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, 1-exp(-deltaCost/temp)); // prev log
                    UndoPerturb(snsClusterAssn);
                    // CopyOver(snsClusterAssn, oldClusterAssn);
                    nochangeiter++;
                }
                else { // accept bad move
                    CurrentCost = NewCost;
                    CurrentCostClus = NewCostClus;
                    // fprintf(stderr, "BAD MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, exp(-deltaCost/temp)); // prev log
                    nochangeiter = 0;
                    UndoLogSize = 0;
                }
            }
            else { // always accept good moves
                CurrentCost = NewCost;
                CurrentCostClus = NewCostClus;
                // fprintf(stderr, "GOOD MOVE\t%g to %g at t %g\n", OldCost, NewCost, temp); // prev log
                if (deltaCost > 1) { // any improvement less than this is not counted as an improvement
                    nochangeiter = 0;
                }
                UndoLogSize = 0;
                // update BestCost
                if (CurrentCost + 1 < BestCost){
                    BestCost = CurrentCost;
                    bestClusterAssn = snsClusterAssn;
                    for (int i=0; i<numSpecies; i++){
                        for (int j=0; j<snsNumNodes[i];j++){
                            bestClusterAssn[i][j] = snsClusterAssn[i][j];
                        }
                    }
                }
                
            }
//            if (nochangeiter>500) break;
        }
    }
//    m.unlock();
}

double Cluster::DeltaCost(int **ClusterAssn) {
    
    double cost = 0.0;
    
    int node = UndoLog[UndoLogSize-1].node;
    int spc = UndoLog[UndoLogSize-1].spc;
    int oldstate = UndoLog[UndoLogSize-1].oldstate;
    int curstate = ClusterAssn[spc][node];
    
    std::vector<int> new_cluster_node;
    std::vector<int> old_cluster_node;
    
    int new_cluster_size = 0;
    int old_cluster_size = 0;
    
    for (int i=0; i<snsNumNodes[spc]; i++) {
        
        if (ClusterAssn[spc][i] == curstate){
            new_cluster_size ++;
            new_cluster_node.push_back(i);
        }
        if (ClusterAssn[spc][i] == oldstate){
            old_cluster_size ++;
            old_cluster_node.push_back(i);
        }
    }
    old_cluster_node.push_back(node);
    old_cluster_size ++;
//    fprintf(stderr, "oldClusterSize %d, newClusterSize %d\n", old_cluster_size, new_cluster_size);
    /********* new cluster ************/
    int in_count_before_perturb = 0;
    int in_count_after_perturb = 0;
    int in_edge_before_perturb = 0;
    int in_edge_after_perturb = 0;
    int out_count_before_perturb = 0;
    int out_count_after_perturb = 0;
    int out_edge_before_perturb = 0;
    int out_edge_after_perturb = 0;
    
    for (int i = 0; i < new_cluster_node.size(); i++){
        // orig new cluster without added node
        if (new_cluster_node[i] != node){
            in_count_before_perturb ++;
            out_count_before_perturb += (snsNumNodes[spc] - new_cluster_size + 1);
            out_edge_before_perturb += degree[spc][new_cluster_node[i]];
        }
        // new cluster with added node
        in_count_after_perturb ++;
        out_count_after_perturb += (snsNumNodes[spc] - new_cluster_size);
        out_edge_after_perturb += degree[spc][new_cluster_node[i]];
        
        for (list<int>::iterator it = sns[spc]->adjacencyList[new_cluster_node[i]].begin(); it!=sns[spc]->adjacencyList[new_cluster_node[i]].end(); it++){
            if (ClusterAssn[spc][*it] == curstate){
                if (*it != node and new_cluster_node[i] != node){
//                    fprintf(stderr, "spe %d, node1 %d, node2 %d\n",spc, new_cluster_node[i],*it);
                    in_edge_before_perturb ++;
                }
                in_edge_after_perturb ++;
            }
        }
        
//        for (int j = 0; j < new_cluster_node.size(); j++)
//        {
//            // only count old nodes
//            if (new_cluster_node[i] != node && new_cluster_node[j] != node)
//            {
//                in_count_before_perturb ++;
//                if (sns[spc]->IsEdge(new_cluster_node[i],new_cluster_node[j]))
//                {
//                    in_edge_before_perturb ++;
//                }
//            }
//            // count new + old nodes
//            in_count_after_perturb ++;
//            if (sns[spc]->IsEdge(new_cluster_node[i],new_cluster_node[j]))
//            {
//                in_edge_after_perturb ++;
//            }
//        }
        
    }
//    fprintf(stderr, "out_degree_before_perturb %d, in_edge_before_perturb %d\n", out_edge_before_perturb, in_edge_before_perturb);
    out_edge_before_perturb = (out_edge_before_perturb) - in_edge_before_perturb;
    
    out_edge_after_perturb = (out_edge_after_perturb) - in_edge_after_perturb;
    out_count_before_perturb = out_count_before_perturb;
    out_count_after_perturb = out_count_after_perturb;
    
//     fprintf(stderr, "out_edge_before_perturb: %d\t out_edge_after_perturb: %d\t out_count_before_perturb: %d\t out_count_after_perturb: %d\t\n", out_edge_before_perturb, out_edge_after_perturb, out_count_before_perturb, out_count_after_perturb);
    
    // fprintf(stderr, "in_edge_before_perturb: %d\t in_edge_after_perturb: %d\t in_count_before_perturb: %d\t in_count_after_perturb: %d\t\n", in_edge_before_perturb, in_edge_after_perturb, in_count_before_perturb, in_count_after_perturb);
    double in_density_before_perturb = 0.0;
    if (in_count_before_perturb >= 1 and curstate != 0){ // new_cluster is not 0
        int in_count_before_perturb_sqr = in_count_before_perturb * in_count_before_perturb;
        in_density_before_perturb = double(in_edge_before_perturb) / (in_count_before_perturb * in_count_before_perturb);
        double out_density_before_perturb = double(out_edge_before_perturb) / out_count_before_perturb;
        double inOutRatio_before_perturb = (in_density_before_perturb / (in_density_before_perturb + out_density_before_perturb)) * double(new_cluster_size - 1);
        cost = cost - inOutRatio_before_perturb;
//        fprintf(stderr, "in_count_before_perturb %d, in_density_before_perturb %g, out_density_before_perturb %g, inOutRatio_before_perturb %g\n", in_count_before_perturb_sqr, in_density_before_perturb, out_density_before_perturb, inOutRatio_before_perturb);
    }
    
//     fprintf(stderr, "in_density_before_perturb: %g\t out_density_before_perturb: %g\t inOutRatio_before_perturb: %g\t new_cluster_size: %d\t\n", in_density_before_perturb, out_density_before_perturb, inOutRatio_before_perturb, new_cluster_size);
    
    
    double in_density_after_perturb = 0.0;
    if (in_count_after_perturb >=1 and curstate != 0){
        int in_count_after_perturb_sqr = in_count_after_perturb * in_count_after_perturb;
        in_density_after_perturb = double(in_edge_after_perturb) / (in_count_after_perturb * in_count_after_perturb);
        double out_density_after_perturb = double(out_edge_after_perturb) / out_count_after_perturb;
        double inOutRatio_after_perturb = (in_density_after_perturb / (in_density_after_perturb + out_density_after_perturb)) * double(new_cluster_size);
        cost = cost + inOutRatio_after_perturb;
//        fprintf(stderr, "in_count_after_perturb %d, in_density_after_perturb %g, out_density_after_perturb %g, inOutRatio_after_perturb %g\n", in_count_after_perturb_sqr, in_density_after_perturb, out_density_after_perturb, inOutRatio_after_perturb);
    }
    
    
    /********* old cluster ************/
    in_count_before_perturb = 0;
    in_count_after_perturb = 0;
    in_edge_before_perturb = 0;
    in_edge_after_perturb = 0;
    out_count_before_perturb = 0;
    out_count_after_perturb = 0;
    out_edge_before_perturb = 0;
    out_edge_after_perturb = 0;
    
    for (int i = 0; i < old_cluster_node.size(); i++){
        // orig old cluster without removed node
        if (old_cluster_node[i] != node){ // exclude new node
            in_count_after_perturb ++;
            out_count_after_perturb += (snsNumNodes[spc] - old_cluster_size + 1);
            out_edge_after_perturb += degree[spc][old_cluster_node[i]];
        }
        // old cluster with removed node
        in_count_before_perturb ++;
        out_count_before_perturb += (snsNumNodes[spc] - old_cluster_size);
        out_edge_before_perturb += degree[spc][old_cluster_node[i]];
        
        for (list<int>::iterator it = sns[spc]->adjacencyList[old_cluster_node[i]].begin(); it!=sns[spc]->adjacencyList[old_cluster_node[i]].end(); it++){
//            fprintf(stderr, "node1 %d, node2 %d\n", old_cluster_node[i], *it);
            if (ClusterAssn[spc][*it] == oldstate || *it == node){
//                fprintf(stderr, "oldstate %d, node1 %d, node2 %d\n", oldstate, old_cluster_node[i], *it);
                if (*it != node and old_cluster_node[i] != node){
                    in_edge_after_perturb ++;
                }
                in_edge_before_perturb ++;
            }
        }
        
//        for (int j = 0; j < old_cluster_node.size(); j++)
//        {
//            // only count old nodes
//            if (old_cluster_node[i] != node && old_cluster_node[j] != node)
//            {
//                in_count_after_perturb ++;
//                if (sns[spc]->IsEdge(old_cluster_node[i],old_cluster_node[j]))
//                {
//                    in_edge_after_perturb ++;
//                }
//            }
//            // count new + old nodes
//            in_count_before_perturb ++;
//            if (sns[spc]->IsEdge(old_cluster_node[i],old_cluster_node[j]))
//            {
//                in_edge_before_perturb ++;
//            }
//        }
        
    }
    
    out_edge_before_perturb = (out_edge_before_perturb) - in_edge_before_perturb;
    out_edge_after_perturb = (out_edge_after_perturb) - in_edge_after_perturb;
    out_count_before_perturb = out_count_before_perturb;
    out_count_after_perturb = out_count_after_perturb;
//    fprintf(stderr, "out_edge_before_perturb: %d\t out_edge_after_perturb: %d\t out_count_before_perturb: %d\t out_count_after_perturb: %d\t\n", out_edge_before_perturb, out_edge_after_perturb, out_count_before_perturb, out_count_after_perturb);
    
    in_density_before_perturb = 0.0;
    if (in_count_before_perturb >= 1 and oldstate != 0){
        int in_count_before_perturb_sqr = in_count_before_perturb * in_count_before_perturb;
        in_density_before_perturb = double(in_edge_before_perturb) / (in_count_before_perturb * in_count_before_perturb);
        
        // in_density_before_perturb = double(in_edge_before_perturb) / in_count_before_perturb;
        
        double out_density_before_perturb = double(out_edge_before_perturb) / out_count_before_perturb;
        double inOutRatio_before_perturb = (in_density_before_perturb / (in_density_before_perturb + out_density_before_perturb)) * double(old_cluster_size);
        cost = cost - inOutRatio_before_perturb;
//        fprintf(stderr, "in_count_before_perturb %d, in_density_before_perturb %g, out_density_before_perturb %g, inOutRatio_before_perturb %g, in_edge_before_perturb %d\n", in_count_before_perturb_sqr, in_density_before_perturb, out_density_before_perturb, inOutRatio_before_perturb, in_edge_before_perturb);
    }
    
    in_density_after_perturb = 0.0;
    if (in_count_after_perturb >=1 and oldstate != 0){
        int in_count_after_perturb_sqr = in_count_after_perturb * in_count_after_perturb;
        in_density_after_perturb = double(in_edge_after_perturb) / (in_count_after_perturb * in_count_after_perturb);
        
        // in_density_after_perturb = double(in_edge_after_perturb) / in_count_after_perturb;
        double out_density_after_perturb = double(out_edge_after_perturb) / out_count_after_perturb;
        double inOutRatio_after_perturb = (in_density_after_perturb / (in_density_after_perturb + out_density_after_perturb)) * double(old_cluster_size-1);
        cost = cost + inOutRatio_after_perturb;
//        fprintf(stderr, "in_count_after_perturb %d, in_density_after_perturb %g, out_density_after_perturb %g, inOutRatio_after_perturb %g\n", in_count_after_perturb_sqr, in_density_after_perturb, out_density_after_perturb, inOutRatio_after_perturb);
    }
    // fprintf(stderr, "Delta Cost %g\n", -cost);
    return -cost;
}

double Cluster::Cost(int **ClusterAssn) {
    double cost = 0.0;
    // for each species
    for (int spc = 0; spc < numSpecies; spc++)
    {
        std::unordered_map<int, int> total_cluster_outEdges;
        std::unordered_map<int, int> total_cluster_edge;
        std::unordered_map<int, int> total_cluster_inCount;
        std::unordered_map<int, int> total_cluster_outCount;
        std::unordered_map<int, int> total_cluster_nodes;
        // for node i in a species
        for (int i=0; i<snsNumNodes[spc]; i++)
        {
            // initialize in-cluster total edge = 0
            total_cluster_edge[ClusterAssn[spc][i]] += 0;
            total_cluster_nodes[ClusterAssn[spc][i]] ++;
            //            fprintf(stderr, "cluster nodes: %d \n", total_cluster_nodes[ClusterAssn[spc][i]]);
            for (int j=0; j<snsNumNodes[spc]; j++)
            { // should this be j=i+1 ...?
                // if (i == j)
                // continue;
                // i and j in the same cluster
                if (ClusterAssn[spc][i] == ClusterAssn[spc][j])
                {
                    total_cluster_inCount[ClusterAssn[spc][i]] ++;
                    
                    // i and j have an edge
                    if (sns[spc]->IsEdge(i,j)){
                        // sum in-cluster edges
                        total_cluster_edge[ClusterAssn[spc][i]] ++;
                        // fprintf(stderr, "cluster edges: %d \n", total_cluster_edge[ClusterAssn[spc][i]]);
                    }
                }
                else{
                    total_cluster_outCount[ClusterAssn[spc][i]] ++;
                    if (sns[spc]->IsEdge(i,j))
                        total_cluster_outEdges[ClusterAssn[spc][i]] ++;
                    
                }
            }
        }
        for (std::unordered_map<int, int>::iterator it=total_cluster_nodes.begin(); it != total_cluster_nodes.end(); it++)
        {
//            fprintf(stderr, "cluster %d cluster size\n", it->first, it->second);
            if (it->second >= 1 and it->first != 0 ){ // skip '0' cluster
                double total_cluster_inDensity = double(total_cluster_edge[it->first]) / double(total_cluster_inCount[it->first]);
                // double total_cluster_outDensity = double(total_cluster_outEdges[it->first]) / (snsNumNodes[spc] - total_cluster_nodes[it->first]);
                double total_cluster_outDensity = double(total_cluster_outEdges[it->first]) / double(total_cluster_outCount[it->first]);
                
                cost += (total_cluster_inDensity / (total_cluster_inDensity + total_cluster_outDensity)) * double(it->second);
                //                 fprintf(stderr, "cluster size %d in-density %g out-density %g cost %g \n", it->second, total_cluster_inDensity, total_cluster_outDensity, cost);
                // cost += double(cluster_modularity[it->first]) / total_cluster_nodes[it->first];
            }
            else{
                cost += 0.0;
            }
        }
    }
    return - cost;
}

double Cluster::OrthCost(int **ClusterAssn){
    double orthterms = 0.0;
    for (int spc1=0; spc1 < numSpecies; spc1++) {
        for (int spc2=spc1+1; spc2 < numSpecies; spc2++) {
            for (hash_map<NodePair *, int, NodePairHasher, eqNodePair>::iterator iter = orth->orthtable[spc1][spc2].begin();
                 iter != orth->orthtable[spc1][spc2].end(); ++iter) {
                NodePair *np = (NodePair *)(iter->first);
                int n1 = np->id1;
                int n2 = np->id2;
                // fprintf(stderr, "orth\t%d\t%d\t%d\t%d\t%d\t%d\n", spc1, spc2, n1, n2, ClusterAssn[spc1][n1], ClusterAssn[spc2][n2]);
                // fprintf(stderr, "orth effect %f\n", orth->weighted_orth[spc1][spc2][n1]);
                // if (ClusterAssn[spc1][n1] == ClusterAssn[spc2][n2]) orthterms = orthterms+1;
                if (ClusterAssn[spc1][n1] == ClusterAssn[spc2][n2] and ClusterAssn[spc1][n1] != 0)
                    orthterms += (orth->weighted_orth[spc1][spc2][n1] + orth->weighted_orth[spc1][spc2][n2]) / 2.0;
                // fprintf(stderr, "sum of orth effect %f\n", orthterms);7
            }
        }
    }
//        fprintf(stderr, "cc = %g\tcccost = %g\n",couplingConstant, couplingConstant*orthterms); // prev log
    return -orthterms;
}

int **Cluster::Copy(int **ClusterAssn) {
    
    int **ca = new int *[numSpecies];
    for (int i=0; i<numSpecies; i++) {
        ca[i] = new int[snsNumNodes[i]];
        for (int j=0; j<snsNumNodes[i]; j++) ca[i][j] = ClusterAssn[i][j];
    }
    return ca;
}

void Cluster::Delete(int **ClusterAssn) {
    
    for (int i=0; i<numSpecies; i++) {
        delete [] ClusterAssn[i];
    }
    delete [] ClusterAssn;
}

void Cluster::CopyOver(int **ClusterAssn, int **oldClusterAssn) {
    
    for (int i=0; i<numSpecies; i++) {
        for (int j=0; j<snsNumNodes[i]; j++) ClusterAssn[i][j] = oldClusterAssn[i][j];
    }
}

void Cluster::UndoPerturb(int **ClusterAssn) {
    
    for (int i = UndoLogSize; i>0; i--) {
        ClusterAssn[UndoLog[i-1].spc][UndoLog[i-1].node] = UndoLog[i-1].oldstate;
    }
    UndoLogSize = 0;
}

void Cluster::Perturb(int **ClusterAssn, int numchanges = 1) {
    
    double deltaCostAtPerturb = 0.0;
    for (int i=0; i<numchanges; i++) {
        // choose a random species
        int spc = rand() % numSpecies;
        // chose a random node
        int node = rand() % (snsNumNodes[spc]);
        while (1) {
            int oldstate = ClusterAssn[spc][node];
            // choose a randome new cluster
            int newstate = rand() % numClusters;
            
            if (newstate == oldstate) continue;
//             fprintf(stderr,"Per\t%d\t%s\t%d\t%d\n",spc,sns[spc]->nodeId2Name[node].c_str(),oldstate,newstate); // prev log
            // assign node to new cluster
            ClusterAssn[spc][node] = newstate;
            if (UndoLogSize > maxLogSize) { // maxLogSize = 1000
                printf("Error: Undo Log Size Limit reached\n");
                exit(1);
            }
            // record the change for later undo
            UndoLog[UndoLogSize].spc = spc;
            UndoLog[UndoLogSize].node = node;
            UndoLog[UndoLogSize].oldstate = oldstate;
//             fprintf(stderr, "Pertub - species: %d, node: %d, old cluster: %d, new cluster: %d\n", spc, node, oldstate, newstate);
            UndoLogSize++;
            // the reason to move deltaCost here is because it captures state of a node before and right after update
            // deltaCostAtPerturb += DeltaCost(snsClusterAssn);
            break;
        }
    }
}
