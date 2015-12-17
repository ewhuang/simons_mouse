#ifndef __CLUSTER__
#define __CLUSTER__

#include "cs-grn.h"
#define maxLogSize 1000

struct LogEntry {
    int spc;
    int node;
    int oldstate;
};

struct Cluster {
//    Cluster(){};
    Cluster(SpeciesNetwork **specnws, int numS, Orthology *orthology, int numC, double cC);
    void LearnGroundState();
    
    SpeciesNetwork **sns;
    Orthology *orth;
    int numSpecies;
    int numClusters;
    double couplingConstant;
    double maxtemp;
    
    int ** snsClusterAssn;	// assignment of each node in each species to a cluster
    int ** bestClusterAssn;
    int * snsNumNodes; // number of nodes in each species network
    int * totalNumEdges;
    int ** degree;
    
//    std::map<std::pair<int, int>,std::vector<int>> clusterContainer; // a hash table to store nodes for a cluster in a species
    
    double CurrentCost;
    double CurrentCostClus;
    
    double BestCost;
    
    struct LogEntry UndoLog[maxLogSize];
    int UndoLogSize;
    
    double OrthCost(int **ClusterAssn);
    double DeltaCost(int **ClusterAssn);
    double Cost(int **ClusterAssn);
    int **Copy(int **ClusterAssn);
    void TotalNumEdges();
    void Delete(int **ClusterAssn);
    void CopyOver(int **ClusterAssn, int **oldClusterAssn);
    void Preset(char *filename);
    void Perturb(int **ClusterAssn, int numchanges);
    void UndoPerturb(int **ClusterAssn);
    void Print();
    void PrintBrief();
    void SetMaxTemp(double mt);
    
};


#endif
