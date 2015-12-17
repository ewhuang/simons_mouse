#ifndef __CSGRN__
#define __CSGRN__

#define linelen 1024
#define MAXTEMP 10
#define NUM_THREADS 2

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <string>
#include <ext/hash_map>
#include <ctime>
#include <vector>
#include <map>
#include <unordered_map>
#include <thread>
#include <cstring>
#include <list>
// #include <algorithm>

using namespace std;
using namespace __gnu_cxx;

double abs(double x);


struct eqstr
{
    bool operator()(const char* s1, const char* s2) const
    {
        
        return strcmp(s1, s2) == 0;
    }
};

struct eqint
{
    bool operator()(int i1, int i2) const
    {
        if(i1==i2)
            return true ;
        else
            return false ;
    }
};

struct SpeciesNetwork {
    int numNodes;
    char spcname[linelen];
    double **network;
    // adjacency list to store edges for each node
    std::vector<std::list<int>> adjacencyList;
    
    //    hash_map<const char*, int,hash<const char*>, eqstr> nodeName2Id;
    std::unordered_map<std::string, int> nodeName2Id;
    //    hash_map<int, const char*, hash<int>, eqint> nodeId2Name;
    std::unordered_map<int, std::string> nodeId2Name;
    
    bool IsEdge(int i, int j);
    double GetEdgeWeight(int i, int j);
    static const float edgeWeightThreshold = 0.1; // TO DO: make this configurable
};

struct NodePair {
    int id1;
    int id2;
};

struct eqNodePair
{
    bool operator()(NodePair* n1, NodePair* n2) const
    {
        return (n1->id1 == n2->id1 && n1->id2 == n2->id2);
    }
};

class NodePairHasher {
public:
    size_t operator()(const NodePair *r) const
    {
        char line[linelen];
        sprintf(line,"%d,%d",r->id1, r->id2);
        return h(line);
    };
    
private:
    __gnu_cxx::hash<char*> h;
};

struct Orthology {
    int numSpecies;
    hash_map<NodePair *, int, NodePairHasher, eqNodePair> **orthtable;
    //    hash_map<const char *, int, hash<const char *>, eqstr> spcName2Id;
    std::unordered_map<std::string, int> spcName2Id;
    // hash_map<int , double, NodePairHasher, eqNodePair> **weighted_orth;
    hash_map<int, double> **weighted_orth;
    char **spcId2Name;
};



#endif
