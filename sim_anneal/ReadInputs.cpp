#include "cs-grn.h"
#include "ReadInputs.h"
//#include <cmath>

bool SpeciesNetwork::IsEdge(int i, int j)
{
    if (i >= numNodes || j >= numNodes) {
        return false;
    }
    if (std::abs(network[i][j]) > edgeWeightThreshold) return true;
    return false;
}

double SpeciesNetwork::GetEdgeWeight(int i, int j)
{
    if (i >= numNodes || j >= numNodes) {
        return false;
    }
    return network[i][j];
}

SpeciesNetwork *ReadNetworkFile(char *fn)
{
    FILE *F = fopen(fn, "r");
    if (F==NULL) return NULL;
    
    // create object to store network
    SpeciesNetwork *sn = new SpeciesNetwork;
    
    char line[linelen];
    // first line: species name
    fgets(line, linelen, F);
    sscanf(line, "%s", (char *)&sn->spcname);
    
    // second line: number of nodes
    fgets(line, linelen, F);
    sscanf(line, "%d", &sn->numNodes);
    
    // create space for network
    sn->network = new double*[sn->numNodes];
    for (int i=0; i<sn->numNodes; i++) {
        sn->network[i] = new double[sn->numNodes];
        for (int j=0; j<sn->numNodes; j++) {
            sn->network[i][j] = 0;
        }
    }
    
    // store all nodes of a species into a vector
    for (int i=0; i<sn->numNodes;i++){
        list<int> edges;
        sn->adjacencyList.push_back(edges);
    }
    
    
    // following lines: edges and weights
    int uniqueNodes = 0;
    while (fgets(line, linelen, F)) {
        char *id1 = new char[linelen];
        char *id2 = new char[linelen];
        double wt;
        // wt = 1
        sscanf(line, "%s\t%s\t%lf",id1, id2, &wt);
        string str_id1(id1);
        string str_id2(id2);
        if (sn->nodeName2Id.find(str_id1)==sn->nodeName2Id.end()) { // if id1 is not found
            fprintf(stderr,"Assigning id %d to node %s in species %s\n",uniqueNodes, str_id1.c_str(), fn);
            sn->nodeId2Name[uniqueNodes] = str_id1;
            sn->nodeName2Id[str_id1] = uniqueNodes++;
            if (uniqueNodes > sn->numNodes) {
                printf("Error: In file %s, read more unique node ids (%d) than the %d expected\n", fn, uniqueNodes, sn->numNodes);
                exit(1);
            }
        }
        if (sn->nodeName2Id.find(str_id2)==sn->nodeName2Id.end()) {
            fprintf(stderr,"Assigning id %d to node %s in species %s\n",uniqueNodes, str_id2.c_str(), fn);
            sn->nodeId2Name[uniqueNodes] = str_id2;
            sn->nodeName2Id[str_id2] = uniqueNodes++;
            if (uniqueNodes > sn->numNodes) {
                printf("Error: In file %s, read more unique node ids (%d) than the %d expected\n", fn, uniqueNodes, sn->numNodes);
                exit(1);
            }
        }
        // wt = 1
        sn->network[sn->nodeName2Id[str_id1]][sn->nodeName2Id[str_id2]] = wt;
        sn->adjacencyList[sn->nodeName2Id[str_id1]].push_back(sn->nodeName2Id[str_id2]);
    }
    fclose(F);
    
    return sn;
}

Orthology *ReadOrthologyFile(char *orthfn, SpeciesNetwork **sns, int numSpc)
{
    // create space for the table
    Orthology *o = new Orthology;
    o->numSpecies = numSpc;
    o->orthtable = new hash_map<NodePair *, int, NodePairHasher, eqNodePair> *[numSpc];
    // o->weighted_orth = new hash_map<int , double, NodePairHasher, eqNodePair> *[numSpc]; // hash table to store 1/d(#ortho genes from spe1 to spe2)
    o->weighted_orth = new hash_map<int, double> *[numSpc];
    for (int i=0; i<numSpc; i++) {
        o->orthtable[i] = new hash_map<NodePair *, int, NodePairHasher, eqNodePair> [numSpc];
        // o->weighted_orth[i] = new hash_map<int , double, NodePairHasher, eqNodePair> [numSpc];
        o->weighted_orth[i] = new hash_map<int, double> [numSpc];
    }
    
    // create numerical ids for species names
    o->spcId2Name = new char *[numSpc];
    for (int i=0; i<numSpc; i++) {
        o->spcId2Name[i] = new char[linelen];
        strcpy(o->spcId2Name[i], sns[i]->spcname);
        o->spcName2Id[(char *)sns[i]->spcname] = i;
    }
    
    // read the file and populate table
    FILE *F = fopen(orthfn, "r");
    if (F==NULL) return NULL;
    
    char line[linelen];
    while (fgets(line, linelen, F)) {
        char spc1[linelen];
        char spc2[linelen];
        char id1[linelen];
        char id2[linelen];
        sscanf(line, "%s %s %s %s",(char *)&spc1, (char *)&spc2, (char *)&id1, (char *)&id2);
        string str_spc1(spc1);
        string str_spc2(spc2);
        string str_id1(id1);
        string str_id2(id2);
        // convert spc ids to numerics
        if (o->spcName2Id.find(str_spc1) == o->spcName2Id.end()) {
            printf("Error: orthology file %s has species name that wasnt seen in network files\n", orthfn);
            exit(1);
        }
        int spc1id = o->spcName2Id[str_spc1];
        if (o->spcName2Id.find(str_spc2) == o->spcName2Id.end()) {
            printf("Error: orthology file %s has species name that wasnt seen in network files\n", orthfn);
            exit(1);
        }
        int spc2id = o->spcName2Id[str_spc2];
        
        // convert gene ids to numerics
        if (sns[spc1id]->nodeName2Id.find(str_id1) == sns[spc1id]->nodeName2Id.end()) { continue; }
        if (sns[spc2id]->nodeName2Id.find(str_id2) == sns[spc2id]->nodeName2Id.end()) { continue; }
        int gene1id = sns[spc1id]->nodeName2Id[str_id1];
        int gene2id = sns[spc2id]->nodeName2Id[str_id2];
        
        NodePair *np = new NodePair;
        np->id1 = gene1id;
        np->id2 = gene2id;
        o->orthtable[spc1id][spc2id][np] = 1;
        o->weighted_orth[spc1id][spc2id][gene1id] += 1.0;
        // fprintf(stderr, "gene id: %d,  orth effect: %f\n", gene1id,  o->weighted_orth[spc1id][spc2id][gene1id]);
    }
    fclose(F);
    
    for (int spc1=0; spc1 < numSpc; spc1++) {
        for (int spc2=spc1+1; spc2 < numSpc; spc2++) {
            for (hash_map<int, double>::iterator iter = o->weighted_orth[spc1][spc2].begin();
                 iter != o->weighted_orth[spc1][spc2].end(); ++iter) {
                double orth_effect = 1.0 / iter->second;
                // iter->second = orth_effect;
                o->weighted_orth[spc1][spc2][iter->first] = orth_effect;
                // fprintf(stderr, "orth effect %f\n", o->weighted_orth[spc1][spc2][iter->first]);
            }
        }
    }
    
    return o;
}


