#ifndef COMIC_GRAPH_H
#define COMIC_GRAPH_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include "head.h"
using namespace std;

typedef double (*pf)(int,int);

class Graph
{
public:
    int n, m, k;
    vector<int> inDeg;
    vector<int> outDeg;
    vector<vector<int>> gT; // transpose graph, maintaining in-neighbors
    vector<vector<double>> probT;
	vector<vector<int>> gO; // non-transpose graph, maintaining out-neighbors
	vector<vector<double>> probO;

	vector<vector<int>>  outEdgeStatus;
	//unordered_map<int, int> edgeStatusMap;
	vector<unordered_map<int, int>> localEdgeStatus;

    string folder;
    string graph_file;
    vector<bool> hasNode;

    Graph(string folder, string graph_file):folder(folder), graph_file(graph_file) {
	    // read # nodes and # edges
	    readNM();
	    // init vectors
	    for (int i = 0; i < n; i++) {
		    hasNode.push_back(false);
		    gT.push_back(vector<int>());  // transpose graph
		    probT.push_back(vector<double>()); // incoming edge weight
		    inDeg.push_back(0); // in-degree of each node

		    gO.push_back(vector<int>()); // normal graph
		    probO.push_back(vector<double>()); // outgoing edge weight
		    outEdgeStatus.push_back(vector<int>()); 
		    outDeg.push_back(0); // out-degree of each node

		    localEdgeStatus.push_back(unordered_map<int, int>());
	    }
	    // read edge list from file
	    readGraph_v2();
    }

    ~Graph(){}

	void readNM();
    void add_in_edge(int a, int b, double p);  // add edge (a,b): a as b's in-neighbour
    void add_out_edge(int a, int b, double p); // add edge (a,b): b as a's out-neighbour
    void readGraph();
    void readGraph_v2();
	void printGraph(int num_nodes);
	vector<int> findHighDegree(int k);
	vector<int> findPageRank(int k);

	inline int cantor(int x, int y) {
		return ((x+y) * (x+y+1)) / 2 + y;
	}

	inline vector<int> inv_cantor(int z) {
		vector<int> ret;
		int w = (int) floor( (sqrt(8*z+1) - 1) / 2 );
		int t = (w * w + w) / 2;
		int y = z - t;
		int x = w - y;
		ret.push_back(x);
		ret.push_back(y);
		return ret;
	}

	inline void reset_out_edge_status() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < outDeg[i]; j++) {
				outEdgeStatus[i][j] = INACTIVE;
			}
		}
	}

/*
	inline void reset_edge_status_map() {
		for (auto kv : edgeStatusMap) {
			kv.second = INACTIVE;
		}
	}
*/

	// reset all edge status to INACTIVE
	inline void reset_local_edge_status() {
		for (int i = 0; i < n; i++) {
			unordered_map<int, int> mymap = localEdgeStatus[i];
			for (auto kv : mymap) {
				kv.second = INACTIVE;
			}
		}
	}

	// reset only edges coming from nodes that are in the input vector<int>
	inline void reset_local_edge_status(vector<int> touched) {
		for (auto i : touched) {
			unordered_map<int, int> mymap = localEdgeStatus[i];
			for (auto kv : mymap) {
				kv.second = INACTIVE;
			}
		}
	}

};

#endif