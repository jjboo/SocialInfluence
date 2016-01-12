#ifndef COMIC_INFGRAPH_H
#define COMIC_INFGRAPH_H

#include "graph.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <deque>
#include <cstring>

using namespace std;

class InfGraph //: public Graph
{
public:
	vector<vector<int>> hyperG;
	vector<vector<int>> hyperGT;  // all RR sets
	int64 hyperId;
	vector<int> seedSet; // seed set to be mined (A-seeds for Self-Inf-Max and B-seeds for Comp-Inf-Max)
	vector<int> aSeeds;  // S_A from input files (for Comp-Inf-Max)
	vector<int> bSeeds;  // S_B from input files (for Self-Inf-Max)
	double qao;
	double qab;
	double qbo;
	double qba;
	int k;
	double *alpha_A;
	double *alpha_B;
	int *status_A;
	int *status_B;
	deque<int> Q;
	//bool *drawn; // TRUE if its alpha values are drawn in RR-set generation process
	bool *visited;  // TRUE if being visited in the backward BFS
	string b_seeds_file_name;
	string output_file_name;
	string dataset;
	bool ignore_B;
	bool *used; // used[v] is true if v has been added to RR-set
	bool *discovered;

	Graph *graph;
	int n, m;

	InfGraph(string folder, string graph_file) //:Graph(folder, graph_file) 
	{
		dataset = folder;

		graph = new Graph(folder, graph_file);
		this->n = graph->n;
		this->m = graph->m; 

		hyperId = 0;

		hyperG.clear();
		for(int i = 0; i < n; i++)
			hyperG.push_back(vector<int>());

		alpha_A = new double[n];
		alpha_B = new double[n];
		status_A = new int[n];
		status_B = new int[n];
		visited = new bool[n];
		used = new bool[n];
		discovered = new bool[n];
	}

	~InfGraph()
	{
		delete[] alpha_A;
		delete[] alpha_B;
		delete[] status_A;
		delete[] status_B;
		delete[] visited;
		delete[] used;
		delete[] discovered;
		delete graph;
	}

	void setParametersIG(int k, vector<double> qq, bool ignore, string folder, string bseeds);
	void BuildHypergraphR(int64);
	int BuildSingleRRSet(int uStart, int rrSetID, bool addHyperEdge); // same as the above function, but with diff. data structure
	void BuildSeedSetIG();
	double computeInfluenceHyperGraph();  // compute influence
	void readBSeedsIG();
	void printSeeds();

	// for RR-SIM-Fast
	void BuildHypergraphR_fast(int64);
	int BuildSingleRRSet_fast(int uStart, int rrSetID, bool addHyperEdge);

	// for RR-CIM
	void readASeedsIG_CIM(string str);
	void BuildHypergraphR_CIM(int64);
	int BuildSingleRRSet_CIM(int uStart, int rrSetID, bool addHyperEdge);

	void printQQ();
};

#endif
