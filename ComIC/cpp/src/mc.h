#ifndef __MC_H
#define __MC_H

#include "graph.h"
#include <fstream>
#include <sstream>

class MonteCarlo // : public Graph
{
public:
	vector<int>     aSeeds; // seed set to be found (for item A)
	vector<int>     bSeeds; // input: the other company's seed set
	vector<double>  mg;  	// marginal gain of each seed added (in greedy sequence)
	double          qao;
	double          qab;
	double          qbo;
	double          qba;
	int             k;
	bool            ignore_B;
	string          output_file_name;
	string 			b_seeds_file_name;
	string			a_seeds_file_name;

	Graph *graph;
	int n;
	int m;

	MonteCarlo(string folder, string graph_file) 
	{
		ignore_B = false;
		k = 0;
		qao = qab = qbo = qba = 0;
		output_file_name = folder + "/output/default_output.txt";

		graph = new Graph(folder, graph_file);
		this->n = graph->n;
		this->m = graph->n;
	}

	~MonteCarlo()
	{
		delete graph;
	}

	void setParametersMC(int k_A, vector<double> qq, bool ignore, string folder, string bseeds);
	void readBSeedsMC();
	double mineSeedsMC(); 
	double compute_coverage(int *set_A, int size_A);
	/* if v becomes X-adopted (X = A or X = B), we examine v's out-neighbours */
	void examine_out_neighbors(int v, deque<int> *p_list, int *p_next, int *status); 
	

	void readASeedsMC_comp(string str);
	double mineSeedsMC_comp(double baseSpread);
	double compute_coverage_comp(vector<int> set_B, int size_B);

	void estimate_spread_from_file(string seed_file_name);
	void write_seeds_to_file(bool isB);
};


#endif
