#ifndef COMIC_TIMGRAPH_H
#define COMIC_TIMGRAPH_H

#include "infgraph.h"

class TimGraph: public InfGraph {
public:
	string output_file_name;

	TimGraph(string folder, string graph_file):InfGraph(folder, graph_file){}
	~TimGraph(){}

	double EstimateEPT();
	void BuildHyperGraph2(double epsilon, double ept);
	void BuildHyperGraph3(double epsilon, double opt);
	void mineSeeds_TIM_SIM(double epsilon);
	
	double EstimateEPT_fast();
	void BuildHyperGraph2_fast(double epsilon, double ept);
	void BuildHyperGraph3_fast(double epsilon, double opt);
	void mineSeeds_TIM_SIM_fast(double epsilon);
	void writeSeedsToFile_SIM(bool is_fast_verion, double epsilon);

	// for CIM
	void mineSeeds_CIM(double epsilon);
	double EstimateEPT_CIM();
	void BuildHyperGraph2_CIM(double epsilon, double ept);
	void BuildHyperGraph3_CIM(double epsilon, double opt);
	void writeSeedsToFile_CIM(double epsilon);


	inline double logcnk(int n, int k) {
		double ans=0;
		for(int i=n-k+1; i<=n; i++){
			ans+=log(i);
		}
		for(int i=1; i<=k; i++){
			ans-=log(i);
		}
		return ans;
	}

};

#endif //COMIC_TIMGRAPH_H
