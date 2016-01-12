#include "timgraph.h"


double TimGraph::EstimateEPT()
{
	double lb = 0.5, c = 0;
	//int64 lastR = 0;
	int loopEstimateEPT = 0;
	while(true) {
		int loop = (6 * log(n)  +  6 * log(log(n)/ log(2)) ) * 1 / lb;
		c = 0;
		//lastR = loop;
		IF_TRACE(int64 now = rdtsc());
		double sumMgTu = 0;
		for(int i = 0; i < loop; i++)
		{
			int u = rand() % n; // randomly sample a target node
			double MgTu = (double) BuildSingleRRSet(u, 0, false); // generate an RR-set rooted at node u
			double pu = MgTu/m;
			sumMgTu += MgTu;
			c += 1 - pow((1-pu), k);
			loopEstimateEPT++;
		}
		c /= loop;
		if(c > lb)
			break;
		lb /= 2;
	}
	return 0.5 * c * n;
}

void TimGraph::BuildHyperGraph2(double epsilon, double ept)
{
	ASSERT(ept > 0);
	//int64 R =  (8+2 * epsilon) * ( n * log(n) +  n * log(2)  ) / ( epsilon * epsilon * ept) ;
	int64 R = (2 + epsilon) * ( n * log(n) ) / ( epsilon * epsilon * ept);
	cout << "[info] # of RR sets for computing KPT = " << R << endl;
	BuildHypergraphR(R);
}


void TimGraph::BuildHyperGraph3(double epsilon, double opt){
	ASSERT(opt > 0);
	int64 R = (8+2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, k) ) / ( epsilon * epsilon * opt);
	//int64 R = (8+2 * epsilon) * ( log(n) + log(2) +  n * logcnk(n, k) ) / ( epsilon * epsilon * opt);
	cout << "[info] # of RR sets for mining seeds = " << R << endl;
	BuildHypergraphR(R);
}

void TimGraph::mineSeeds_TIM_SIM(double epsilon)
{
	double eps_prime = 5 * pow( pow(epsilon, 2) / (k), 1.0 / 3.0 );

	double kpt_star = EstimateEPT();
	BuildHyperGraph2(eps_prime, kpt_star);
	BuildSeedSetIG();
	double kpt = computeInfluenceHyperGraph();
	kpt /= (1 + eps_prime);
	double kpt_plus = max(kpt, kpt_star);

	BuildHyperGraph3(epsilon, kpt_plus);
	BuildSeedSetIG();
	double spread = computeInfluenceHyperGraph();

	printSeeds();
	cout << "[info] final spread by RR-set model = " << spread << endl;
}



/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


void TimGraph::BuildHyperGraph2_fast(double epsilon, double ept)
{
	ASSERT(ept > 0);
	//int64 R = (8+2 * epsilon) * ( n * log(n) +  n * log(2)  ) / ( epsilon * epsilon * ept);
	int64 R = (2 + epsilon) * ( n * log(n) ) / ( epsilon * epsilon * ept);
	cout << "[info] # of RR sets for computing KPT = " << R << endl;
	BuildHypergraphR_fast(R);
}

void TimGraph::BuildHyperGraph3_fast(double epsilon, double opt)
{
	ASSERT(opt > 0);
	int64 R = (8+2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, k) ) / ( epsilon * epsilon * opt);
	cout << "[info] # of RR sets for mining seeds = " << R << endl;
	BuildHypergraphR_fast(R);
}

double TimGraph::EstimateEPT_fast()
{
	double lb = 0.5;
	double c = 0;
	//int64 lastR = 0;
	int loopEstimateEPT = 0;
	while(true) {
		int loop = (6 * log(n)  +  6 * log(log(n)/ log(2)) ) * 1 / lb;
		c = 0;
		//lastR = loop;
		IF_TRACE(int64 now = rdtsc());
		double sumMgTu = 0;
		for(int i = 0; i < loop; i++)
		{
			int u = rand() % n; // randomly sample a target node
			double MgTu = (double) BuildSingleRRSet_fast(u, 0, false); // generate an RR-set rooted at node u
			double pu = MgTu/m;
			sumMgTu += MgTu;
			c += 1 - pow((1-pu), k);
			loopEstimateEPT++;
		}
		c /= loop;
		if(c > lb)
			break;
		lb /= 2;
	}
	return 0.5 * c * n;
}


void TimGraph::mineSeeds_TIM_SIM_fast(double epsilon)
{
	double eps_prime = 5 * pow( pow(epsilon, 2) / (k), 1.0 / 3.0 );

	double kpt_star = EstimateEPT_fast();
	BuildHyperGraph2_fast(eps_prime, kpt_star);
	BuildSeedSetIG();
	double kpt = computeInfluenceHyperGraph();
	kpt /= (1 + eps_prime);
	double kpt_plus = max(kpt, kpt_star);

	BuildHyperGraph3_fast(epsilon, kpt_plus);
	BuildSeedSetIG();
	double spread = computeInfluenceHyperGraph();

	printSeeds();
	cout << "[info] final spread by RR-set model = " << spread << endl;
}



/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/

/*
void TimGraph::mineSeeds_CIM(double epsilon)
{
	double ep_step2, ep_step3;
	ep_step3 = epsilon;
	ep_step2 = 5 * pow( pow(ep_step3, 2) / k, 1.0 / 3.0 );

	double ept = EstimateEPT_CIM();
	BuildSeedSetIG();
	BuildHyperGraph2_CIM(ep_step2, ept);
	ept = computeInfluenceHyperGraph();
	ept /= 1 + ep_step2;

	cout << "[info]: KPT used = " << ept << endl;

	BuildHyperGraph3_CIM(ep_step3, ept);
	BuildSeedSetIG();
	ept = computeInfluenceHyperGraph();

	printSeeds();
	cout << "[info] final spread by RR-set model = " << ept << endl;
}
*/


void TimGraph::mineSeeds_CIM(double epsilon)
{
	double eps_prime = 5 * pow( pow(epsilon, 2) / (k), 1.0 / 3.0 );

	double kpt_star = EstimateEPT_CIM();
	BuildHyperGraph2_CIM(eps_prime, kpt_star);
	BuildSeedSetIG();
	double kpt = computeInfluenceHyperGraph();
	kpt /= (1 + eps_prime);
	double kpt_plus = max(kpt, kpt_star);

	cout << "[info]: KPT used = " << kpt_plus << endl;

	BuildHyperGraph3_CIM(epsilon, kpt_plus);
	BuildSeedSetIG();
	double spread = computeInfluenceHyperGraph();

	printSeeds();
	cout << "[info] final spread by RR-set model = " << spread << endl;
}


double TimGraph::EstimateEPT_CIM()
{
	double lb = 0.5, c = 0;
	int64 lastR = 0;
	int loopEstimateEPT = 0;
	while (true) {
		int loop = (6 * log(n)  +  6 * log(log(n)/ log(2)) ) * 1 / lb;
		c = 0;
		lastR = loop;
		IF_TRACE(int64 now = rdtsc());
		double sumMgTu = 0;
		for(int i = 0; i < loop; i++)
		{
			int u = rand() % n; // randomly sample a target node
			double MgTu = (double) BuildSingleRRSet_CIM(u, 0, false); // generate an RR-set rooted at node u
			double pu = MgTu/m;
			sumMgTu += MgTu;
			c += 1 - pow((1-pu), k);
			loopEstimateEPT++;
		}
		c /= loop;
		if(c > lb)
			break;
		lb /= 2;
	}
	//BuildHypergraphR_CIM(lastR);
	return 0.5 * c * n;
}

void TimGraph::BuildHyperGraph2_CIM(double epsilon, double ept)
{
	ASSERT(ept > 0);
	//int64 R = (8+2 * epsilon) * ( n * log(n) +  n * log(2)  ) / ( epsilon * epsilon * ept) /10;
	int64 R = (2 + epsilon) * ( n * log(n) ) / ( epsilon * epsilon * ept);
	cout << "[info] # of RR sets for computing KPT = " << R << endl;
	BuildHypergraphR_CIM(R);
}
	

void TimGraph::BuildHyperGraph3_CIM(double epsilon, double opt)
{
	ASSERT(opt > 0);
	//int64 R = (8+2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, k) ) / ( epsilon * epsilon * opt)/4;
	int64 R = (8+2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, k) ) / ( epsilon * epsilon * opt);
	cout << "[info] # of RR sets for mining seeds = " << R << endl;
	BuildHypergraphR_CIM(R);
}

/******************************************/
/******************************************/
/******************************************/
/******************************************/

void TimGraph::writeSeedsToFile_CIM(double epsilon)
{
	output_file_name = dataset + "/output/seeds_rrCIM_" + to_string(k);
	output_file_name += "_" + double_to_string(epsilon);
	output_file_name += "_" + double_to_string(qao) + "_" + double_to_string(qab) + "_" + double_to_string(qbo) + "_" + double_to_string(qba);
	if (ignore_B)
		output_file_name += "_1";
	else
		output_file_name += "_0";
	output_file_name += to_string((long)std::time(nullptr)) + ".txt";

	ofstream my_file;
	my_file.open(output_file_name.c_str());
	double spread = 0;
	if (my_file.is_open()) {
		for (int i = 0; i < seedSet.size(); ++i) {
			//my_file << i+1 << "\t" << seedSet.at(i) << endl;
			my_file << seedSet.at(i) << endl;
		}
		my_file.close();
	} else {
		cerr << "[error] unable to open output file; no output was written" << output_file_name << endl;
		exit(1);
	}
}


void TimGraph::writeSeedsToFile_SIM(bool is_fast_verion, double epsilon)
{
	if (!is_fast_verion)
		output_file_name = dataset + "/output/seeds_rrSIM_" + to_string(k);
	else
		output_file_name = dataset + "/output/seeds_rrSIMfast_" + to_string(k);
	output_file_name += "_" + double_to_string(epsilon);
	output_file_name += "_" + double_to_string(qao) + "_" + double_to_string(qab) + "_" + double_to_string(qbo) + "_" + double_to_string(qba);
	if (ignore_B)
		output_file_name += "_1_";
	else
		output_file_name += "_0_";
	output_file_name += to_string((long)std::time(nullptr)) + ".txt";

	ofstream my_file;
	my_file.open(output_file_name.c_str());
	double spread = 0;
	if (my_file.is_open()) {
		for (int i = 0; i < seedSet.size(); ++i) {
			//my_file << i+1 << "\t" << seedSet.at(i) << endl;
			my_file << seedSet.at(i) << endl;
		}
		my_file.close();
	} else {
		cerr << "[error] unable to open output file; no output was written" << output_file_name << endl;
		exit(1);
	}
}


