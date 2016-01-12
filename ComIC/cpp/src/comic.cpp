#include <unordered_set>
#include "mc.h"
#include "memoryusage.h"
#include "anyoption.h"
#include "timgraph.h"


AnyOption* read_options(int argc, char* argv[]) 
{
	// read the command line options
	AnyOption *opt = new AnyOption();

	// ignore POSIX style options
	opt->noPOSIX(); 
	opt->setVerbose(); /* print warnings about unknown options */
	opt->autoUsagePrint(true); /* print usage for bad options */

	/* 3. SET THE USAGE/HELP   */
	opt->addUsage("" );
	opt->addUsage("Usage: ");
	opt->addUsage("");
	opt->addUsage(" -help: print this help ");
	opt->addUsage(" -c <config_file>: specify config file ");
	opt->addUsage(" -dataset <path_to_dataset>: specify the dataset name");
	opt->addUsage(" -aseeds <path_to_A_seeds>: specify the path to A seed set file");
	opt->addUsage(" -bseeds <path_to_B_seeds>: specify the path to B seed set file");
	opt->addUsage("");

	/* 4. SET THE OPTION STRINGS/CHARACTERS */
	opt->setOption("phase");
	opt->setCommandOption("dataset");
	opt->setOption("aseeds");
	opt->setOption("qao");
	opt->setOption("qab");
	opt->setOption("qbo");
	opt->setOption("qba");
	opt->setOption("kA");
	opt->setOption("kB");
	opt->setOption("selectBSeeds");
	opt->setOption("ignoreBSeeds");
	opt->setOption("bseeds");
	opt->setOption("algo");
	opt->setOption("epsilon");
	opt->setOption("pr");

	/* for options that will be checked only on the command line not in option/resource file */
	opt->setCommandFlag("help");
	opt->setCommandOption("c");

	/* go through the command line and get the options  */
	opt->processCommandArgs(argc, argv);

	/* 6. GET THE VALUES */
	if(opt->getFlag( "help" )) {
		opt->printUsage();
		delete opt;
		exit(0);
	}

	const char *configFile = opt->getValue("c");
	if (configFile == NULL) {
		cout << "[error] config file not specified!" << endl;
		opt->printUsage();
		delete opt;
		exit(1);
	}

	opt->processFile(configFile);
	opt->processCommandArgs( argc, argv );

	cout << "[info] config file processed: " << configFile << endl;
	return opt;
}


void select_B_seeds(string folder, int kb, int n) 
{
	unordered_set<int> bset;
	bset.clear();
	while (bset.size() < kb) {
		int node = rand() % n;
		bset.insert(node);
	}
	ofstream my_file;
	string filename = folder + "/bseeds.txt";
	my_file.open(filename.c_str());
	for (auto it = bset.begin(); it != bset.end(); ++it) {
		my_file << *it << endl; // write each B seed to bseeds.txt
	}
	my_file.close();
}


void computeTrueSpreadByMC_CIM(int k, string aseeds_file_name, string dataset, vector<double> qq, string bseeds, bool ignore_B)
{
	cout << "[info] compute true spread using Monte Carlo (CompInfMax)" << endl;
	MonteCarlo mc(dataset, "/graph.txt");
	mc.setParametersMC(k, qq, ignore_B, dataset, bseeds);
	// read A-seeds
	mc.aSeeds.clear();
	ifstream afile(aseeds_file_name.c_str(), ios::in);
	if (afile.is_open()) {
		int idx = 0;
		while (!afile.eof()) {
			string line;
			getline (afile, line);
			if (line == "") continue;
			int u = atoi(line.c_str());
			//seedArr[idx++] = u;
			mc.aSeeds.push_back(u);
		}
		afile.close();
		cout << "[info] A-seeds read from file: " << endl;
		for (int i = 0; i < mc.aSeeds.size(); i++)
			cout << mc.aSeeds[i] << " ";
		cout << endl;
	} else {
		cout << "[error] unable to open file: " << aseeds_file_name << endl;
		ASSERT(false);
	}

	// read B-seeds
	mc.readBSeedsMC();
	// computing \sigma_A(S_A, \emptyset)
	double cov_0 = mc.compute_coverage_comp(mc.bSeeds, 0);
	cout << "[info] spread of A seeds ALONE = " << cov_0 << endl;	
	// estimate spread
	for (int i = 0; i < mc.bSeeds.size(); i++) {
		if (i == 0 ||  (i+1) % 10 == 0) {
			double cov = mc.compute_coverage_comp(mc.bSeeds, i+1) - cov_0;
			cout << cov << endl; 
		}
	}
}


void computeTrueSpreadByMC(int k, int *seedSet, string dataset, vector<double> qq, string bseeds, bool ignore_B)
{
	cout << "[info] compute true spread using Monte Carlo (SelfInfMax)" << endl;
	MonteCarlo mc(dataset, "/graph.txt");
	mc.setParametersMC(k, qq, ignore_B, dataset, bseeds);
	mc.readBSeedsMC();
	printf("\n");
	for (int i = 0; i < k; i++) {
		if (i == 0 || (i+1) % 10 == 0) {
			double cov = mc.compute_coverage(seedSet, i+1);
			cout << i << "\t" << cov << endl; 
		}
	}	
}


int main(int argc, char ** argv) 
{
 	cout << endl;
    AnyOption *opt = read_options(argc, argv);
    
    string dataset = string(opt->getValue("dataset"));
    cout << "[param] *** dataset = " << dataset << " ***" << endl;
    int kA = atoi(opt->getValue("kA"));
    int kB = atoi(opt->getValue("kB"));
    cout << "[param] kA = " << kA << ", kB = " << kB << endl;
    double qao = atof(opt->getValue("qao"));
    double qab = atof(opt->getValue("qab"));
    double qbo = atof(opt->getValue("qbo"));
    double qba = atof(opt->getValue("qba"));
    vector<double> qq;
    qq.push_back(qao);
    qq.push_back(qab);
    qq.push_back(qbo);
    qq.push_back(qba);
    cout << "[param] Q-probability: " << qq[0] << ", " << qq[1] << ", " << qq[2] << ", " << qq[3] << endl;

    bool ignore_B = false, select_B = false;
    ignore_B = (atoi(opt->getValue("ignoreBSeeds"))) == 1 ? true : false;
    select_B = (atoi(opt->getValue("selectBSeeds"))) == 1 ? true : false;
    cout << "[param] ignore_B = " << ignore_B << endl;
    cout << "[param] select_B = " << select_B << endl;

    string bseeds = string(opt->getValue("bseeds"));
    cout << "[param] B seed set file = " << bseeds << endl;

    int algo = atoi(opt->getValue("algo"));
    int phase = atoi(opt->getValue("phase"));

    clock_t begin;
    clock_t end;

    // read seeds from file (instead of mining) for computing \sigma_A(S_A, S_B)
    if (algo == 100) {
    	cout << "[info] algorithm: (100) reading seeds and compute coverage SelfInfMax" << endl;

		// read A-seeds from file aseeds:
		string aseeds = string(opt->getValue("aseeds"));
		int *seedArr = new int[kA];
		ifstream afile(aseeds.c_str(), ios::in);
		if (afile.is_open()) {
			int idx = 0;
			while (!afile.eof()) {
				string line;
				getline (afile, line);
				if (line == "") 
					continue;
				int u = atoi(line.c_str());
				seedArr[idx++] = u;
			}
			afile.close();
			cout << "[info] A-seeds read from file: " << endl;
			for (int i = 0; i < kA; i++)
				cout << seedArr[i] << " ";
			cout << endl;
		} else {
			cout << "[error] unable to open file: " << aseeds << endl;
			ASSERT(false);
		}
		
		computeTrueSpreadByMC(kA, seedArr, dataset, qq, bseeds, ignore_B);
		delete seedArr;
    }

    // reading seeds and compute coverage CompInfMax
    if (algo == 104) {
    	cout << "[info] algorithm: (104) reading seeds and compute coverage CompInfMax" << endl;
    	string aseeds_file_name = (opt->getValue("aseeds"));
    	computeTrueSpreadByMC_CIM(kA, aseeds_file_name, dataset, qq, bseeds, ignore_B);
    }

    // 1: Greedy + CELF + Monte Carlo for SIM
    if (algo == 1)  {
    	cout << "[info] algorithm: (1) Monte Carlo for SelfInfMax" << endl;

    	MonteCarlo mc(dataset, "/graph.txt");
    	ASSERT(kA > 0);
	    ASSERT(kA <= mc.n);

	    if (select_B) {
			cout << "[info] select B seeds: YES" << endl;
			select_B_seeds(dataset, kB, mc.n);
		} else
			cout << "[info] select B seeds: NO" << endl;

		mc.setParametersMC(kA, qq, ignore_B, dataset, bseeds);
	    mc.readBSeedsMC();
	    
	    if (phase == 1) {
	    	cout << endl << "[info] phase: mining seeds for A ..." << endl;
		    begin = clock();
		    mc.mineSeedsMC();
		    end = clock();
		    disp_mem_usage("");
			mc.write_seeds_to_file(false);
	    } else if (phase == 2) {
	    	string aSeedsFile = string(opt->getValue("aseeds")); 
	    	cout << endl << "[info] phase: estimating true spread of seeds in " << aSeedsFile << endl;
	    	mc.estimate_spread_from_file(aSeedsFile);
	    } 
    } 

    // RR-SIM
    if (algo == 2) 
    {
    	cout << "[info] algorithm: (2) RR-SIM" << endl;

    	TimGraph tim(dataset, "/graph.txt");
    	ASSERT(kA > 0);
    	ASSERT(kA <= tim.n);
    	tim.setParametersIG(kA, qq, ignore_B, dataset, bseeds);
    	tim.readBSeedsIG();
		
		double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;

		begin = clock();
		tim.mineSeeds_TIM_SIM(epsilon);
		end = clock();
		disp_mem_usage("");
		tim.writeSeedsToFile_SIM(false, epsilon);

	    int *seedArr = new int[kA];
	    for (int i = 0; i < kA; ++i)
	    	seedArr[i] = tim.seedSet.at(i);
	    computeTrueSpreadByMC(kA, seedArr, dataset, qq, bseeds, ignore_B);
	    delete seedArr;  
    }	


    // RR-SIM-FAST
    if (algo == 3) 
    {
    	cout << "[info] algorithm: (3) RR-SIM-FAST" << endl;

    	TimGraph tim(dataset, "/graph.txt");
    	ASSERT(kA > 0);
    	ASSERT(kA <= tim.n);
    	tim.setParametersIG(kA, qq, ignore_B, dataset, bseeds);
    	tim.readBSeedsIG();

    	double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;

		begin = clock();
		tim.mineSeeds_TIM_SIM_fast(epsilon);
		end = clock();
		disp_mem_usage("");
		
		tim.writeSeedsToFile_SIM(true, epsilon);

		int *seedArr = new int[kA];
	    for (int i = 0; i < kA; ++i)
	   		seedArr[i] = tim.seedSet.at(i);
	    computeTrueSpreadByMC(kA, seedArr, dataset, qq, bseeds, ignore_B);
	    delete seedArr;
    }

/*
    // RR-SIM-FAST, sandwich
    if (algo == 33) 
    {
    	cout << "[info] algorithm: (33) RR-SIM-FAST with Sandwich Approximation" << endl;

    	TimGraph tim(dataset, "/graph.txt");
    	double epsilon = atof(opt->getValue("epsilon"));
    	cout << "[param] epsilon = " << epsilon << endl;

    	tim.setParametersIG(50, qq, ignore_B, dataset, bseeds);
    	tim.readBSeedsIG();

		// true setting, as in config file
		tim.printQQ();
		tim.mineSeeds_TIM_SIM_fast(epsilon);
		int *seedArr = new int[50];
	    for (int i = 0; i < 50; ++i)
	   		seedArr[i] = tim.seedSet.at(i);
	    computeTrueSpreadByMC(50, seedArr, dataset, qq, bseeds, ignore_B);
	    
	    // upper bound
	    tim.qbo = qq[3];  // increase qbo to qba
	    tim.printQQ();
	    tim.seedSet.clear();
	    tim.mineSeeds_TIM_SIM_fast(epsilon);
	    for (int i = 0; i < 50; ++i)
	   		seedArr[i] = tim.seedSet.at(i);
	   	computeTrueSpreadByMC(50, seedArr, dataset, qq, bseeds, ignore_B);

	   	// lower bound
	   	tim.qbo = qq[2];  // decrease qba = qbo
	   	tim.qba = qq[2];
	   	tim.printQQ();
	   	tim.seedSet.clear();
	    tim.mineSeeds_TIM_SIM_fast(epsilon);
	    for (int i = 0; i < 50; ++i)
	   		seedArr[i] = tim.seedSet.at(i);
	   	computeTrueSpreadByMC(50, seedArr, dataset, qq, bseeds, ignore_B);

	    delete seedArr;
    }
*/

    // 11: Greedy + CELF + Monte Carlo for CIM
    if (algo == 11) 
    {
    	cout << "[info] algorithm: (11) Monte Carlo for CompInfMax" << endl;
    	MonteCarlo mc(dataset, "/graph.txt");
    	ASSERT(kB > 0);
	    ASSERT(kB <= mc.n);

	    mc.setParametersMC(kB, qq, ignore_B, dataset, "dummy.txt");
	    string aSeedsFile = string(opt->getValue("aseeds")); 
	    mc.readASeedsMC_comp(aSeedsFile);
	    kA = (int)mc.aSeeds.size();
	    int *ss = new int[kA];
	    for (int i = 0; i < kA; i++) 
	    	ss[i] = mc.aSeeds.at(i);
	    double baseSpread = mc.compute_coverage(ss, kA);
	    cout << "[info] spread of A seeds alone = " << baseSpread << endl;
	    delete ss;

	    begin = clock();
	    mc.mineSeedsMC_comp(baseSpread);
	    end = clock();
		disp_mem_usage("");
		//mc.write_seeds_to_file(true);
    }

    // 14: RR-CIM
    if (algo == 14) 
    {
    	cout << "[info] algorithm: (14) RR-CIM" << endl;

    	TimGraph tim(dataset, "/graph.txt");
    	ASSERT(kB > 0); 
    	ASSERT(kB <= tim.n);

	    tim.setParametersIG(kB, qq, ignore_B, dataset, "dummy.txt");
	    string str = string(opt->getValue("aseeds"));
	    tim.readASeedsIG_CIM(str);
	    double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;
	    disp_mem_usage("");

	    MonteCarlo mc(dataset, "/graph.txt");
	    kA = (int) tim.aSeeds.size();
	    int *ss = new int[kA];
	    for (int i = 0; i < kA; i++) {
	    	int z = tim.aSeeds.at(i);
	    	ss[i] = z; 
	    	mc.aSeeds.push_back(z);
	    }
	    
	    mc.setParametersMC(0, qq, ignore_B, dataset, "dummy.txt");
	    double baseSpread = mc.compute_coverage(ss, kA);
	    cout << "[info] spread of A seeds alone = " << baseSpread << endl;
	    delete ss;
	    disp_mem_usage("");


	    begin = clock();
	    tim.mineSeeds_CIM(epsilon);
	    end = clock();
	    disp_mem_usage("");
	    tim.writeSeedsToFile_CIM(epsilon);

	    // evalute true spread using Monte Carlo simulations
	    vector<int> vb = tim.seedSet;
	    for (int i = 0; i < vb.size(); i++) {
	    	double d = mc.compute_coverage_comp(vb, i+1);
	    	cout << "\t round = " << i+1 << ", node = " << vb.at(i) << ", total = " << d << endl; 
	    }
    }

/*
    // 34: RR-CIM with Sandwich
    if (algo == 34) 
    {
    	cout << "[info] algorithm: (34) RR-CIM with SANDWICH" << endl;

    	TimGraph tim(dataset, "/graph.txt");

	    tim.setParametersIG(50, qq, ignore_B, dataset, "dummy.txt");
	    string str = string(opt->getValue("aseeds"));
	    tim.readASeedsIG_CIM(str);
	    double epsilon = atof(opt->getValue("epsilon"));
		cout << "[param] epsilon = " << epsilon << endl;
	   
	    MonteCarlo mc(dataset, "/graph.txt");
	    kA = (int) tim.aSeeds.size();
	    int *ss = new int[kA];
	    for (int i = 0; i < kA; i++) {
	    	int z = tim.aSeeds.at(i);
	    	ss[i] = z; 
	    	mc.aSeeds.push_back(z);
	    }
	    
	    mc.setParametersMC(0, qq, ignore_B, dataset, "dummy.txt");
	    double baseSpread = mc.compute_coverage(ss, kA);
	    cout << "[info] spread of A seeds alone = " << baseSpread << endl;
	    delete ss;

	    // true_spread
	    tim.printQQ();
	    tim.mineSeeds_CIM(epsilon);
	    // evalute true spread using Monte Carlo simulations
	    vector<int> vb = tim.seedSet;
	    for (int i = 0; i < vb.size(); i++) {
	    	double d = mc.compute_coverage_comp(vb, i+1);
	    	cout << "\t round = " << i+1 << ", node = " << vb.at(i) << ", total = " << d << endl; 
	    }

	    // upper bound
	    tim.qba = 1;  // increase qbo to qba
	    tim.printQQ();
	    tim.seedSet.clear();
	    tim.mineSeeds_CIM(epsilon);
	    vb = tim.seedSet;
	    for (int i = 0; i < vb.size(); i++) {
	    	double d = mc.compute_coverage_comp(vb, i+1);
	    	cout << "\t round = " << i+1 << ", node = " << vb.at(i) << ", total = " << d << endl; 
	    }
    }
*/

    // random (for both SIM and CIM)
    if (algo == 99) 
    {
    	cout << "[info] algorithm: Random for both SIM and CIM" << endl;
    	MonteCarlo mc(dataset, "/graph.txt");
	    mc.setParametersMC(50, qq, ignore_B, dataset, bseeds);
	    mc.readBSeedsMC();

	    // Just randomly choose 50 seeds
	    unordered_set<int> seedSet;
	    while (seedSet.size() < 50) {
	    	int x = rand() % mc.graph->n;
	    	seedSet.insert(x);
	    }

	    // push it into a vector (for later use)
	    int *seedArr = new int[50];
	    vector<int> vb;
	    int idx = 0;
	    for (auto it = seedSet.begin(); it != seedSet.end(); ++it) {
	    	seedArr[idx++] = *it;
	    	vb.push_back(*it);
	    }

	    // spread for SIM
	 	for (int i = 0; i < 50; i++) {
			if (i == 0 || (i+1) % 10 == 0) {
				double cov = mc.compute_coverage(seedArr, i+1);
				cout << (i+1) << "\t" << seedArr[i] << "\t" << cov << endl; 
			}
		}	

/*
	    // spread for CIM
	    string aseeds_file_name = string(opt->getValue("aseeds"));
	    mc.readASeedsMC_comp(aseeds_file_name);
	    //for (int i = 0; i < kA; i++)
	    //	seedArr[i] = mc.aSeeds.at(i);
	    //double baseSpread = mc.compute_coverage(seedArr, kA);
	    //cout << "[info] spread of A-seeds ALONE = " << baseSpread << endl;

	    for (int i = 0; i < vb.size(); i++)  {
	    	if (i == 0 || (i+1) % 10 == 0) {
	    		double cov = mc.compute_coverage_comp(vb, i+1);
	    		cout << (i+1) << "\t" << cov << endl; 
	    	}
	    }
*/

	    delete seedArr;
    }

    // high degree
    if (algo == 98) {
    	cout << "[info] algorithm: Highest Degree for both SIM and CIM!" << endl;
    	MonteCarlo mc(dataset, "/graph.txt");
	    mc.setParametersMC(50, qq, ignore_B, dataset, bseeds);
	    mc.readBSeedsMC();
		// find high-degree seeds!	    
	    //begin = clock(); //end = clock();
	    vector<int> vec = mc.graph->findHighDegree(kA);
	    
	    // spread for SIM
	    int *seedArr = new int[kA];
	    for (int i = 0; i < kA; i++) 
	    	seedArr[i] = vec.at(i);

	    printf("\n");
		for (int i = 0; i < kA; i++) {
			if (i == 0 || (i+1) % 10 == 0) {
				double cov = mc.compute_coverage(seedArr, i+1);
				cout << i << "\t" << cov << endl; 
			}
		}

/*
	    // spread for CIM
	    string aseeds_file_name = string(opt->getValue("aseeds"));
	    //mc.aSeeds.clear();
	    mc.readASeedsMC_comp(aseeds_file_name);

	    for (int i = 0; i < vec.size(); i++)  {
	    	if (i == 0 || (i+1) % 10 == 0) {
	    		double cov = mc.compute_coverage_comp(vec, i+1);
	    		cout << (i+1) << "\t" << cov << endl; 
	    	}
	    }
*/
	    delete seedArr;
    }

    // page rank
    if (algo == 97) 
    {
    	cout << "[info] algorithm: PageRank for both CIM and SIM!" << endl;
    	MonteCarlo mc(dataset, "/graph.txt");
	    mc.setParametersMC(50, qq, ignore_B, dataset, bseeds);
	    mc.readBSeedsMC();

    	string prSeedFileName = string(opt->getValue("pr"));
    	vector<int> prSeeds;

    	// read seed set 
    	FILE *fin = fopen(prSeedFileName.c_str(), "r");
    	ASSERT(fin != NULL);
    	for (int i = 0; i < 50; i++) {
    		int x;
    		fscanf(fin, "%d", &x);
    		prSeeds.push_back(x);
    	}
    	fclose(fin);

    	// spread for SIM
    	int *seedArr = new int[50];
		for (int i = 0; i < 50; ++i)
	    	seedArr[i] = prSeeds.at(i);
	   	for (int i = 0; i < 50; i++) {
			if (i == 0 || (i+1) % 10 == 0) {
				double cov = mc.compute_coverage(seedArr, i+1);
				cout << (i+1) << "\t" << seedArr[i] << "\t" << cov << endl; 
			}
		}	
    	delete seedArr;

		/*
    	// spread for CIM
    	string aseeds_file_name = string(opt->getValue("aseeds"));
	    mc.aSeeds.clear();
	    mc.readASeedsMC_comp(aseeds_file_name);

	    for (int i = 0; i < prSeeds.size(); i++)  {
	    	if (i == 0 || (i+1) % 10 == 0) {
	    		double cov = mc.compute_coverage_comp(prSeeds, i+1);
	    		cout << (i+1) << "\t" << cov << endl; 
	    	}
	    }
	    */
    }


 	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "[info] *** running time for seed mining = " << elapsed_secs << " sec ***" << endl << endl;

    delete opt;
    cout << endl;
    return 0;
}



