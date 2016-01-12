#include "graph.h"

void Graph::readNM()
{
	ifstream cin((folder+"/attribute.txt").c_str());
	ASSERT(!cin == false);
	string s;
	while(cin >> s)
	{
		if(s.substr(0,2) == "n=")
		{
			n=atoi(s.substr(2).c_str());
			continue;
		}
		if(s.substr(0,2) == "m=")
		{
			m=atoi(s.substr(2).c_str());
			continue;
		}
		ASSERT(false);
	}
	cin.close();
	cout << "[info] n = " << n << ", m = " << m << endl;
}

void Graph::add_in_edge(int a, int b, double p)
{
	gT[b].push_back(a);
	probT[b].push_back(p);
	inDeg[b]++;
}

void Graph::add_out_edge(int a, int b, double p)
{
	gO[a].push_back(b);
	probO[a].push_back(p);
	outDeg[a]++;

	outEdgeStatus[a].push_back(INACTIVE);
	localEdgeStatus[a][b] = 0;
}

void Graph::readGraph_v2()
{
	string fileName = folder + graph_file;
	ifstream myfile (fileName.c_str(), ios::in);
	string delim = " \t";
	if (myfile.is_open()) {
		while (! myfile.eof() ) {
			std::string line;
			getline (myfile,line);
			if (line.empty()) continue;
			std::string::size_type pos = line.find_first_of(delim);
			int	prevpos = 0;
			
			string str = line.substr(prevpos, pos-prevpos);
			int a = std::stoi(str);
			//cout << a << endl;
			prevpos = line.find_first_not_of(delim, pos);
			pos = line.find_first_of(delim, prevpos);
			int b = std::stoi(line.substr(prevpos, pos-prevpos));
			//cout << b << endl;
			
			double p = 0;
			prevpos = line.find_first_not_of(delim, pos);
			pos = line.find_first_of(delim, prevpos);
			if (pos == string::npos) 
				p = std::stod(line.substr(prevpos));
			else
				p = std::stod(line.substr(prevpos, pos-prevpos));

			//cout << a << ", " << b << ", " << p << endl;

			ASSERT(a < n);
			ASSERT(b < n);

			hasNode[a]=true;
			hasNode[b]=true;	


			add_in_edge(a, b, p);
			add_out_edge(a, b, p);	
		}

		myfile.close();
	
	} else {
		cout << "[error] can't open graph file " << fileName << endl;
		exit(1);
	} 

	cout << "[info] finish reading graph data" << endl;
}

void Graph::readGraph()
{
	FILE *fin = fopen((folder+graph_file).c_str(), "r");
	ASSERT(fin != NULL);

	int readCnt = 0;
	for(int i=0; i<m; i++)
	{
		readCnt++;
		//cout << readCnt << endl;
		int a, b;
		double p;
		int c = fscanf(fin, "%d%d%lf", &a, &b, &p);
		//cout << a << ", " << b << ", " << p << endl;
		ASSERT(c == 3);
		ASSERTT(c == 3, a, b, p, c);

		ASSERT(a < n);
		ASSERT(b < n);
		hasNode[a]=true;
		hasNode[b]=true;
		add_in_edge(a, b, p);
		add_out_edge(a, b, p);
	}

	if(readCnt != m)
		ExitMessage("[error] m != number of edges in file " + graph_file);

	fclose(fin);
	cout << "[info] finish reading graph data" << endl;
}

void Graph::printGraph(int num_nodes)
{
	if (num_nodes > n)
		num_nodes = n;
	for (int u = 0; u < num_nodes; u++) {
		printf("%d: ", u);
		for (int i = 0; i < outDeg[u]; i++) {
			printf("(%d, %g) ", gO[u][i], probO[u][i]);
		}
		printf("\n");
	}
}


// find top-k highest degree nodes
vector<int> Graph::findHighDegree(int k) 
{
	cout << "[info] ranking nodes by degrees." << endl;
	vector<int> ret;
	vector<std::pair<int, int>> degVec;
	for (int i = 0; i < n; i++) {
		std::pair<int, int> p;
		p = std::make_pair(i, outDeg[i]);
		degVec.push_back(p);
	}

	std::sort(degVec.begin(), degVec.end(), [](const std::pair<int,int> &left, const std::pair<int,int> &right) {
    	return left.second < right.second;
	});

	for (int i = 0; i < k; i++) {
		std::pair<int, int> p = degVec.at(degVec.size() - 1 - i);
		//cout << "degree = " << p.second << endl;
		ret.push_back(p.first);
	}

	return ret;
}

// I might just use ICDM'11 code to compute PageRank score
vector<int> Graph::findPageRank(int k) 
{
	vector<int> retVec;
	//
	//
	//
	return retVec;
}
