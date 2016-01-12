#include "infgraph.h"
#include <algorithm>

void InfGraph::setParametersIG(int k, vector<double> qq, bool ignore, string folder, string bseeds)
{
	ASSERT((int)qq.size() == 4)
	this->k = k;
	this->qao = qq[0];
	this->qab = qq[1];
	this->qbo = qq[2];
	this->qba = qq[3];
	ignore_B = ignore;
	b_seeds_file_name = bseeds;

	output_file_name = folder + "/output/seeds_tim_" + to_string(k);
	output_file_name += "_" + double_to_string(qao) + "_" + double_to_string(qab) + "_" + double_to_string(qbo) + "_" + double_to_string(qba);
	if (ignore_B)
		output_file_name += "_1";
	else
		output_file_name += "_0";
	output_file_name += ".txt";
}

void InfGraph::printQQ()
{
	cout << "[info] QQ = " << qao << ", " << qab << ", " << qbo << ", " << qba << endl;
}


void InfGraph::readBSeedsIG()
{
	bSeeds.clear();
	if (ignore_B) {
		cout << "[info] config file indicates that B seeds should be ignored!" << endl;
		return;
	}

	ifstream myfile(b_seeds_file_name.c_str(), ios::in);
	if (myfile.is_open()) {
		while (! myfile.eof()) {
			string line;
			getline (myfile, line);
			if (line == "")
				continue;
			int u = atoi(line.c_str());
			bSeeds.push_back(u);
		}
		myfile.close();
	} else {
		cout << "[error] unable to open file: " << graph->folder << "/bseeds.txt" << endl;
		exit(1);
	}

	std::sort(bSeeds.begin(), bSeeds.end());

	cout << "[info] B-seeds are: ";
	for (auto i : bSeeds)
		cout << i << " ";
	cout << endl;
}



/**
 * Generate RR sets
 */
void InfGraph::BuildHypergraphR(int64 theta)
{
	hyperId = theta;

	hyperG.clear();
	for(int i = 0; i < n; i++)
		hyperG.push_back(vector<int>());

	hyperGT.clear();
	while((int)hyperGT.size() <= theta)
		hyperGT.push_back(vector<int>());

	for(int i = 0; i < theta; i++) {
		BuildSingleRRSet(rand() % n, i, true);
	}

	int totAddedElement = 0;
	for(int i = 0; i < theta; i++){
		for(int t : hyperGT[i]) {
			hyperG[t].push_back(i);
			totAddedElement++;
		}
	}

	ASSERT(hyperId == theta);
}


void InfGraph::BuildHypergraphR_fast(int64 theta)
{
	hyperId = theta;

	hyperG.clear();
	for(int i = 0; i < n; i++)
		hyperG.push_back(vector<int>());

	hyperGT.clear();
	while((int)hyperGT.size() <= theta)
		hyperGT.push_back(vector<int>());

	for(int i = 0; i < theta; i++) {
		BuildSingleRRSet_fast(rand() % n, i, true);
	}

	int totAddedElement = 0;
	for(int i = 0; i < theta; i++) {
		for(int t : hyperGT[i]) {
			hyperG[t].push_back(i);
			totAddedElement++;
		}
	}

	ASSERT(hyperId == theta);
}


void InfGraph::BuildSeedSetIG()
{
	vector<int> degree;
	vector<int> visit_local(hyperGT.size());

	seedSet.clear();
	for(int i = 0; i < n; i++)
		degree.push_back(hyperG[i].size());

	ASSERT(k > 0);
	ASSERT(k < degree.size());
	for(int i = 0; i < k; i++) {
		auto t = max_element(degree.begin(), degree.end());
		int id = t - degree.begin();
		seedSet.push_back(id);
		degree[id] = 0;
		for(int t : hyperG[id]) {
			if(!visit_local[t]) {
				visit_local[t]=true;
				for(int item:hyperGT[t]) {
					degree[item]--;
				}
			}
		}
	}
}


// compute influence spread
double InfGraph::computeInfluenceHyperGraph()
{
	unordered_set<int> s;
	for(auto t : seedSet){
		for(auto tt : hyperG[t]){
			s.insert(tt);
		}
	}
	double inf = (double) n * s.size() / hyperId;
	return inf;
}


// print seed sets
void InfGraph::printSeeds() 
{
	int count = 1;
	for (auto s : seedSet) {
		cout << "[seed] rank = " << count << ", id = " << s << endl;
		count += 1;
	}
}



int InfGraph::BuildSingleRRSet(int uStart, int rrSetID, bool addHyperEdge)
{
	vector<int> touched;
	int n_visit_edge = 1;

	std::fill(alpha_A, alpha_A + n, -1.0);  // init them to all -1
	std::fill(alpha_B, alpha_B + n, -1.0);  // init them to all -1	
	memset(status_B, 0, n*sizeof(int));
	
	// step 1: forward simulation from S_B
	Q.clear();
	for (auto x : bSeeds) {
		Q.push_back(x);
		status_B[x] = ADOPTED;
		alpha_A[x] = (double) rand() / (double) RAND_MAX;
	}

	while (!Q.empty()) {
		int u = Q.front();
		touched.push_back(u);  // conservative, but safe

		for(int j = 0; j < graph->gO[u].size(); j++) {
			int v = graph->gO[u][j];
			if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
				continue;  
			
			if (alpha_A[v] == -1.0) {
				alpha_A[v] = (double) rand() / (double) RAND_MAX;
				alpha_B[v] = (double) rand() / (double) RAND_MAX;
			}

			// test if edge is live
			double coin = (double) rand() / (double) RAND_MAX;
			if (coin <= graph->probO[u][j]) {
				graph->localEdgeStatus[u][v] = LIVE;
				if (alpha_B[v] < qbo) {
					status_B[v] = ADOPTED;
					Q.push_back(v); // v becomes B-adopted, insert into queue
				} else {
					status_B[v] = REJECTED;
				}
			} else {
				graph->localEdgeStatus[u][v] = BLOCKED;
			}
		}

		Q.pop_front();
	}

	// step 2: generate the RR set
	Q.clear();
	memset(visited,  0, n*sizeof(bool));
	Q.push_back(uStart);
	visited[uStart] = true;  // put into queue means visited
	
	if (addHyperEdge) {
		ASSERT(hyperGT.size() > rrSetID);
		hyperGT[rrSetID].push_back(uStart); // add to RR-set!
	}

	while (!Q.empty()) {
		int u = Q.front(); 
		//touched.push_back(u);
		bool flag = false;  // TRUE if this node is qualified and we need to examine its in-neighbours
		if (status_B[u] == ADOPTED) {
			if (alpha_A[u] <= qab) {
				flag = true;
			}
		} else { 
			if (alpha_A[u] == -1) {
				alpha_A[u] = (double) rand() / (double) RAND_MAX;
			}
			if (alpha_A[u] <= qao) {
				flag = true;
			}
		}

		if (flag) {
			int nn = (int) graph->gT[u].size();
			for(int j = 0; j < nn; j++) {
				int w = graph->gT[u][j]; // edge w-->u
				n_visit_edge++;
				if (visited[w])
					continue;

				int prevEdgeStatus = graph->localEdgeStatus[w][u];
				if (prevEdgeStatus == BLOCKED)
					continue;     // blocked edge, skip this node
				else if (prevEdgeStatus == INACTIVE) {
					double coin = (double) rand() / (double) RAND_MAX;
					if (coin > graph->probT[u][j])  {
						continue; // blocked edge, skip this node
					}
				} 

				// at this point, (w,u) must be live, so add w to BFS queue
				visited[w] = true;
				Q.push_back(w);  

				if (addHyperEdge) {
					ASSERT((int)hyperGT.size() > rrsetID);
					hyperGT[rrSetID].push_back(w); // add to RR-set!
				}
			} // end-for			
		}

		Q.pop_front();
	}

	graph->reset_local_edge_status(touched);

	return n_visit_edge;
}



int InfGraph::BuildSingleRRSet_fast(int uStart, int rrSetID, bool addHyperEdge)
{
	vector<int> touched;
	memset(visited,  0, n*sizeof(bool));
	
	// first round of backaward BFS from random start node
	Q.clear();
	Q.push_back(uStart);
	visited[uStart] = true;
	vector<int> tvec;
	while (!Q.empty()) {
		int u = Q.front();
		tvec.push_back(u);
		// examine all in-neighbours
		for (int j = 0; j < graph->gT[u].size(); j++)  {
			int w = graph->gT[u][j]; // edge is w --> u
			if (visited[w])
				continue;
			// test if edge (w,u) is live or blocked
			touched.push_back(w);
			double coin = (double) rand() / (double) RAND_MAX;
			if (coin <= graph->probT[u][j]) {
				graph->localEdgeStatus[w][u] = LIVE;
				Q.push_back(w);
				visited[w] = true;
			} else {
				graph->localEdgeStatus[w][u] = BLOCKED;
			}
		}
		Q.pop_front();
	}

	// get T1 \cap S_B
	std::sort(tvec.begin(), tvec.end());
	vector<int> joined;
	joined.resize(tvec.size()); 
	vector<int>::iterator it = std::set_intersection(tvec.begin(), tvec.end(), bSeeds.begin(), bSeeds.end(), joined.begin());
	joined.resize(it - joined.begin()); 

	memset(status_B, 0, n*sizeof(int));
	std::fill(alpha_A, alpha_A + n, -1.0);  // initialize them to all -1
	std::fill(alpha_B, alpha_B + n, -1.0);  // initialize them to all -1

	// forward simulation from T1 \cap S_B (only when non-empty)
	if (joined.size() > 0) {
		Q.clear();
		for (auto x : joined) {
			Q.push_back(x);
			status_B[x] = ADOPTED;
			alpha_A[x] = (double) rand() / (double) RAND_MAX;
		}

		while (!Q.empty()) {
			int u = Q.front();
			touched.push_back(u);
			for(int j = 0; j < graph->gO[u].size(); j++) {
				int v = graph->gO[u][j]; // edge is u --> v
				if (status_B[v] == ADOPTED || status_B[v] == REJECTED)
					continue;  
				int prevEdgeStatus = graph->localEdgeStatus[u][v];
				if (prevEdgeStatus == BLOCKED)
					continue;
				
				if (prevEdgeStatus == INACTIVE) {
					double coin = (double) rand() / (double) RAND_MAX;
					if (coin <= graph->probO[u][j]) {
						graph->localEdgeStatus[u][v] = LIVE;
					} else {
						graph->localEdgeStatus[u][v] = BLOCKED;
						continue;
					}
				}
				// at this point, (u,v) must be live!
				// if its alpha values are not drawn, sample it
				//if (alpha_A[v] == -1.0) {
					alpha_A[v] = (double) rand() / (double) RAND_MAX;
					alpha_B[v] = (double) rand() / (double) RAND_MAX;
				//}

				if (alpha_B[v] < qbo) {
					status_B[v] = ADOPTED;
					Q.push_back(v); // v becomes B-adopted, insert into queue
				} else {
					status_B[v] = REJECTED;
				}
			}
			Q.pop_front();
		} // endwhile
	} // endif

	// second round of BFS from the random start node
	Q.clear();
	memset(visited,  0, n*sizeof(bool));
	Q.push_back(uStart);
	visited[uStart] = true;  
	if (addHyperEdge) {
		ASSERT(hyperGT.size() > rrSetID);
		hyperGT[rrSetID].push_back(uStart); // add to RR-set!
	}
	int n_visit_edge = 1;

	while (!Q.empty()) {
		int u = Q.front(); 
		//touched.push_back(u);
		bool flag = false;  // TRUE if this node is qualified and we need to examine its in-neighbours
		
		if (status_B[u] == ADOPTED) {
			if (alpha_A[u] <= qab) {
				flag = true;
			}
		} else { 
			if (alpha_A[u] == -1.0) {
				alpha_A[u] = (double) rand() / (double) RAND_MAX;
			}
			if (alpha_A[u] <= qao) {
				flag = true;
			}
		}

		if (flag) {
			int nn = (int) graph->gT[u].size();
			for(int j = 0; j < nn; j++) {
				int w = graph->gT[u][j]; // edge w-->u
				n_visit_edge++;
				if (visited[w])
					continue;

				int prevEdgeStatus = graph->localEdgeStatus[w][u];
				if (prevEdgeStatus == BLOCKED)
					continue;     // blocked edge, skip this node

				// at this point, (w,u) must be live, so add w to BFS queue
				visited[w] = true;
				Q.push_back(w);  

				if (addHyperEdge) {
					ASSERT((int)hyperGT.size() > rrsetID);
					hyperGT[rrSetID].push_back(w); // add to RR-set!
				}
			} // end-for			
		} // end-if

		Q.pop_front();
	} // end-while

	// finally
	graph->reset_local_edge_status(touched);
	return n_visit_edge;
}


/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/
/**********************************//**********************************//**********************************/


void InfGraph::BuildHypergraphR_CIM(int64 theta)
{
	hyperId = theta;

	hyperG.clear();
	for(int i = 0; i < n; i++)
		hyperG.push_back(vector<int>());

	hyperGT.clear();
	while((int)hyperGT.size() <= theta)
		hyperGT.push_back(vector<int>());

	for(int i = 0; i < theta; i++) {
		BuildSingleRRSet_CIM(rand() % n, i, true);
	}

	int totAddedElement = 0;
	for(int i = 0; i < theta; i++) {
		for(int t : hyperGT[i]) {
			hyperG[t].push_back(i);
			totAddedElement++;
		}
	}

	ASSERT(hyperId == theta);
}


void InfGraph::readASeedsIG_CIM(string str) 
{
	aSeeds.clear();

	ifstream myfile(str.c_str(), ios::in);
	if (!myfile.is_open()) {
		cout << "[error] unable to open file: " << graph->folder << "/" << str << endl;
		exit(1);
	}

	while (!myfile.eof()) {
		string line;
		getline (myfile, line);
		if (line == "")
			continue;
		int u = atoi(line.c_str());
		aSeeds.push_back(u);
	}
	myfile.close();

	std::sort(aSeeds.begin(), aSeeds.end());

	cout << "[info] A-seeds are: ";
	for (auto i : aSeeds)
		cout << i << " ";
	cout << endl;
}


int InfGraph::BuildSingleRRSet_CIM(int uStart, int rrSetID, bool addHyperEdge) 
{
	memset(status_A, 0, sizeof(int) * n);
	memset(status_B, 0, sizeof(int) * n);
	for (int v = 0; v < n; v++) {
		alpha_A[v] = (double) rand() / (double) RAND_MAX;
		alpha_B[v] = (double) rand() / (double) RAND_MAX;
	}
	
	if (alpha_A[uStart] > qab) {
		return 1; // cannot be A-adopted, no matter what B-seeds are
	}

	// forward labeling process from A-seed set
	vector<int> touched;
	Q.clear();
	for (auto x : aSeeds) {
		status_A[x] = ADOPTED;
		Q.push_back(x);
	}
	while (!Q.empty()) 
	{
		int u = Q.front();
		Q.pop_front();
		touched.push_back(u);

		for(int j = 0; j < graph->outDeg[u]; j++) 
		{
			int v = graph->gO[u][j];
			if (status_A[v] == ADOPTED || status_A[v] == REJECTED || status_A[v] == SUSPENDED)
				continue; 
			if (status_A[u] == ADOPTED && status_A[v] == POTENTIAL) {
				status_A[v] = SUSPENDED; // in this case, v changes from Potential to Suspended
				continue;
			}

			double coin = (double) rand() / (double) RAND_MAX;
			if (coin <= graph->probO[u][j]) 
			{
				graph->localEdgeStatus[u][v] = LIVE;

				if (status_A[u] == ADOPTED) {
					if (alpha_A[v] <= qao)
						status_A[v] = ADOPTED;
					else {
						if (alpha_A[v] <= qab)
							status_A[v] = SUSPENDED;
					}
				} else if (status_A[u] == SUSPENDED || status_A[u] == POTENTIAL) {
					if (alpha_A[v] <= qab)
						status_A[v] = POTENTIAL;
				}

				/*
				if (alpha_A[v] <= qao) {
					status_A[v] = ADOPTED;
					Q.push_back(v);
				} else {
					if (alpha_A[v] <= qab) {
						if (status_A[u] == ADOPTED) {
							status_A[v] = SUSPENDED;
						}
						else if (status_A[u] == SUSPENDED || status_A[u] == POTENTIAL) {
							status_A[v] = POTENTIAL;
						}

						Q.push_back(v);

					} else {
						status_A[v] = REJECTED; // alphaA[v] > qab, no chance
					}
				}
				 */

			} else {
				graph->localEdgeStatus[u][v] = BLOCKED;
			}
		} // end-for
	} // end-while

	// check if we need to proceed at all
	if (status_A[uStart] != SUSPENDED && status_A[uStart] != POTENTIAL) {
		return 1;
	}

	// _primary_ backward BFS from uStart;
	Q.clear();
	memset(visited, 0, sizeof(bool) * n); // if node is discovered already in primary search
	memset(used, 0, sizeof(bool) * n);    // if node has been added to RR-set already (for safety)
	Q.push_back(uStart);
	visited[uStart] = true;
	int n_visit_edge = 1;

	while (!Q.empty()) 
	{
		int u = Q.front();
		Q.pop_front();

		if (status_A[u] == SUSPENDED) 
		{
			// A-suspended node, qualify for RR-set
			if (addHyperEdge && !used[u]) {
				hyperGT[rrSetID].push_back(u); 
				used[u] = true;
				//n_visit_edge += graph->inDeg[u];
			}
			// If u is further diffusible, launch a secondary BACKWARD search from u to find all B-diffusible nodes, until we hit non-B-diffusible nodes
			if (alpha_B[u] <= qbo) 
			{
				memset(discovered, 0, sizeof(bool) * n);
				deque<int> QQQ;
				QQQ.push_back(u);
				discovered[u] = true;
				
				while (!QQQ.empty()) 
				{
					int x = QQQ.front();
					QQQ.pop_front();

					for (int j = 0; j < graph->inDeg[x]; ++j) 
					{
						n_visit_edge++;  // because x is in the RR-set, so this edge counts for EPT
						int y = graph->gT[x][j]; // edge is y --> x
						if (discovered[y]) 
							continue;
						
						bool live_flag = false;
						int prevEdgeStatus = graph->localEdgeStatus[y][x];
						if (prevEdgeStatus == BLOCKED) {
							continue;
						} else if (prevEdgeStatus == LIVE) {
							live_flag = true;
						} else {
							double coin = (double) rand() / (double) RAND_MAX;
							if (coin <= graph->probT[x][j]) 
								live_flag = true;	
						}
						
						// if edge y-->x is LIVE, y is qualified to RR-set
						if (live_flag) {
							discovered[y] = true;
							if (addHyperEdge && !used[y]) {
								hyperGT[rrSetID].push_back(y);
								used[y] = true;
								//n_visit_edge += graph->inDeg[y];
							} 
							// if y is further B-diffusible, add it to queue, so its in-neighbors will also be examined
							if (status_A[y] == ADOPTED || alpha_B[y] <= qbo) {
								QQQ.push_back(y);
							}
						} // END-IF
					} // END-FOR
				} // END-WHILE (secondary search)
			}  // END-IF
		} 

		else if (status_A[u] == POTENTIAL) 
		{
			// A-potential & diffusible, add its live in-neighbors into FIFO queue
			if (alpha_B[u] <= qbo) {	
				for (int j = 0; j < graph->inDeg[u]; j++) {
					int w = graph->gT[u][j]; // edge is w --> u
					// n_visit_edge++;  // u is not added to RR-set, so this edge w->u does not count for EPT
					if (graph->localEdgeStatus[w][u] == LIVE) {
						if (!visited[w] && !used[w]) {
							Q.push_back(w);
							visited[w] = true;
						}
					} 
				} 
			}  

			// A-potential, but non-diffusible
			else  
			{
				bool good_flag = false; // indicating whether u should be put into RR-set
				unordered_set<int> q_fwd;
				unordered_set<int> q_bwd;

				// secondary FORWARD search to find all B-DIFFUSIBLE nodes reachable via live-edges
				deque<int> QQQ;
				memset(discovered, 0, sizeof(bool) * n);
				QQQ.push_back(u);
				discovered[u] = true;
				
				while (!QQQ.empty()) 
				{
					int x = QQQ.front();
					QQQ.pop_front();

					for (int j = 0; j < graph->outDeg[x]; ++j) {
						int y = graph->gO[x][j]; // edge x --> y
						if (discovered[y]) continue;
						
						bool live_flag = false;
						int prevEdgeStatus = graph->localEdgeStatus[x][y];
						if (prevEdgeStatus == BLOCKED)
							continue;
						else if (prevEdgeStatus == LIVE)
							live_flag = true;
						else {
							double coin = (double) rand() / (double) RAND_MAX;
							if (coin <= graph->probO[x][j])
								live_flag = true;
 						}

 						if (live_flag) {
 							discovered[y] = true; // regardless of B-diffusible or not, we already discovered this node
 							if (status_A[y] == ADOPTED || alpha_B[y] <= qbo) {
	 							QQQ.push_back(y);
	 							q_fwd.insert(y);
	 						}
 						}
					}
				} // END-WHILE

				if (!q_fwd.empty()) 
				{
					// secondary BACKWARD search from u, to find all A-adopted or diffusible A-suspended nodes reachable via live-edges
					QQQ.clear();
					memset(discovered, 0, sizeof(bool) * n);
					QQQ.push_back(u);
					discovered[u] = true;

					while (!QQQ.empty()) 
					{
						int x = QQQ.front();
						QQQ.pop_front();
						
						for (int j = 0; j < graph->inDeg[x]; ++j)
						{
							int y = graph->gT[x][j]; // edge is y --> x
							if (discovered[y]) continue;

							bool live_flag = false;
							int prevEdgeStatus = graph->localEdgeStatus[y][x];
							if (prevEdgeStatus == BLOCKED)
								continue;
							else if (prevEdgeStatus == LIVE)
								live_flag = true;
							else {
								double coin = (double) rand() / (double) RAND_MAX;
								if (coin <= graph->probT[x][j])
									live_flag = true;
							}

							if (live_flag) {
								discovered[y] = true;
								if ( (status_A[y] == ADOPTED) || (status_A[y] == SUSPENDED && alpha_B[y] <= qbo) ) {
									QQQ.push_back(y);
									q_bwd.insert(y);
								}							
							}
						}
					} // END-WHILE

					if (!q_bwd.empty()) {
						for (auto z : q_fwd) {
							unordered_set<int>::const_iterator iter = q_bwd.find(z);
							if (iter != q_bwd.end() && status_A[z] == SUSPENDED) {
								good_flag = true;
								break;
							}
						}

						if (good_flag && addHyperEdge && !used[u]) {
							hyperGT[rrSetID].push_back(u); 
							used[u] = true;
							//n_visit_edge += graph->inDeg[u];
						}
					}
				}

			} // END-ELSE

		} // END-ELSE-IF

/*
		else 
		{
			if (status_A[u] == ADOPTED)
				cout << "[warning]: A-adopted node encountered in the primary search of RR-CIM" << endl;
			else if (status_A[u] == REJECTED)
				cout << "[warning]: A-rejected node encountered in the primary search of RR-CIM" << endl;
			//exit(1);
		} 
*/

	} // END-WHILE (primary search)


	graph->reset_local_edge_status(touched);
	return n_visit_edge;
}



