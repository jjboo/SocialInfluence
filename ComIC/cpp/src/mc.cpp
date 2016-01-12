#include "mc.h"
#include <algorithm>

void MonteCarlo::setParametersMC(int k_A, vector<double> qq, bool ignore, string folder, string bseeds) 
{
	ASSERT((int)qq.size() == 4);
	k = k_A;
	qao = qq[0];
	qab = qq[1];
	qbo = qq[2];
	qba = qq[3];
	ignore_B = ignore;
	b_seeds_file_name = bseeds;
	output_file_name = folder + "/output/seeds_mc_" + to_string(k);
	output_file_name += "_" + double_to_string(qao) + "_" + double_to_string(qab) + "_" + double_to_string(qbo) + "_" + double_to_string(qba);
	if (ignore_B)
		output_file_name += "_1";
	else
		output_file_name += "_0";
	output_file_name += ".txt";
}


void MonteCarlo::readBSeedsMC()
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
			if (line == "") continue;
			int u = atoi(line.c_str());
			bSeeds.push_back(u);
		}
		myfile.close();
	} else {
		cout << "[error] unable to open file: " << graph->folder << "/bseeds.txt" << endl;
		exit(1);
	}

	// we should not need sorting in any case
	//std::sort(bSeeds.begin(), bSeeds.end());

	cout << "[info] B-seeds (from input) are: ";
	for (auto i : bSeeds)
		cout << i << " ";
	cout << endl;

	//ASSERT((int)aSeeds.size() == 0);
}


void MonteCarlo::readASeedsMC_comp(string str)
{
	aSeeds.clear();
	ifstream myfile(str.c_str(), ios::in);
	if (myfile.is_open()) {
		while (! myfile.eof()) {
			string line;
			getline (myfile, line);
			if (line == "") continue;
			int u = atoi(line.c_str());
			aSeeds.push_back(u);
		}
		myfile.close();
	} else {
		cout << "[error] unable to open file: " << graph->folder << "/" << str << endl;
		exit(1);
	}

	//std::sort(aSeeds.begin(), aSeeds.end());

	cout << "[info] A-seeds (from input) are: ";
	for (auto i : aSeeds)
		cout << i << " ";
	cout << endl;
	cout << "[info] size of A-seeds: " << aSeeds.size() << endl;

	//ASSERT((int)bSeeds.size() == 0);
}

/**
 *  baseSpread: \sigma_A(S_A, \emptyset)
 */
double MonteCarlo::mineSeedsMC_comp(double baseSpread)
{
	double  *improve = new double[graph->n];
	int     *last_update = new int[graph->n];
	int     *heap = new int[graph->n];
	vector<int> tmp_set;
	tmp_set.resize(k);

	for (int i = 0; i < n; i++) {
		heap[i] = i;
		last_update[i] = -1;
		improve[i] = (double)(n+1);
	}

	double old = 0;
	srand(time(NULL));
	bSeeds.clear();
	mg.clear();

	for (int i = 0; i < k; i++) {
		while (last_update[heap[0]] != i) {
			last_update[heap[0]] = i;
			tmp_set[i] = heap[0];
			improve[heap[0]] = compute_coverage_comp(tmp_set, i+1) - old - baseSpread;

			int x = 0;
			while (x*2 + 2 <= n-i) {
				int newx = x*2 + 1;
				if ((newx+1 < n-i) && (improve[heap[newx]] < improve[heap[newx+1]]))
					newx++;
				if (improve[heap[x]] < improve[heap[newx]]) {
					int t = heap[x];
					heap[x] = heap[newx];
					heap[newx] = t;
					x = newx;
				} else {
					break;
				}
			} // end-while
		} // end-while

		bSeeds.push_back(heap[0]);
		tmp_set[i] = heap[0];
		mg.push_back(improve[heap[0]]);
		old += improve[heap[0]];

		cout << "\tround " << i+1 << ": node = " << bSeeds[i] << ", mg = " << mg[i] << ", total = " << (old + baseSpread) <<  endl;

		heap[0] = heap[n-i-1];
		int x = 0;
		while (x*2 + 2 <= n - i) {
			int newx = x*2 + 1;
			if ((newx+1 < n-i) && (improve[heap[newx]] < improve[heap[newx+1]]))
				newx++;
			if (improve[heap[x]]<improve[heap[newx]]) {
				int t = heap[x];
				heap[x] = heap[newx];
				heap[newx] = t;
				x = newx;
			} else { 
				break;
			}
		} // endwhile	

	} // end-for

	int rep = 3;
	double final_spread = 0;
	for (int j = 0; j < rep; j++)
		final_spread += compute_coverage_comp(bSeeds, k);
	
	delete[] heap;
	delete[] last_update;
	delete[] improve;

	return final_spread / (double) rep;
}


double MonteCarlo::mineSeedsMC()
{
	double  *improve = new double[graph->n];
	int     *last_update = new int[graph->n];
	int     *heap = new int[graph->n];
	int     *tmp_set = new int[k];
	memset(tmp_set, 0, sizeof(int)*k);

	for (int i = 0; i < n; i++) {
		heap[i] = i;
		last_update[i] = -1;
		improve[i] = (double)(n+1);
	}

	double old = 0;
	srand(time(NULL));

	aSeeds.clear();
	mg.clear();

	int count_nodes = 0;
	for (int i = 0; i < k; i++) {
		while (last_update[heap[0]] != i) {
			last_update[heap[0]] = i;
			tmp_set[i] = heap[0];
			improve[heap[0]] = compute_coverage(tmp_set, i+1) - old;
			if (i==0) {
				count_nodes++;
				if (count_nodes % 1000 == 0)
					printf("1st iteration with %d nodes done...\n", count_nodes);
			}

			int x = 0;
			while (x*2 + 2 <= n-i) {
				int newx = x*2 + 1;
				if ((newx+1 < n-i) && (improve[heap[newx]] < improve[heap[newx+1]]))
					newx++;
				if (improve[heap[x]] < improve[heap[newx]]) {
					int t = heap[x];
					heap[x] = heap[newx];
					heap[newx] = t;
					x = newx;
				} else {
					break;
				} 
			} //endwhile
		} //endwhile

		aSeeds.push_back(heap[0]);
		tmp_set[i] = heap[0];
		mg.push_back(improve[heap[0]]);
		old += improve[heap[0]];

		cout << "\tround " << i+1 << ": node = " << aSeeds[i] << ", mg = " << mg[i] << ", total = " << old <<  endl;

		heap[0] = heap[n-i-1];
		int x = 0;
		while (x*2 + 2 <= n - i) {
			int newx = x*2 + 1;
			if ((newx+1 < n-i) && (improve[heap[newx]] < improve[heap[newx+1]]))
				newx++;
			if (improve[heap[x]]<improve[heap[newx]]) {
				int t = heap[x];
				heap[x] = heap[newx];
				heap[newx] = t;
				x = newx;
			} else { 
				break;
			}
		} // endwhile
	} //endfor

	printf("\n\n");

	if (aSeeds.size() < k)
		cout << "[warning] less than " << k << " seeds were selected." << endl;
	int *seed_set = new int[k];
	for (int i = 0; i < k; i++)
		seed_set[i] = aSeeds.at(i);

	int rep = 3;
	double final_spread = 0;
	for (int j = 0; j < rep; j++) {
		final_spread += compute_coverage(seed_set, k);
	}

	delete[] heap;
	delete[] last_update;
	delete[] improve;
	delete[] tmp_set;
	delete[] seed_set;

	return final_spread / (double) rep;
}



double MonteCarlo::compute_coverage_comp(vector<int> bSeedSet, int size_B)
{
	ASSERT((int)bSeedSet.size() >= size_B);
	//ASSERT((int)aSeeds.size() == 50);

	double  cov         = 0;
	double  *alpha_A    = new double[graph->n];
	double  *alpha_B    = new double[graph->n];
	int     *status_A   = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int     *status_B   = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++)  {
			alpha_A[i] = (double) rand() / (double) RAND_MAX;
			alpha_B[i] = (double) rand() / (double) RAND_MAX;
		}
		graph->reset_out_edge_status();

		// scan all A-seeds
		for (int i = 0; i < (int)aSeeds.size(); ++i) {
			int u = aSeeds.at(i);
			status_A[u] = ADOPTED;
			cov++;
			// iterate over its out-neighbors
			for(int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = (double) rand() / (double) RAND_MAX;
				if (coin <= graph->probO[u][j]) {
					graph->outEdgeStatus[u][j] = LIVE;
					if (status_A[v] != ADOPTED)
						list_A.push_back(v);
				} else {
					graph->outEdgeStatus[u][j] = BLOCKED;
				}
			}
		}

		// scan all B-seeds
		//for (auto it = bSeedSet.begin(); it != bSeedSet.end(); ++it) 
		for (int i = 0; i < size_B; i++)
		{
			int u = bSeedSet.at(i);
			status_B[u] = ADOPTED;

			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				if (graph->outEdgeStatus[u][j] == INACTIVE) {
					double coin = (double) rand() / (double) RAND_MAX;
					if (coin <= graph->probO[u][j]) {
						graph->outEdgeStatus[u][j] = LIVE;  // edge is live
						if (status_B[v] != ADOPTED)
							list_B.push_back(v);
					} else {
						graph->outEdgeStatus[u][j] = BLOCKED; // edge is blocked
					}
				}  else if (graph->outEdgeStatus[u][j] == LIVE && status_B[v] != ADOPTED) {
					list_B.push_back(v);
				}

			}
		}

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;

		while (curr_A > 0 || curr_B > 0) {
			// A-adoption test
			for (int i = 0; i < curr_A; i++) {
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == SUSPENDED || status_A[v] == ADOPTED) 
					continue;
				
				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						cov++;
						if (status_B[v] == SUSPENDED && alpha_B[v] <= qba) {
							status_B[v] = ADOPTED; // reconsider to adopt B
							examine_out_neighbors(v, &list_B, &next_B, status_B);
						}
					} else {
						status_A[v] = SUSPENDED; // A-suspended
					}
					
				} else {
					// v is already B-adopted, test with q_A|B
					if (alpha_A[v] <= qab) {
						status_A[v] = ADOPTED;
						cov++;
					} else {
						status_A[v] = SUSPENDED;
					}
				}

				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) {
					examine_out_neighbors(v, &list_A, &next_A, status_A);
				} // END-IF
			} // ENDFOR

			// B adoption test
			for (int i = 0; i < curr_B; i++) {
				int v = list_B.front();
				list_B.pop_front();
				if (status_B[v] == SUSPENDED || status_B[v] == ADOPTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						if (status_A[v] == SUSPENDED && alpha_A[v] <= qab) {
							status_A[v] = ADOPTED; // reconsideration for A!
							cov++;
							examine_out_neighbors(v, &list_A, &next_A, status_A);
						}
					} else {
						status_B[v] = SUSPENDED;
					}
					
				} else {
					status_B[v] = (alpha_B[v] <= qba) ? ADOPTED : SUSPENDED; // already A-adopted
				}

				if (status_B[v] == ADOPTED) {
					examine_out_neighbors(v, &list_B, &next_B, status_B);
				} // END-IF
			} // END-FOR

			curr_A = next_A;
			curr_B = next_B;
			next_A = next_B = 0;

		} // END-WHILE
	}

	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;

	return cov / (double) MC_RUNS;
}


double MonteCarlo::compute_coverage(int *set_A, int size_A)
{
	double  cov         = 0;
	double  *alpha_A    = new double[graph->n];
	double  *alpha_B    = new double[graph->n];
	int     *status_A   = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)
	int     *status_B   = new int[graph->n]; // 0: inactive, 1: informed, 2: suspended, 3: adopted (active)

	deque<int> list_A;  // hold the nodes currently informed of A
	deque<int> list_B;  // hold the nodes currently informed of B

	for (int r = 0; r < MC_RUNS; r++) {
		list_A.clear();
		list_B.clear();
		memset(status_A, 0, sizeof(int) * n);
		memset(status_B, 0, sizeof(int) * n);
		for (int i = 0; i < n; i++)  {
			alpha_A[i] = (double) rand() / (double) RAND_MAX;
			alpha_B[i] = (double) rand() / (double) RAND_MAX;
		}
		graph->reset_out_edge_status();

		// scan all A-seeds
		for (int i = 0; i < size_A; ++i) {
			int u = set_A[i];
			status_A[u] = ADOPTED;
			cov++;
			// iterate over its out-neighbors
			for(int j = 0; j < graph->outDeg[u]; j++) {
				int v = graph->gO[u][j];
				double coin = (double) rand() / (double) RAND_MAX;
				if (coin <= graph->probO[u][j]) {
					graph->outEdgeStatus[u][j] = LIVE;
					if (status_A[v] != ADOPTED)
						list_A.push_back(v);
				} else {
					graph->outEdgeStatus[u][j] = BLOCKED;
				}
			}
		}

		// scan all B-seeds
		for (auto it = bSeeds.begin(); it != bSeeds.end(); ++it) {
			int u = *it;
			status_B[u] = ADOPTED;

			for (int j = 0; j < graph->outDeg[u]; j++) { // iterate over its out-neighbors
				int v = graph->gO[u][j];
				if (graph->outEdgeStatus[u][j] == INACTIVE) {
					double coin = (double) rand() / (double) RAND_MAX;
					if (coin <= graph->probO[u][j]) {
						graph->outEdgeStatus[u][j] = LIVE;  // edge is live
						if (status_B[v] != ADOPTED)
							list_B.push_back(v);
					} else {
						graph->outEdgeStatus[u][j] = BLOCKED; // edge is blocked
					}
				}  else if (graph->outEdgeStatus[u][j] == LIVE && status_B[v] != ADOPTED) {
					list_B.push_back(v);
				}

			}
		}

		int curr_A = list_A.size();
		int curr_B = list_B.size();
		int next_A = 0, next_B = 0;

		while (curr_A > 0 || curr_B > 0) {
			// A-adoption test
			for (int i = 0; i < curr_A; i++) {
				int v = list_A.front();
				list_A.pop_front();
				if (status_A[v] == SUSPENDED || status_A[v] == ADOPTED) 
					continue;
				
				if (status_B[v] != ADOPTED) {
					// v is NOT B-adopted, test with q_A|0
					if (alpha_A[v] <= qao) {
						status_A[v] = ADOPTED;  // A-adopted
						cov++;
						if (status_B[v] == SUSPENDED && alpha_B[v] <= qba) {
							status_B[v] = ADOPTED; // reconsider to adopt B
							examine_out_neighbors(v, &list_B, &next_B, status_B);
						}
					} else {
						status_A[v] = SUSPENDED; // A-suspended
					}
					
				} else {
					// v is already B-adopted, test with q_A|B
					if (alpha_A[v] <= qab) {
						status_A[v] = ADOPTED;
						cov++;
					} else {
						status_A[v] = SUSPENDED;
					}
				}

				// if v adopts a product for the first time, we test its outgoing edges
				if (status_A[v] == ADOPTED) {
					examine_out_neighbors(v, &list_A, &next_A, status_A);
				} // END-IF
			} // ENDFOR

			// B adoption test
			for (int i = 0; i < curr_B; i++) {
				int v = list_B.front();
				list_B.pop_front();
				if (status_B[v] == SUSPENDED || status_B[v] == ADOPTED)
					continue;

				// B adoption test for v
				if (status_A[v] != ADOPTED) { // not A-adopted
					if (alpha_B[v] <= qbo) {
						status_B[v] = ADOPTED;
						if (status_A[v] == SUSPENDED && alpha_A[v] <= qab) {
							status_A[v] = ADOPTED; // reconsideration for A!
							cov++;
							examine_out_neighbors(v, &list_A, &next_A, status_A);
						}
					} else {
						status_B[v] = SUSPENDED;
					}
					
				} else {
					status_B[v] = (alpha_B[v] <= qba) ? ADOPTED : SUSPENDED; // already A-adopted
				}

				if (status_B[v] == ADOPTED) {
					examine_out_neighbors(v, &list_B, &next_B, status_B);
				} // END-IF
			} // END-FOR

			curr_A = next_A;
			curr_B = next_B;
			next_A = next_B = 0;

		} // END-WHILE
	}

	delete[] status_A;
	delete[] status_B;
	delete[] alpha_A;
	delete[] alpha_B;

	return cov / (double) MC_RUNS;
}


void MonteCarlo::examine_out_neighbors(int v, deque<int> *p_list, int *p_next, int *status)
{
	for (int j = 0; j < graph->outDeg[v]; j++) {
		int w = graph->gO[v][j];
		
		if (graph->outEdgeStatus[v][j] == LIVE && status[w] != ADOPTED && status[w] != SUSPENDED) {
			p_list->push_back(w);
			(*p_next)++;
		} else if (graph->outEdgeStatus[v][j] == INACTIVE) {
			double coin = (double) rand() / (double) RAND_MAX;
			
			if (coin <= graph->probO[v][j]) {
				graph->outEdgeStatus[v][j] = LIVE;
				if (status[w] != ADOPTED && status[w] != SUSPENDED) {
					p_list->push_back(w);
					(*p_next)++;
				}
			} else {
				graph->outEdgeStatus[v][j] = BLOCKED;
			}
		} // ENDIF
	} // ENDFOR
}


void MonteCarlo::estimate_spread_from_file(string seed_file_name) 
{
	FILE *fin = fopen(seed_file_name.c_str(), "r");
    ASSERT(fin != NULL);

    int *seed_set = new int[k];
    for (int i = 0; i < k; i++) {
    	int j;
    	double t1, t2;
    	fscanf(fin, "%d\t%d\t%lf\t%lf", &j, &(seed_set[i]), &t1, &t2);
    }
    fclose(fin);

    int found = (int) seed_file_name.find_last_of(".");
    string ofile_name = seed_file_name.substr(0, found) + "_TRUE.txt";
    ofstream myfile;
    myfile.open(ofile_name.c_str());
    if (myfile.is_open()) {
    	for (int i = 0; i < k; i++) {
    		if (i == 0 || (i+1) % 10 == 0) {
    			double spread = compute_coverage(seed_set, i+1);
    			myfile << (i+1) << "\t" << spread << endl;
    			cout << (i+1) << "\t" << spread << endl;	
    		}
    	}
    	myfile.close();
    } else {
    	cerr << "Unable to open output file; no output was written in: " << ofile_name << endl;
		exit(1);
    }

    delete[] seed_set;
}


void MonteCarlo::write_seeds_to_file(bool isB) 
{
	ofstream myfile;
	myfile.open(output_file_name.c_str());
	double spread = 0;
	if (myfile.is_open()) {
		for (int i = 0; i < aSeeds.size(); ++i) {
			spread += mg.at(i);
			if (!isB)
				myfile << i+1 << "\t" << aSeeds.at(i) << "\t" << mg.at(i) << "\t" << spread << endl;
			else
				myfile << i+1 << "\t" << bSeeds.at(i) << "\t" << mg.at(i) << "\t" << spread << endl;
		}
		myfile.close();
	} else {
		cerr << "[error] unable to open output file; no output was written" << output_file_name << endl;
		exit(1);
	}
}


