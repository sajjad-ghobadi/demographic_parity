class VV
{
    public:
        vector<int> head;
        vector<int> next;
        vector<int> data;
        vector<int> vsize;
        void clear()
        {
            head.clear();
            next.clear();
            data.clear();
            vsize.clear();
        }
        // trick for not change code
        void push_back(vector<int> x){
            ASSERT(x.size()==0);
            addVector();
        }
        void addVector()
        {
            head.push_back(-1);
            vsize.push_back(0);
        }
        int size(int t){
            return vsize[t];
        }
        //vv[a].push_back(b)
        void addElement( int a, int b)
        {
            //a.push_back(b);
            vsize[a]++;
            data.push_back(b);
            next.push_back(head[a]);
            head[a]=next.size()-1;

        }
};

#include <queue>	
#include <utility> 

struct CompareBySecond {
	bool operator()(pair<int, int> a, pair<int, int> b)
	{
		return a.second < b.second;
	}
};

class InfGraph:public Graph
{
    public:
        set<int> seedSet_final;
        int idx = 0.0;
        int64 R = 0.0;
        double spread = 0.0;
        vector<vector<int>> SS;
        vector<double> SIGMA;
        map<string, map<int, double>> sigma_CS;
        vector<vector<int>> hyperGT;  // RR sets
        map<int, vector<int>> hyperG;  
        map<int, vector<int>> rootR; 
        
        InfGraph(string folder, string graph_file):Graph(folder, graph_file){
            hyperG.clear();

            for(int i=0; i<n; i++){
                hyperG.insert(pair<int,vector<int>>(i, vector<int>())); 
            }
            for(int i=0; i<12; i++)
                sfmt_init_gen_rand(&sfmtSeed, i+1234); 
        }

        enum ProbModel {TR, WC, TR001};
        ProbModel probModel;


        void hyper_GT_G(int64 R, vector<int>nodes_or_comm_nodes){
            // Create RR sets (hyperGT), hyperG and rootR based on nodes_or_comm_nodes
            // Parameters
            // ----------
            // nodes_or_comm_nodes : list of nodes
            // R : int64
            //    number showing the number of RR sets

            hyperG.clear();
            rootR.clear();   
            for(int i:nodes){
                hyperG.insert(pair<int,vector<int>>(i, vector<int>()));    
                rootR.insert(pair<int,vector<int>>(i, vector<int>()));  
            }
                
            int R_per_node = 0.0;
            R_per_node = (int) (ceil((double) R/nodes_or_comm_nodes.size()));  
            int up_R = R_per_node * nodes_or_comm_nodes.size();
            hyperId = up_R;

            hyperGT.clear();
            while((int)hyperGT.size() < up_R)
                hyperGT.push_back(vector<int>());

            int index = 0.0;
            for (int v:nodes_or_comm_nodes){
                for (int t = 0; t < R_per_node; t++){
                    rootR[v].push_back(index);        
                    BuildHypergraphNode(v, index, true);
                    index += 1;
                }
            }

            int totAddedElement=0; 
            for(int i=0; i<up_R; i++){     
                for(int t:hyperGT[i]){
                    hyperG[t].push_back(i);
                    totAddedElement++;
                }
            }
        }

        void BuildHypergraphR(int64 R, string kpt_or_others){
            // Execution of the experiments
            // Parameters
            // ----------
            // R : int64
            //    number showing the number of RR sets
            // kpt_or_others : str
            //                string based on the setting

            if (kpt_or_others == "kpt_star" or kpt_or_others == "node_dp" or kpt_or_others == "tim"){
                // for independent_lp and greedy_im
                k_dp = k;
                hyper_GT_G(R, nodes);
            }
            else if (kpt_or_others == "maximin"){ // for multi_weight of Becker et al 2022.
                k_dp = k;
                hyperG.clear();
                rootR.clear();   
                for(int i=0; i<n; i++){
                    hyperG.insert(pair<int,vector<int>>(i, vector<int>()));    
                    rootR.insert(pair<int,vector<int>>(i, vector<int>()));  
                }
                random_device rd;
                mt19937 generator{rd()}; 
                discrete_distribution<> distribution(weights.begin(),weights.end());

                hyperGT.clear();
                    while((int)hyperGT.size() < R)
                        hyperGT.push_back(vector<int>());

                hyperId = R; 
                    
                for(int i=0; i<R; i++){
                    int number = distribution(generator);
                    rootR[index_node[number]].push_back(i); 
                    BuildHypergraphNode(index_node[number], i, true); 
                }
                int totAddedElement=0; 
                for(int i=0; i<R; i++){     
                    for(int t:hyperGT[i]){
                        hyperG[t].push_back(i);
                        totAddedElement++;
                    }
                }
                ASSERT(hyperId == R);    
            }
            else if (kpt_or_others == "set_dp2"){ // for greedy_group_lp
                idx = 0.0;
                for (auto C:communities){
                    hyper_GT_G(R, communities_node[C]);
                    k_dp = k;
                    vector<int> set_solution;
                    if (k >= communities_node[C].size()){
                        for (int i:communities_node[C]){
                            set_solution.push_back(i);
                        }
                    }else{
                        BuildSeedSet();
                        set_solution = seedSet;
                    }
                    SS.push_back(set_solution);
                    idx += 1;
                }

                hyper_GT_G(R, nodes);          
                map<string, double> pro_C;
                for (int j = 0; j < SS.size(); j++){
                    pro_C.clear();
                    double set_value = 0.0;
                    set_value = sigma_S(SS[j]);
                    SIGMA.push_back(set_value);
                    pro_C = sigma_C(SS[j]);
                    for(auto C:communities){
                        sigma_CS[C][j] = pro_C[C];
                    }
                }

                k_dp = 2*k;
                BuildSeedSet(); 

                for (int j=0 ; j < Ss.size(); j++){
                    SS.push_back(Ss[j]);
                    SIGMA.push_back(S_sigma[j]);
                }
                for(auto C:communities){
                    sigma_CS[C][idx] = 0.0;
                }
                idx += 1;

                map<string, double> prob_C;
                for (auto ss:Ss){
                     if (ss.first != 0.0){
                        prob_C.clear();
                        prob_C = sigma_C(ss.second);
                        for(auto C:communities){
                            sigma_CS[C][idx] = prob_C[C];
                        }
                        idx += 1;
                        if (ss.first == k){
                            auto min_C = *min_element(prob_C.begin(), prob_C.end(),
                                [](decltype(prob_C)::value_type& l, decltype(prob_C)::value_type& r) -> bool { return l.second < r.second; });

                            auto max_C = *max_element(prob_C.begin(), prob_C.end(),
                                [](decltype(prob_C)::value_type& l, decltype(prob_C)::value_type& r) -> bool { return l.second < r.second; });
                            
                            double diff = max_C.second - min_C.second;
                            vector<double> range = {16.0, 8.0, 4.0, 2.0};
                            for(auto x:range){
                                violations.push_back((double) diff / x);
                            } 
                        }   
                    }
                } 
            }
            
            else if (kpt_or_others == "dp3"){ // for maximin_lp
                idx = 0.0;
                for (auto C:communities){
                    hyper_GT_G(R, communities_node[C]);
                    vector<int> set_solution;
                    if (k >= communities_node[C].size()){
                        for (int i:communities_node[C]){
                            set_solution.push_back(i);
                        }
                    }else{
                        BuildSeedSet();
                        set_solution = seedSet;
                    }
                    SS.push_back(set_solution);
                    idx += 1;
                }
                hyper_GT_G(R, nodes);          
                map<string, double> pro_C;
                for (int j = 0; j < SS.size(); j++){
                    pro_C.clear();
                    double set_value = 0.0;
                    set_value = sigma_S(SS[j]);
                    SIGMA.push_back(set_value);
                    pro_C = sigma_C(SS[j]);
                    for(auto C:communities){
                        sigma_CS[C][j] = pro_C[C];
                    }
                }
            } 

            else if (kpt_or_others == "greedy"){ // for greedy_prop
                int sum_k = 0.0;
                map<string, int> comm_k;

                double summ_comms_size = 0.0;
                for (auto C:communities){
                    summ_comms_size += communities_node[C].size();
                }

                for (auto C:communities){
                    comm_k[C] = (int) (floor((double) k * ((double)communities_node[C].size()/summ_comms_size)));
                    sum_k += comm_k[C];
                }
                
                if (sum_k < k){
                    for (auto C:communities){
                        comm_k[C] += 1;
                        sum_k += 1;
                        if (sum_k == k){
                            break;
                        }
                    }
                }

                for (auto C:communities){
                    hyper_GT_G(R, communities_node[C]);
                    k_dp = comm_k[C];
                    BuildSeedSet();
                    for (int v:seedSet){
                        seedSet_final.insert(v);
                    }    
                }

                if(seedSet_final.size() < sum_k){
                    hyper_GT_G(R, nodes);
                    k_dp = sum_k;
                    BuildSeedSet();
                    for(int v:seedSet){
                        seedSet_final.insert(v);
                        if (seedSet_final.size() == k){
                            break;
                        }
                    }  
                }  
            }
            else {
                cout <<"Error: Unknown option (kpt_or_others)."<< endl;
                ASSERT(0); 
            }
        }

        int BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge){
            int n_visit_edge=1;
            if(addHyperEdge)
            {
                ASSERT((int)hyperGT.size() > hyperiiid);
                hyperGT[hyperiiid].push_back(uStart);
            }

            int n_visit_mark=0;
            q.clear();
            q.push_back(uStart);
            ASSERT(n_visit_mark < n);
            visit_mark[n_visit_mark++]=uStart;
            visit[uStart]=true;
            while(!q.empty()) {
                int expand=q.front();
                q.pop_front();
                if(influModel==IC){
                    int i=expand;
                    for(int j=0; j<(int)gT[i].size(); j++){
                        int v=gT[i][j];
                        n_visit_edge++;
                        double randDouble=double(sfmt_genrand_uint32(&sfmtSeed))/double(RAND_MAX)/2; // 
                        if(randDouble > probT[i][j])
                            continue;
                        if(visit[v])
                            continue;
                        if(!visit[v])
                        {
                            ASSERT(n_visit_mark < n);
                            visit_mark[n_visit_mark++]=v;
                            visit[v]=true;
                        }
                        q.push_back(v);
                        if(addHyperEdge)
                        {
                            ASSERT((int)hyperGT.size() > hyperiiid);
                            hyperGT[hyperiiid].push_back(v);
                        }
                    }
                }
                else if(influModel==LT){
                    if(gT[expand].size()==0)
                        continue;
                    ASSERT(gT[expand].size()>0);
                    n_visit_edge+=gT[expand].size();
                    double randDouble=double(sfmt_genrand_uint32(&sfmtSeed))/double(RAND_MAX)/2;
                    for(int i=0; i<(int)gT[expand].size(); i++){
                        ASSERT( i< (int)probT[expand].size());
                        randDouble -= probT[expand][i];
                        if(randDouble>0)
                            continue;
                        int v=gT[expand][i];

                        if(visit[v])
                            break;
                        if(!visit[v])
                        {
                            visit_mark[n_visit_mark++]=v;
                            visit[v]=true;
                        }
                        q.push_back(v);
                        if(addHyperEdge)
                        {
                            ASSERT((int)hyperGT.size() > hyperiiid);
                            hyperGT[hyperiiid].push_back(v);
                        }
                        break;
                    }
                }
                else
                    ASSERT(false);
            }
            for(int i=0; i<n_visit_mark; i++)
                visit[visit_mark[i]]=false;
            return n_visit_edge;
        }

        //return the number of edges visited
        int64 hyperId = 0;
        deque<int> q;
        sfmt_t sfmtSeed;
		vector<int> seedSet;
        vector<double> sigma;
        map<int, vector<int>> Ss; 
        map<int, double> S_sigma; 
		double BuildSeedSet(){
            // Compute seed set (based on k) and expected spread from seed set
			priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
            map<int, int>coverage;
        
            for (int i=0; i<n; i++){
				pair<int, int>tep(make_pair(i, (int)hyperG[i].size()));
				heap.push(tep);
				coverage[i] = (int)hyperG[i].size();
			}
            if (k_dp == 2*k){
                Ss.clear();
                S_sigma.clear();
                Ss[0.0] = {};          
                S_sigma[0.0] = 0.0; 
            }
			int maxInd;

			long long influence = 0;
			long long numEdge = hyperGT.size();

			// check if an edge is removed
			vector<bool> edgeMark(numEdge, false);
			// check if an node is remained in the heap
            map<int, bool> nodeMark;
            for (int v=0; v<n; v++){
                nodeMark[v] = true;
            }

			seedSet.clear();            
			while ((int)seedSet.size()<k_dp) 
			{
				pair<int, int>ele = heap.top();
				heap.pop();
				if (ele.second > coverage[ele.first])
				{
					ele.second = coverage[ele.first];
					heap.push(ele);
					continue;
				}

				maxInd = ele.first;
				vector<int>e = hyperG[maxInd];  //the edge influence
				influence += coverage[maxInd];
				seedSet.push_back(maxInd);
				nodeMark[maxInd] = false;

                if (k_dp == 2*k){
                    Ss[seedSet.size()] = seedSet;              
                    S_sigma[seedSet.size()] = double(1.0*influence / hyperGT.size()*nodes.size()); //
                    sigma.push_back(1.0*influence / hyperGT.size()*nodes.size());
                }
            
				for (unsigned int j = 0; j < e.size(); ++j){
					if (edgeMark[e[j]])continue;
					vector<int>nList = hyperGT[e[j]];
					for (unsigned int l = 0; l < nList.size(); ++l){
						if (nodeMark[nList[l]]){
                            coverage[nList[l]]--;
                        }
					}
					edgeMark[e[j]] = true;
				}
			}

            spread = 1.0*influence / hyperGT.size()*nodes.size();
			return  spread;
		}

        map<int, double>sigma_v(vector<int> S){
            // Compute sigma_v(S).
            // Parameters
            // ----------
            // S : list of nodes
            // Returns
            // -------
            // pro_v : dict
            //       Dictionary with keys node and values probabilities of being reached of nodes from S. 

           sort(S.begin(), S.end());
           map<int, double> pro_v;
           for(int index=0; index<nodes.size(); index++){
               pro_v[nodes[index]] = 0.0;
               for(int r:rootR[nodes[index]]){
                   sort(hyperGT[r].begin(), hyperGT[r].end());
                   vector<int> intersect;   
                    set_intersection(hyperGT[r].begin(), hyperGT[r].end(), S.begin(),
                    S.end(), back_inserter(intersect));
                    if(intersect.size()>0){
                       pro_v[nodes[index]] += 1.0;
                   }
               }
               pro_v[nodes[index]] = pro_v[nodes[index]] * (double) n/(hyperGT.size() * weights[index]);
            }
           return pro_v;
        }
       
        map<string, double>sigma_C(vector<int> S){
            // Compute sigma_C(S).
            // Parameters
            // ----------
            // S : list of nodes
            // Returns
            // -------
            // comm_probs : dict
            //       Dictionary with keys communities and values average probability of being reached of nodes in each community from S. 
           map<int, double>node_probs = sigma_v(S);
           map<string, double>comm_probs;
           for(auto C:communities){
               comm_probs[C] = 0.0;
           }
           for (int v:nodes){
               for (auto C:node_communities[v]){
                       comm_probs[C] += (double)(node_probs[v] / communities_node[C].size());
               }
           }
           return comm_probs;
        }

        double sigma_S(vector<int> S){
            // Compute sigma(S).
            // Parameters
            // ----------
            // p : list of nodes
            // Returns
            // -------
            // infl : float
            //       number showing the expected number of nodes reached from S.         
           map<int, double>node_probs = sigma_v(S);           
           double infl = 0.0;
           for (int v:nodes){
               infl += node_probs[v];
           }
           return infl;
        }
};



