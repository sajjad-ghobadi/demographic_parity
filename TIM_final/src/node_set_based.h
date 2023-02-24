//typedef double (*pf)(int,int);
#include <math.h> 
#include "head.h"
#include <map>
#include <random>
#include <numeric>
#include <algorithm>
#include "gurobi_c++.h"

class node_set_based: public TimGraph
{
    public:
        int i_multi = 1;
        node_set_based(string folder, string graph_file):TimGraph(folder, graph_file){
                }

        map<vector<int>, double> maximin_lp(double eps_node_set){
            // Implement the maximin_lp algorithm in the paper.
            // Parameters
            // ----------
            // eps_node_set : float
            //     violation in fairness
            // Returns
            // -------
            // dict
            //     Dictionary with keys sets or the returned solution by function set_based_maximin
            //     and values probabilities.     

            map<vector<int>, double> Pp;
            try{
                GRBEnv env = GRBEnv(true);
                env.set("LogFile", "set_dp3.log");
                env.start();

                GRBModel model = GRBModel(env);

                map<vector<int>, double>Pp; 
                int S_size = SS.size();
                GRBVar p[S_size]; 
                for (int i = 0; i < SS.size(); i++){
                    ostringstream vname ;
                    vname << "p" << i;
                    p[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, vname.str());
                }

                GRBVar gamma = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "gamma");
                model.update();
                GRBLinExpr sum1;
                for (auto C:communities){
                    sum1 = 0.0;
                    for (int i = 0; i < SS.size(); i++){
                        sum1 += p[i] * sigma_CS[C][i];  
                    }
                    model.addConstr(sum1 <= gamma + eps_node_set);
                    model.addConstr(sum1 >= gamma - eps_node_set);
                }

                GRBLinExpr sum2 = 0.0;
                for (int i = 0; i < SS.size(); i++){
                    sum2 += p[i];
                }
                model.addConstr(sum2 == 1);
                GRBLinExpr sum4 = 0.0;
                for (int i = 0; i < SS.size(); i++){
                    sum4 += p[i] * SIGMA[i];
                }
                model.update();
                model.setObjective(sum4, GRB_MAXIMIZE);
                model.optimize();
                for (int i = 0; i < SS.size(); i++){
                    if (p[i].get(GRB_DoubleAttr_X) > 0.00000001){       
                        Pp[SS[i]] = p[i].get(GRB_DoubleAttr_X);
                    }
                }
                return Pp; 

            } catch(GRBException e) {
                cout << "Error code = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
            } catch(...) {
                cout << "Exception during optimization" << endl;
            }

        }

        map<vector<int>, double> greedy_group_lp(double eps_node_set){
            // Implement the greedy_group_lp algorithm in the paper.
            // Parameters
            // ----------
            // eps_node_set : float
            //     violation in fairness
            // Returns
            // -------
            // dict
            //     Dictionary with keys sets and values probabilities.   
            map<vector<int>, double> Pp;
            try{
                GRBEnv env = GRBEnv(true);
                env.set("LogFile", "set_dp.log");
                env.start();
                GRBModel model = GRBModel(env);

                map<vector<int>, double>Pp;
                int S_size = SS.size();
                GRBVar p[S_size]; //
                for (int i = 0; i < SS.size(); i++){
                    ostringstream vname ;
                    vname << "p" << i;
                    p[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, vname.str());
                }

                GRBVar gamma = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "gamma");
                model.update();
                GRBLinExpr sum1;
                for (auto C:communities){
                    sum1 = 0.0;
                    for (int i = 0; i < SS.size(); i++){
                        sum1 += p[i] * sigma_CS[C][i];  
                    }
                    model.addConstr(sum1 <= gamma + eps_node_set);
                    model.addConstr(sum1 >= gamma - eps_node_set);
                }

                GRBLinExpr sum2 = 0.0;
                for (int i = 0; i < SS.size(); i++){
                    sum2 += p[i];
                }
                model.addConstr(sum2 == 1);
                GRBQuadExpr sum3 = 0.0;
                for (int i = 0; i < SS.size(); i++){
                    sum3 += p[i] * SS[i].size();
                }
                model.addQConstr(sum3 <= k);

                GRBLinExpr sum4 = 0.0;
                for (int i = 0; i < SS.size(); i++){
                    sum4 += p[i]* SIGMA[i];
                }
                model.update();
                model.setObjective(sum4, GRB_MAXIMIZE);
                model.optimize();

                for (int i = 0; i < SS.size(); i++){
                    if (p[i].get(GRB_DoubleAttr_X) > 0.00000001){       
                        Pp[SS[i]] = p[i].get(GRB_DoubleAttr_X);
                    }
                }
                return Pp; 

            } catch(GRBException e) {
                cout << "Error code = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
            } catch(...) {
                cout << "Exception during optimization" << endl;
            } 
        }

        map<vector<int>, double> independent_lp(double eps_node_set){
            // Implement the independent_lp algorithm in the paper.
            // Parameters
            // ----------
            // eps_node_set : float
            //     violation in fairness
            // Returns
            // -------
            // dict
            //     Dictionary with keys nodes and values probabilities.  
            map<vector<int>, double> Pp;
            
            GRBVar** y = 0;
            try{
                GRBEnv env = GRBEnv(true);
                env.set("LogFile", "node_dp.log");
                env.start();
                GRBModel model = GRBModel(env);
        
                GRBVar x[n];
                for (int i = 0; i < n; i++){
                    ostringstream vname ;
                    vname << "x" << i;
                    x[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, vname.str());
                }

                int T = hyperGT.size();
                y = new GRBVar* [T];
                for (int t = 0; t < T; t++){
                    y[t] = model.addVars(n);
                    for (int j = 0; j < n; j++){
                        ostringstream vname;
                        vname << "y" << t << "," << j;
                        y[t][j].set(GRB_DoubleAttr_LB, 0.0);
                        y[t][j].set(GRB_DoubleAttr_UB, 1.0);
                        y[t][j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                        y[t][j].set(GRB_StringAttr_VarName, vname.str());
                    }
                }

                GRBVar gamma = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "gamma");
                model.update();

                GRBLinExpr sum1 = 0.0;
                for (int i = 0; i < n; i++){
                    sum1 += x[i];
                }
                model.addConstr(sum1 <= k);
                GRBLinExpr sum2;
                for (int i = 0; i < n; i++){
                    for (int t = 0; t < T; t++){
                        sum2 = 0.0;
                        if (find(rootR[i].begin(), rootR[i].end(), t) != rootR[i].end()){
                            for (int node:hyperGT[t]){
                                sum2 += x[node];
                            }
                        }
                        model.addConstr(sum2 >= y[t][i]);
                        model.addConstr(sum2 - eps_node_set <= y[t][i]);                       
                    }
                }
        
                GRBLinExpr sum3;
                for (auto C:communities){
                    sum3 = 0.0;
                    for (int t = 0; t < T; t++){
                        for (int v:communities_node[C]){
                            sum3 += y[t][v] * 1/T * n;
                        }
                    }
                    model.addConstr(sum3 == communities_node[C].size() * gamma);
                }

                GRBLinExpr sum4 = 0.0;
                for (int t = 0; t < T; t++){
                    for (int v = 0; v < n; v++){
                        sum4 += y[t][v] * 1/T * n;
                    }
                }
                
                model.setObjective(sum4, GRB_MAXIMIZE);
                model.optimize();

                for (int i = 0; i < n; i++){
                    if (x[i].get(GRB_DoubleAttr_X) > 0.0){
                        Pp[{i}] = x[i].get(GRB_DoubleAttr_X);
                    }
                }

            } catch(GRBException e) {
                cout << "Error code = " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
            } catch(...) {
                cout << "Exception during optimization" << endl;
            }
            return Pp; 
        }

        double sigma_p(map<vector<int>, double> p){
            // Compute sigma(p).
            // Parameters
            // ----------
            // p : dict
            //     Dictionary with keys sets and values probabilities.
            // Returns
            // -------
            // spread_p : float
            //     the expected spread when sampling sets according to p.   
            
            double spread_p = 0.0;
            for(auto item:p){
               spread_p += sigma_S(item.first) * p[item.first];
            }
            return spread_p;
        }

        map<string, double> sigma_C_p(map<vector<int>, double> p){
            // Compute the average probability that communities are reached according to p.
            // Parameters
            // ----------
            // p : dict
            //     Dictionary with keys sets and values probabilities.
            // Returns
            // -------
            // comm_prob : flot
            //     Dictionary with keys communities and values average probabilities of being reaced
            //     of nodes in communities.

           map<vector<int>, map<string, double>> comm_probs_for_sets;
           for(auto item:p){
               comm_probs_for_sets[item.first] = sigma_C(item.first);
           }
           map<string, double> comm_prob;
           for(auto C:communities_node){
               for(auto item:p){
                   comm_prob[C.first] += p[item.first] * comm_probs_for_sets[item.first][C.first];
               }
           }
           return comm_prob;
        }

        map<vector<int>, double> set_based_maximin(double epsilon, double eps=0.1, string set_or_node = "set"){
            // Call multi_weight for the set-based problem of Becker et al (2022).
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation).
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the set-based problem of Becker et al.

            return multi_weight(epsilon, eps, set_or_node = "set");
        }

        map<vector<int>, double> multi_weight(double epsilon, double eps, string set_or_node = "set"){
            // Implement multi_weight for the set-based algorithm of Becker et al (2022).
            // Parameters
            // ----------
            // epsilon : float
            //      epsilon of (eps, delta)-approximation (for tim implementation)     
            // eps : float
            //      eps of multi_weight.
            // set_or_node: str
            //      string showing that multi_weight should be executed for the set-based problem of Becker et al.
            // Returns
            // -------
            // p : dict
            //      Dictionary with keys sets and values probabilities.
            
            map<string, double>z;
            map<string, double>s;

            for(auto C:communities){
                z[C] = 1.0;
                s[C] = 0.0;
            }

            map<vector<int>, double>p;
            i_multi = 1;
            double primal = -INFINITY;
            double dual = INFINITY;

            while(true)
            {
                double sum_z = 0;
                for(auto item: z){
                    sum_z += item.second;
                }

                double sum_weights = 0.0;
                for(int index=0; index<nodes.size(); index++){
                    double temp = 0.0;
                    for(auto c: node_communities[nodes[index]]){
                        temp += (double) (z[c] / communities_node[c].size());
                    }
                    weights[index] = (double) temp;
                    sum_weights += (double)weights[index];
                }
 
                for(int v:nodes){
                    weights[v] *= (double) (n/sum_weights); //
                }
                EstimateOPT("maximin", epsilon, true);
                vector<int>oracle_solution = seedSet;
                double oracle_value = spread;

                if((double)(oracle_value * (sum_weights / nodes.size()) / sum_z) > 1.0000001){ 
                    cout << "_____________________________"<<endl;
                    for(auto item: z){
                        cout<< item.second<<endl;
                    }
                    cout<< "oracle_value, sum(z.values())" <<endl;
                    cout << oracle_value<< "," << sum_z<< endl;
                    cout << (double)(oracle_value * (sum_weights / nodes.size()) / sum_z) <<endl;
                    cout << "_____________________________"<<endl;
                }
                ASSERT((double)(oracle_value * (sum_weights/nodes.size()) / sum_z) <= 1.0000001); 

                if(set_or_node == "set"){
                    vector<int>difference;
                    sort(oracle_solution.begin(), oracle_solution.end());
                    int ii = 0;
                    for(auto item:p){
                        if(oracle_solution == item.first)
                        {
                            p[oracle_solution] += 1;
                            ii += 1;
                            break;
                        }
                    }
                    if (ii == 0)
                    {
                        p[oracle_solution] = 1;
                    }
                }
                else if(set_or_node == "node"){
                    int jj = 0;
                    for(int v:oracle_solution){
                        jj = 0;
                        for(auto item:p){
                            if(find(item.first.begin(), item.first.end(), v) != item.first.end()){
                                p[item.first] += 1;
                                jj += 1;
                                break;
                            }
                        }
                        if(jj == 0){
                            p[{v}] = 1;
                        }
                    }
                }
                else{
                    cout <<"Error: Unknown option."<< set_or_node <<endl;
                }

                dual = min((double)dual, (double)oracle_value * (sum_weights / nodes.size()) / sum_z);

                map<string, double> community_prob;
                community_prob = sigma_C(oracle_solution);

                for (auto C:communities){
                    z[C] *= (double)(1 - (double)(eps * community_prob[C]));
                }

                for (auto C:communities){
                    s[C] = (double)(i_multi - 1) / i_multi * s[C] + (double)1 / i_multi * community_prob[C];
                }

                auto min_s = *min_element(s.begin(), s.end(),
                    [](decltype(s)::value_type& l, decltype(s)::value_type& r) -> bool { return l.second < r.second; });
                primal = min_s.second;
                if(primal >= (1 - eps) * dual){
                    break;
                }
                i_multi += 1;
            }

            for(auto item:p){
                p[item.first] = (double)item.second/i_multi;
            }
            return p; 
        }
};
