//#define HEAD_TRACE
#define HEAD_INFO

#define HEAD_INFO
//#define HEAD_TRACE
#include "sfmt/SFMT.h"
#include "head.h"
#include "memoryusage.h"
#include "graph.h"
#include "node_set_based.h"
#include <set>
#include <chrono>
#include <iomanip>

void run(node_set_based & x, string dataset, int k, double epsilon, string model, string setting, string file_name, string folder){
    //Infrastructure for the execution of the experiments
    // Parameters
    // ----------
    // dataset : str
    //    path to the dataset to be executed
    // k : int
    //    budget (seed set size)
    // epsilon : float
    //    epsilon of (eps, delta)-approximation (for tim implementation).
    // model : str
    //    underlying diffusion model (IC or LT)
    // setting : str
    //    string showing which function should be executed
    // file_name : str 
    //    name of the file the results should be written
    // folder : str
    //    name of the folder the results (using file_name) should be created

    cout << "dataset:" << dataset << " k:" << k << " epsilon:"<< epsilon << " model:" << model << " setting:" << setting << " file_name:" << file_name << " folder: " << folder << endl;
    x.k=k;
    x.setting = setting;

    if(model=="IC")
        x.setInfuModel(InfGraph::IC);
    else if(model=="LT")
        x.setInfuModel(InfGraph::LT);
    else
        ASSERT(false);

    map<vector<int>, double>p_x;
    vector<int> S;
    set<int> S_greedy;

    auto start = chrono::high_resolution_clock::now();
    ios_base::sync_with_stdio(false);

    int n_time = 0.0;
    double time2 = 0.0;
    double timee = 0.0;

    vector<string> omega;
    if (setting == "node_dp" or setting == "set_dp2"){ // for independent_lp and greedy_group_lp
        if(setting == "node_dp"){
            omega = {"0","1/4", "1/3", "1/2"};
            x.violations = {0, 0.25, 0.33, 0.5};
        }
        if(setting == "set_dp2"){
            omega = {"0", "x/16", "x/8", "x/4", "x/2"};
            x.violations.push_back(0.0); 
        }
        x.EstimateOPT("T", epsilon, true);
        auto end1 = chrono::high_resolution_clock::now();
        double time_taken1 = chrono::duration_cast<chrono::nanoseconds>(end1 - start).count();
        time_taken1 *= 1e-9;
        int count = 0.0;
        for (int i = 0; i < x.violations.size(); i++){
            p_x.clear();
            if (setting == "node_dp"){
                p_x = x.independent_lp(x.violations[i]);
            }
            if(setting == "set_dp2"){
                p_x = x.greedy_group_lp(x.violations[i]);
            }
            
            auto end = chrono::high_resolution_clock::now();
            double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            time_taken *= 1e-9;
            if (n_time == 0.0){
                timee = time_taken;
                time2 = time_taken;
                n_time += 1;
            }else{
                time2 = time_taken - timee + time_taken1;
                timee = time_taken;
            }
            count += 1;
            double rr = (double) rand() / (RAND_MAX);
            string file_name_cpp = to_string(x.violations[i]) + to_string(time_taken)+ "_" + to_string(count) + "_" + to_string(x.k) + to_string(rr) + file_name + ".txt";
            ofstream out_file(folder+file_name_cpp);
            out_file << time2 << endl;
            out_file << omega[i] << endl;
            out_file << x.violations[i] << endl;
            for(auto S_pro:p_x){
                for (int s:S_pro.first){
                    out_file << s << " ";
                }
            out_file << fixed << setprecision(20) << S_pro.second << endl;
            }
            out_file.close();
        }
    }
    else if (setting == "set_dp3"){ // for maximin_lp
        map<vector<int>, double> p_maximin;
        vector<double> vio_dp3;
        vector<string> omega = {"0", "x/16", "x/8", "x/4", "x/2"};
        vio_dp3.push_back(0.0); 
        x.EstimateOPT("dp3", epsilon, true);
        p_maximin = x.set_based_maximin(epsilon);

        x.hyper_GT_G(x.R, x.nodes);  
           
        for(int v:x.nodes){
            x.weights[v] = 1.0; 
        }
   
        map<string, double> C_probs;
        C_probs = x.sigma_C_p(p_maximin);
        for(auto C:x.communities){
            x.sigma_CS[C][x.idx] = C_probs[C];
        }
        x.idx += 1;             

        double infl_p = x.sigma_p(p_maximin);
        x.SS.push_back({99999999});           
        x.SIGMA.push_back(infl_p);

        x.BuildSeedSet();
        map<string, double> c_pros; 
        c_pros = x.sigma_C(x.seedSet);   
        
        auto min_s = *min_element(c_pros.begin(), c_pros.end(),
            [](decltype(c_pros)::value_type& l, decltype(c_pros)::value_type& r) -> bool { return l.second < r.second; });
        double min_C = min_s.second;

        auto max_s = *max_element(c_pros.begin(), c_pros.end(),
            [](decltype(c_pros)::value_type& l, decltype(c_pros)::value_type& r) -> bool { return l.second < r.second; });

        double max_C  = max_s.second;
        double diff = max_C - min_C;

        vector<double> range = {16.0, 8.0, 4.0, 2.0};
        for(auto x:range){
            vio_dp3.push_back((double) diff / x);
        }

        x.SS.push_back({});          
        x.SIGMA.push_back(0.0);
        for (auto C:x.communities){
            x.sigma_CS[C][x.idx] = 0.0;
        }
        
        auto end1 = chrono::high_resolution_clock::now();
        double time_taken1 = chrono::duration_cast<chrono::nanoseconds>(end1 - start).count();
        time_taken1 *= 1e-9;
        int count = 0.0;
        for (int i = 0; i < vio_dp3.size(); i++){
            p_x.clear();
            p_x = x.maximin_lp(vio_dp3[i]);

            auto end = chrono::high_resolution_clock::now();
            double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            time_taken *= 1e-9;
            if (n_time == 0.0){
                timee = time_taken;
                time2 = time_taken;
                n_time += 1;
            }else{
                time2 = time_taken - timee + time_taken1;
                timee = time_taken;
            }
            count += 1;
            double rr = (double) rand() / (RAND_MAX);
            string file_name_cpp = to_string(vio_dp3[i]) + to_string(time_taken)+ "_" + to_string(count) + "_" + to_string(x.k) + to_string(rr) + file_name + ".txt";
            ofstream out_file(folder+file_name_cpp);
            out_file << time2 << endl;
            out_file << omega[i] << endl;
            out_file << vio_dp3[i] << endl;

            for(auto S_pro:p_x){
                int a = 99999999;
                if (S_pro.first.size() == 0.0){
                    out_file << fixed << setprecision(20) << S_pro.second << endl;
                }
                for(int s:S_pro.first){
                    if (s!= 99999999){
                        a = s;
                        out_file << s << " ";
                    }
                }
                if (a != 99999999){
                 out_file << fixed << setprecision(20) << S_pro.second << endl;
                }
            }
            
            for(auto S_pro:p_x){
                for(int s:S_pro.first){
                    if(s == 99999999){
                        // p_maximin
                        out_file << "maximin" << endl;
                        out_file << fixed << setprecision(20) << S_pro.second << endl; 
                        for (auto s_pro_maximin:p_maximin){
                            for (int v:s_pro_maximin.first){
                                out_file << v << " ";
                            }
                        out_file << fixed << setprecision(20) << s_pro_maximin.second << endl;  
                        } 
                    }
                }
            }
            out_file.close();
        }
    }
    else if(setting == "set"){ // for multi_weight of Becker et al 2022
        p_x = x.set_based_maximin(epsilon);    
    }
    else if (setting == "tim"){ // for greedy_im
        x.EstimateOPT("T", epsilon, true);
        S = x.seedSet;
    }
    else if (setting == "greedy"){ // for greedy_prop
        x.EstimateOPT("T",epsilon, true);
        S_greedy = x.seedSet_final;
    }
    else{
        cout <<"Error: Unknown option (setting)."<< endl;
        ASSERT(0); 
    }

    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    if (setting == "set" or setting == "node"){
        ofstream out_file;
        out_file.open(file_name);
        out_file << time_taken << endl;
        for(auto S_pro:p_x){
            for (int s:S_pro.first){
                out_file << s << " ";
            }
            out_file << fixed << setprecision(20) << S_pro.second << endl;
        }
    }
    if (setting == "tim"){
        ofstream out_file;
        out_file.open(file_name);
        out_file << time_taken << endl;
        int i = 0;
        for(int s:S){
            i += 1;
            if(i < S.size()){
                out_file << s << " ";
            }
            else{
                out_file << s << endl;
            }
        }
    }
    if (setting == "greedy"){
        ofstream out_file;
        out_file.open(file_name);
        out_file << time_taken << endl;
        int i = 0;
        for(int s:S_greedy){
            i += 1;
            if(i < S_greedy.size()){
                out_file << s << " ";
            }
            else{
                out_file << s << endl;
            }
        }
    }
       
    Counter::show();
}
void parseArg(int argn, char ** argv)
{
    string folder;
    string file_name=""; 
    string setting=""; 
    string dataset="";

    double epsilon=0;
    string model="";
    int k=0;

    for(int i=0; i<argn; i++)
    {
        if(argv[i]==string("-folder")) folder=string(argv[i+1]);
        if(argv[i]==string("-file_name")) file_name=string(argv[i+1]);
        if(argv[i]==string("-setting")) setting=string(argv[i+1]);
        if(argv[i]==string("-dataset")) dataset=string(argv[i+1])+"/";
        if(argv[i]==string("-epsilon")) epsilon=atof(argv[i+1]);
        if(argv[i]==string("-k")) k=atoi(argv[i+1]);
        if(argv[i]==string("-model")) {
            if(argv[i+1]==string("LT"))
            {
                model=argv[i+1];
            }
            else if(argv[i+1]==string("IC"))
            {
                model=argv[i+1];
            }
            else
                ExitMessage("model should be IC or LT");
        }
    }
    if (dataset=="")
        ExitMessage("argument dataset missing");
    if (k==0)
        ExitMessage("argument k missing");
    if (epsilon==0)
        ExitMessage("argument epsilon missing");
    if (model=="")
        ExitMessage("argument model missing");
    if (setting=="")                             
        ExitMessage("argument setting missing");
    if (file_name=="")                             
        ExitMessage("argument file_name missing");  

    string graph_file;
    if(model=="IC")
        graph_file=dataset + "graph_ic.txt"; 
    else if(model=="LT")
        graph_file=dataset + "graph_lt.inf";

    node_set_based x(dataset, graph_file);  
   
    run(x, dataset, k, epsilon, model, setting, file_name, folder);
}

int main(int argn, char ** argv)
{
    OutputInfo info(argn, argv);
    parseArg( argn, argv );
}
