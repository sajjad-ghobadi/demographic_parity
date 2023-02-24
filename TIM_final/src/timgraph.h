class TimGraph: public InfGraph
{
    public:
        TimGraph(string folder, string graph_file):InfGraph(folder, graph_file ){
        }
        double MgT(int u){
            //static int i=0;
            //i++;
            //TRACE_LINE("mgt", i);
            ASSERT(u>=0);
            //ASSERT(u<n);
            return (double)BuildHypergraphNode(u, 0, false);
        }
    
        double algo2(){
            double lb=1/2.0;
            double c=0;
            random_device rd;
            mt19937 generator{rd()}; 
            discrete_distribution<> distribution(weights.begin(),weights.end());
            while(true){
                int loop= (6 * log(nodes.size())  +  6 * log(log(nodes.size())/ log(2)) )* 1/lb;  
                c=0;
                IF_TRACE(int64 now=rdtsc());
                double sumMgTu=0;
                for(int i=0; i<loop; i++){               
                    int u = distribution(generator);        
                    double MgTu=MgT(index_node[u]);
                    double pu=MgTu/m;
                    sumMgTu+=MgTu;
                    c+=1-pow((1-pu), k);
                }
                c/=loop;    
                if(c>lb) break;
                lb /= 2;
            }
            return c * nodes.size();      
        }
        double KptEstimation()
        {
            Timer t(1, "step1");
            double ept=algo2();
            ept/=2;
            return ept;
        }
        void RefindKPT(string maxmn_dp3, double epsilon, double ept){
            Timer t(2, "step2" );
            ASSERT(ept > 0);
            R = (2 + epsilon) * ( nodes.size() * log(nodes.size()) ) / ( epsilon * epsilon * ept); 
            if (setting == "node_dp" or setting == "set_dp2" or setting == "tim" or setting == "greedy"){
                 BuildHypergraphR(R, "kpt_star");
            }
            else if (maxmn_dp3 == "maximin"){
                BuildHypergraphR(R, "maximin");
            }
            else if (maxmn_dp3 == "dp3"){
                BuildHypergraphR(R, "kpt_star");
            }
            else {
                cout <<"Error: Unknown option (maxmn_dp3)."<< endl;
                ASSERT(0); 
            }
              
        }
        double logcnk(int n, int k){
            double ans=0;
            for(int i=n-k+1; i<=n; i++){
                ans+=log(i);
            }
            for(int i=1; i<=k; i++){
                ans-=log(i);
            }
            return ans;
        }
        void NodeSelection(string maximin_dp3,double epsilon, double opt, bool build_SeedSet){
            Timer t(3, "step3");
            ASSERT(opt > 0);
            R = 0.0;      
            if(setting == "node_dp"){
                R = 1000;
            }else{
                R = (8+2 * epsilon) * ( log(nodes.size()) + log(2) +  nodes.size() * logcnk(nodes.size(), k)) / ( epsilon * epsilon * opt); //
            }
            if(setting == "set_dp2"){          
                BuildHypergraphR(R, "set_dp2");
            }
            else if (setting == "tim"){
                BuildHypergraphR(R, "tim");
            }
            else if (setting == "greedy"){
                BuildHypergraphR(R, "greedy");
            }
            else if (maximin_dp3 == "maximin"){
                BuildHypergraphR(R, "maximin");
            }
            else if (maximin_dp3 == "dp3"){
                BuildHypergraphR(R, "dp3");
            }
            else if(setting == "node_dp"){
                BuildHypergraphR(R, "node_dp");
            }
            else {
                cout <<"Error: Unknown option (maximin_dp3)."<< endl;
                ASSERT(0); 
            }
            
            if(build_SeedSet == true){
                BuildSeedSet();
            }
        }
        void EstimateOPT(string maximin_dp3, double epsilon, bool build_SeedSet){
            Timer t(100,"EstimateOPT");

            if (setting == "node_dp"){
                double kpt_plus= 1.0;
                //Node Selection
                NodeSelection(maximin_dp3, epsilon, kpt_plus, build_SeedSet);
                disp_mem_usage(""); 
            }else{
                // KPT estimation
                double kpt_star;
                kpt_star=KptEstimation();
                // Refine KPT
                double eps_prime;
                eps_prime=5*pow(sqr(epsilon)/(k+1), 1.0/3.0);
                RefindKPT(maximin_dp3, eps_prime, kpt_star);
                
                double kpt = BuildSeedSet();
                kpt/=1+eps_prime;
                double kpt_plus = max(kpt, kpt_star);

                //Node Selection
                NodeSelection(maximin_dp3, epsilon, kpt_plus, build_SeedSet);
                disp_mem_usage(""); 
            }
        }
};