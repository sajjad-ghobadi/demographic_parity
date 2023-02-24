This is the source code used for the experimental part of the paper 
			On the Cost of Demographic Parity in Influence Maximization published at AAAI-23

For execution of the experiments use
    python main_dp.py experiment_type N

where experiment_type is one of[

	    'ba-singletons-0_0.4-200',
            'ba-singletons-0_0.4-100',
            'ba-singletons-0_0.4-50',
            'tsang-region-gender-0-0.4',
            'irvine_0_0.2_bfs_comm_10',
            'email-Eu-core_0_0.2',
            'ca-GrQc_0_0.2-leiden',
            'ca-HepTh_0_0.2_random-n_10',
            'facebook_combined_0_0.2_bfs_comm_2',
            ]
and N is the number of experiments that are supposed to be run in parallel.

For creating the instances use
    python main_generate.py experiment_type


- For creating the networks used by Fish et al. [1] download the data sets into a subfolder 'fish' and put it into the folder 'data_set'.

- For running tests on the networks used by Tsang et al. [2], download the instances into the subfolder 'networks' and put it into the folder 'tsang'.

- The files have the following content:
    - main_dp.py: contains the infrastructure for the execution of the experiments.
    - generation.py: contains the functions that generate the different instances.
    - influence_max.py: contains the necessary implementations of influence maximization functions: computation of the influence sigma of a set (and related functions sigma_v, sigma_C, see the paper), computation of nodes reachable from a given set, generation of live edge graphs, greedy algorithms for influence maximization and maximin
    - maximin_fish.py: contains the implementations of the myopic routine of Fish et al.
    - gurobi_influence_max.py: contains the implementations of the milp of Farnadi et al. [3]
    - main_generate.py: contains the infrastructure for generating the instances used in the paper.
    - print_functions.py: contains some functions for prettier printing of graphs etc.
    - tsang_maximin.py: contains the function used to call the code of Tsang et al. [2]
    - cpp_proxy.py: contains the functions for executing the c++ implementations, the algorithms proposed in the paper and the multiplicative weight routine of Young [4] for the set-based problem of Becker et al.[5].

- Folder 'TIM_final' contains the c++ implementation of the proposed algorithms in the paper. It uses TIM implementation for influence maximization [6].
- The execution generates an output file within the folder 'results_dp' with the name being the experiment_type.


- The experiment_types have the following meaning, see the paper for the interpretation of the community structure types:
    - 'ba-singletons-0_0.4-200':            + Barabasi albert graph (parameter m=2),
                                            + singleton community structure
                                            + edge weights uniformly at random in [0,0.4]
                                            + k = 25
                                            + n = 200

    - 'ba-singletons-0_0.4-100':            + Barabasi albert graph (parameter m=2),
                                            + singleton community structure
                                            + edge weights uniformly at random in [0,0.4]
                                            + k = 25
                                            + n = 100

    - 'ba-singletons-0_0.4-50':             + Barabasi albert graph (parameter m=2),
                                            + singleton community structure
                                            + edge weights uniformly at random in [0,0.4]
                                            + k = 25
                                            + n = 50

    - 'tsang-region-gender-0-0.4':          + instances of Tsang et al. (2019)
                                            + community structure induced by attributes gender and region
                                            + edge weights uniformly at random in [0,0.4]
                                            + k = 25
                                            + n = 500

    - 'arena_0_0.2-random_overlap_10':      + instances used by Fish et al. (2019)
                                            + random_overlap community structure with 10 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 100
                                            + n = 1133                                

    - 'irvine_0_0.2_bfs_comm_10':           + instances used by Fish et al. (2019)
                                            + BFS community structure with 10 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 50
                                            + n = 1899
                                            
    - 'email-Eu-core_0_0.2':                + instances used by Fish et al. (2019)
                                            + community structure induced by departments
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 100
                                            + n = 1005

    - 'ca-GrQc_0_0.2-leiden':               + instances used by Fish et al. (2019)
                                            + leidenalg community structure [7]
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 100
                                            + n = 5242

    - 'ca-HepTh_0_0.2_random-n_10':         + instances used by Fish et al. (2019)
                                            + random community structure with n/10 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 100
                                            + n = 9877

    - 'facebook_combined_0_0.2_bfs_comm_2': + instances used by Fish et al. (2019)
                                            + BFS community structure with 2 communities
                                            + edge weights uniformly at random in [0,0.2]
                                            + k = 50
                                            + n = 4039


[1] Fish, Bashardoust, Boyd, Friedler, Scheidegger, Venkatasubramanian. Gaps in Information Access in Social Networks. WWW 2019.
    http://sorelle.friedler.net/papers/access_networks_www19.pdf
    
[2] Tsang, Wilder, Rice, Tambe, Zick. Group-Fairness in Influence Maximization. IJCAI 2019.
    https://www.ijcai.org/Proceedings/2019/0831.pdf
    
[3] Farnad, Golnoosh, Behrouz Babaki, and Michel Gendreau. "A unifying framework for fairness-aware influence maximization." Companion Proceedings of the Web Conference 2020.

[4] Young, Neal E. "Randomized Rounding Without Solving the Linear Program." SODA. Vol. 95. 1995.

[5] Becker, D’Angelo, Ghobadi, and Gilbert. Fairness in Influence Maximization through Randomization. J. Artif. Intell. Res., 73: 1251–1283. 2022.

[6] Tang, Youze, Xiaokui Xiao, and Yanchen Shi. "Influence maximization: Near-optimal time complexity meets practical efficiency." Proceedings of the 2014 ACM SIGMOD international conference on Management of data. 2014.

[7] Traag, Vincent A., Ludo Waltman, and Nees Jan Van Eck. "From Louvain to Leiden: guaranteeing well-connected communities." Scientific reports 9.1 (2019): 1-12.
