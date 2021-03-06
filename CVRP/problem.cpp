#include<cassert>
#include<iostream>
#include<unordered_map>
#include<vector>
#include<thread>

#include<problem.hpp>
#include<ilcplex/cplexx.h>


/**********************/
Problem::Problem() = default;
Problem::Problem(int n) : costs(n, std::vector<double>(n, 0.0)), demands(n, 0.0) {}

/**********************/
lazy_constrain_info::lazy_constrain_info(CPXLPptr lp, Problem* prob) : lp(lp), prob(prob) {}

/**********************/
int initialize_structures(CPXENVptr& env, CPXLPptr& lp) {

	int status = 0;
	env = CPXXopenCPLEX(&status);

	//std::cout<<env<<" "<<status<<std::endl;
	if(env == nullptr) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXXgeterrorstring(env, status, errmsg);
		std::cerr << "Could not open CPLEX env" << std::endl << errmsg << std::endl;
		return 1;
	}

	/* Create the problem. */
	std::string lpname = "1e.lp";
	lp = CPXXcreateprob(env, &status, lpname.c_str());

	if(lp == nullptr) {
		std::cerr << "Failed to create LP" << std::endl;
		return 1;
	}

	return 0;
}

void free_structures(CPXENVptr env, CPXLPptr lp) {
	int status = 0;

	if(lp != nullptr) {
		status = CPXXfreeprob(env, &lp);
		if(status != 0) {
			std::cerr << "Problem feeing lp - status: " << status << std::endl;
		}
	}

	if(env != nullptr) {
		status = CPXXcloseCPLEX(&env);
		if(status != 0) {
			std::cerr << "Problem closing CPLEX environment" << std::endl;
			char errmsg[CPXMESSAGEBUFSIZE];
			CPXXgeterrorstring(env, status, errmsg);
			std::cerr << errmsg << std::endl;
		}
	}

}

int set_parameters(CPXENVptr env) {
	int status;

	status = CPXXsetintparam(env, CPXPARAM_Read_DataCheck, CPX_DATACHECK_ASSIST);
	if(status != 0) { return status; }

	/****** Output ******/
	status = CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
	if(status != 0) { return status; }

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Display, 5);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Interval, 1 );
	if(status != 0) { return(status); }*/
	/**********************/

	/****** Presolve ******/
	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
	if(status != 0) { return status; }

	status = CPXXsetintparam(env, CPXPARAM_Preprocessing_Presolve, CPX_ON);
	if(status != 0) { return status; }

	/*status = CPXXsetintparam(env, CPXPARAM_Preprocessing_Reduce, CPX_PREREDUCE_NOPRIMALORDUAL);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_Preprocessing_Linear, 0);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_Preprocessing_Reduce, CPX_PREREDUCE_PRIMALONLY);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_Preprocessing_RepeatPresolve, 0);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Probe, -1);
	if(status != 0) { return(status); }*/
	/**********************/

	/** CPLEX heuristics **/
	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_HeuristicFreq, -1);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_RINSHeur, -1);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_FPHeur, -1);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_LBHeur, CPX_OFF);
	if(status != 0) { return(status); }
	/**********************/

	/***** CPLEX cuts *****/
	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Limits_CutPasses, 100);
	if(status != 0) { return(status); }*/

	status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_Cliques, 0);
	if(status != 0) { return (status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_Covers, 0);
	if(status != 0) { return (status); }

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_Disjunctive, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_FlowCovers, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_PathCut, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_Gomory, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_GUBCovers, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_MCFCut, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_Implied, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_MIRCut, -1);
	if(status != 0) { return (status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Cuts_ZeroHalfCut, -1);
	if(status != 0) { return (status); }*/
	/**********************/

	/** Execution decisions **/
	status = CPXXsetdblparam(env, CPXPARAM_TimeLimit, 3600.0);
	if(status != 0) { return status; }

	/*status = CPXXsetintparam(env, CPXPARAM_MIP_SubMIP_NodeLimit, 1);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPXPARAM_Parallel, 1);
	if(status != 0) { return(status); }*/

	status = CPXXsetintparam(env, CPXPARAM_Threads, std::thread::hardware_concurrency());
	if(status != 0) { return status; }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_DFS);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_BESTBOUND);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_BESTEST);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_BESTEST_ALT);
	if(status != 0) { return status; }

	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Branch, -1);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Branch, 0);
	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Branch, 1);
	if(status != 0) { return status; }

	//status = CPXXsetintparam(env, CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_OPTIMALITY);
	status = CPXXsetintparam(env, CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_BALANCED);
	//status = CPXXsetintparam(env, CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_FEASIBILITY);
	//status = CPXXsetintparam(env, CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_BESTBOUND);
	//status = CPXXsetintparam(env, CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_HIDDENFEAS);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_AUTOMATIC);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_PRIMAL);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_DUAL);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_NET);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_BARRIER);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_SIFTING);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_StartAlgorithm, CPX_ALG_CONCURRENT);
	if(status != 0) { return(status); }

	/*status = CPXXsetintparam(env, CPXPARAM_Barrier_Crossover, 0);
	if(status != 0) { return(status); }*/

	status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_AUTOMATIC);
	//status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_PRIMAL);
	//status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_DUAL);
	//status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_NET);
	//status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_BARRIER);
	//status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_SIFTING);
	//status = CPXXsetintparam(env, CPXPARAM_LPMethod, CPX_ALG_CONCURRENT);
	if(status != 0) { return(status); }

	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_MININFEAS);
	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_DEFAULT);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_MAXINFEAS);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_PSEUDO);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_PSEUDOREDUCED);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Order, CPX_ON);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_AUTO);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_DYNAMIC);
	if(status != 0) { return status; }
	/**********************/

	/** Lazy Constraints **/
	/* Assure linear mappings between the presolved and original models */
	status = CPXXsetintparam(env, CPXPARAM_Preprocessing_Linear, CPX_OFF);
	if(status != 0) { return status; }

	/* Turn on traditional search for use with control callbacks */
	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
	if(status != 0) { return status; }

	/* Let MIP callbacks work on the original model */
	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
	if(status != 0) { return status; }
	/**********************/

	/******* Logging ******/
	//CPXXsetlogfilename(env, "1e.log", "a");

	/**********************/

	return status;
}

int edge_var_number(const Problem* prob, int day, int from, int to) {
	assert(0 <= day && day < prob->schedules);
	assert(from != to);
	assert(0 <= from && from < prob->N);
	assert(0 < to && to <= prob->N);
	assert(from != 0 || to != prob->N);
	//This is a picewise function that maps the edge from 'from' to 'to', taking into consideration
	//that the graph doesn't have loop edges.
	//Depot vertex is divided in depot_src (0) and depot_dst (n), and depot_src don't have incoming edges nor depot_dst outgoing.
	int edge_number;
	if(from == 0) {
		edge_number = day*(prob->N-1)+(to-1);
	} else {
		edge_number = prob->schedules*(prob->N-1)                	//depot_src outgoing edges
									+day*((prob->N-1)*(prob->N-2)+(prob->N-1))	//day offset, depot_src doesn't have outgoing edges
									+(prob->N-1)*(from-1)                       //row_offset
									+(to-1);                                    //col_offset

		if(to > from) { edge_number--; } //take into consideration loop edges not labeled/numbered.
	}

	return edge_number;
}

int vertex_var_number(const Problem* prob, int day, int vertex) {
	assert(0 <= day && day < prob->schedules);
	assert(0 <= vertex && vertex <= prob->N);
	return prob->schedules*prob->N*(prob->N-1)	//edge_vars_offset
				 +day*(prob->N+1)											//row_offset
				 +vertex;															//col_offset
}

static int CPXPUBLIC subtour_constraint_generator(CPXCENVptr env,
																									void* cbdata,
																									int wherefrom,
																									void* cbhandle,
																									int* useraction_p) {
	int status = 0;
	*useraction_p = CPX_CALLBACK_DEFAULT;
	auto* lazyconinfo = static_cast<lazy_constrain_info*>(cbhandle);
	int cur_numcols = CPXXgetnumcols(env, lazyconinfo->lp);
	auto* x = new double[cur_numcols];

	status = CPXXgetcallbacknodex(env, cbdata, wherefrom, x, 0, cur_numcols-1);
	if(status != 0) {
		std::cerr << "Failed to get node solution." << std::endl;
		return (status);
	}

	//Search for all (sub)tours in the integer solution
	for(int day = 0; day < lazyconinfo->prob->schedules; ++day) {
		// TODO: estas dos estructuras de datos deberían ser una clase/struct con sus métodos que mantengan la info consistente
		std::vector<std::vector<int>> subtours;
		std::unordered_map<int, int> client_subtour;

		for(int i = 0; i < lazyconinfo->prob->N+1; ++i) {
			if(x[vertex_var_number(lazyconinfo->prob, day, i)] == 1.0 && client_subtour.count(i) == 0) {
				subtours.emplace_back();
				int cur_client = i;
				do {
					subtours.back().push_back(cur_client);
					client_subtour[cur_client] = subtours.size()-1;

					int j = 1;
					while(j == cur_client
								|| (cur_client < lazyconinfo->prob->N
										&& j < lazyconinfo->prob->N
										&& x[edge_var_number(lazyconinfo->prob, day, cur_client, j)] != 1.0)) { ++j; }
					cur_client = j;

					assert(0 < cur_client && cur_client < lazyconinfo->prob->N+1);
				} while(client_subtour.count(cur_client) == 0);
			}
		}

		if(subtours.size() > 1) {
			auto* cutind = new int[cur_numcols];
			auto* cutval = new double[cur_numcols];
			double rhs = 1.0;
			for(auto subtour : subtours) {
				int cutnz = 0;

				for(auto i : subtour) {
					cutval[cutnz] = 1.0;
					cutind[cutnz] = vertex_var_number(lazyconinfo->prob, day, i);
					cutnz++;

					for(auto j : subtour) {
						if(i != j
							 && (i != 0 || j != lazyconinfo->prob->N)
							 && i != lazyconinfo->prob->N
							 && j != 0) {
							cutval[cutnz] = -1.0;
							cutind[cutnz] = edge_var_number(lazyconinfo->prob, day, i, j);
							cutnz++;
						}
					}
				}
				status = CPXXcutcallbackadd(env, cbdata, wherefrom, cutnz, rhs, 'G', cutind, cutval, CPX_USECUT_FORCE);
				if(status != 0) {
					std::cerr << "Failed to add cut." << std::endl;
					return (status);
				}
			}
			//Tell CPLEX that cuts have been created
			*useraction_p = CPX_CALLBACK_SET;

			delete[] cutind;
			delete[] cutval;
		}
	}

	return status;
}

int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp) {
	//Problem is minimization
	CPXXchgobjsen(env, lp, CPX_MIN);
	int status;

	//Model's variables definition
	int edge_vars_quantity = prob->schedules*prob->N*(prob->N-1); //complete directed graph with split depot, for prob->schedules days
	int vertex_vars_quantity = prob->schedules*(prob->N+1); //+1 because of depot_dst
	int vars_quantity = edge_vars_quantity+vertex_vars_quantity;

	auto* obj_var_coefs = new double[vars_quantity];
	auto* lower_bounds = new double[vars_quantity];
	auto* upper_bounds = new double[vars_quantity];
	auto* column_types = new char[vars_quantity];

	int def_vars = 0;
	for(int day = 0; day < prob->schedules; ++day) {
		//depot_src vars
		obj_var_coefs[vertex_var_number(prob, day, 0)] = 0.0;
		lower_bounds[vertex_var_number(prob, day, 0)] = 0.0;
		upper_bounds[vertex_var_number(prob, day, 0)] = 1.0;
		column_types[vertex_var_number(prob, day, 0)] = 'B';
		def_vars++;

		//depot_dst vars
		obj_var_coefs[vertex_var_number(prob, day, prob->N)] = 0.0;
		lower_bounds[vertex_var_number(prob, day, prob->N)] = 0.0;
		upper_bounds[vertex_var_number(prob, day, prob->N)] = 1.0;
		column_types[vertex_var_number(prob, day, prob->N)] = 'B';
		def_vars++;

		for(int i = 1; i < prob->N; ++i) {
			//Client edges
			for(int j = 1; j < prob->N; ++j) {
				if(i != j) {
					obj_var_coefs[edge_var_number(prob, day, i, j)] = prob->costs[i][j];
					lower_bounds[edge_var_number(prob, day, i, j)] = 0.0;
					upper_bounds[edge_var_number(prob, day, i, j)] = 1.0;
					column_types[edge_var_number(prob, day, i, j)] = 'B';
					def_vars++;
				}
			}

			//Client vars
			obj_var_coefs[vertex_var_number(prob, day, i)] = 0.0;
			lower_bounds[vertex_var_number(prob, day, i)] = 0.0;
			upper_bounds[vertex_var_number(prob, day, i)] = 1.0;
			column_types[vertex_var_number(prob, day, i)] = 'B';
			def_vars++;

			//depot_src edges
			obj_var_coefs[edge_var_number(prob, day, 0, i)] = prob->costs[0][i];
			lower_bounds[edge_var_number(prob, day, 0, i)] = 0.0;
			upper_bounds[edge_var_number(prob, day, 0, i)] = 1.0;
			column_types[edge_var_number(prob, day, 0, i)] = 'B';
			def_vars++;

			//depot_dst edges
			obj_var_coefs[edge_var_number(prob, day, i, prob->N)] = prob->costs[i][0];
			lower_bounds[edge_var_number(prob, day, i, prob->N)] = 0.0;
			upper_bounds[edge_var_number(prob, day, i, prob->N)] = 1.0;
			column_types[edge_var_number(prob, day, i, prob->N)] = 'B';
			def_vars++;
		}
	}
	//Add variables defs/cols to lp
	assert(def_vars == vars_quantity);
	status = CPXXnewcols(env, lp, vars_quantity, obj_var_coefs, lower_bounds, upper_bounds, column_types, nullptr);
	if(status != 0) {
		std::cerr << "Variable definitions/Columns could not be defined." << std::endl;
		return status;
	}

	//Model's constraints
	auto* rhs = new double[1];
	auto* sense = new char[1];
	auto* rmatbeg = new CPXNNZ[1];
	auto* rmatind = new int[vars_quantity];
	auto* rmatind2 = new int[vars_quantity];
	auto* rmatval = new double[vars_quantity];

	//Vehicle capacity constraint
	// Sum_ij e_ijk dj <= Capacity (k=1..prob->schedules)
	rmatbeg[0] = 0; //Just one constraint
	rhs[0] = prob->K;
	sense[0] = 'L';
	for(int day = 0; day < prob->schedules; ++day) {
		int def_vars = 0;
		for(int i = 0; i < prob->N; ++i) {
			for(int j = 1; j < prob->N; ++j) {
				if(i != j && (i != 0 || j != prob->N)) {
					rmatval[def_vars] = prob->demands[j]; //Variable coef
					rmatind[def_vars] = edge_var_number(prob, day, i, j); //Variable number
					def_vars++;
				}
			}
		}
		//Add constraint
		assert(def_vars == prob->N*(prob->N-1)-(prob->N-1));
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}
	}

	//Clients are visited accordingly to their setup (all days, cust_clients_freq times each 'schedule days')
	// Sum_k(u_ik) = prob->schedules  (for all i in every_day_clients)
	// Sum_k(u_ik) = prob->cust_clients_freq  (for all i in cust_clients)
	rmatbeg[0] = 0; //Just one constraint
	sense[0] = 'E';
	for(int i = 0; i < prob->N+1; ++i) {
		int def_vars = 0;
		for(int day = 0; day < prob->schedules; ++day) {
			rmatval[day] = 1.0; //Variable coef
			rmatind[day] = vertex_var_number(prob, day, i); //Variable number
			def_vars++;
		}

		if(prob->depots+prob->every_day <= i && i < prob->N) {
			rhs[0] = static_cast<double>(prob->cust_clients_freq);
		} else {
			rhs[0] = static_cast<double>(prob->schedules);
		}

		//Add constraint
		assert(def_vars == prob->schedules);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) { return status; }
	}

	//Each day a client is visited, vehicle arrives and leaves exactly one time.
	// Sum_i(e_ijk) = u_jk  (k=1..prob->schedules, j=1..n, i=1..n-1, i#j)
	rmatbeg[0] = 0; //Just one constraint
	rhs[0] = 0.0;
	sense[0] = 'E';
	for(int day = 0; day < prob->schedules; ++day) {
		for(int j = 1; j < prob->N+1; ++j) {
			rmatval[0] = -1.0;
			rmatind[0] = vertex_var_number(prob, day, j);
			rmatind2[0] = vertex_var_number(prob, day, j);

			int def_vars = 1;
			for(int i = 0; i < prob->N; ++i) {
				if(i != j && (i != 0 || j != prob->N)) {
					rmatval[def_vars] = 1.0; //Variable coef
					rmatind[def_vars] = edge_var_number(prob, day, i, j); //Variable number
					def_vars++;
				}
			}
			//Add constraints
			assert(def_vars == prob->N-1+1);
			status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added." << std::endl;
				return (status);
			}
		}
	}

	// Sum_j(e_ijk) = u_ik  (k=1..prob->schedules, i=0..n-1, j=1..n, i#j)
	rmatbeg[0] = 0; //Just one constraint
	rhs[0] = 0.0;
	sense[0] = 'E';
	for(int day = 0; day < prob->schedules; ++day) {
		for(int i = 0; i < prob->N; ++i) {
			rmatval[0] = -1.0;
			rmatind[0] = vertex_var_number(prob, day, i);
			rmatind2[0] = vertex_var_number(prob, day, i);

			int def_vars = 1;
			for(int j = 1; j < prob->N+1; ++j) {
				if(i != j && (i != 0 || j != prob->N)) {
					rmatval[def_vars] = 1.0; //Variable coef
					rmatind[def_vars] = edge_var_number(prob, day, i, j); //Variable number
					def_vars++;
				}
			}
			//Add constraints
			assert(def_vars == prob->N-1+1);
			status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added." << std::endl;
				return (status);
			}
		}
	}

	//Subtours elimination
	/* Set up to use MIP lazyconstraint callback.*/
	auto* lazyconinfo = new lazy_constrain_info(lp, prob);
	status = CPXXsetlazyconstraintcallbackfunc(env, subtour_constraint_generator, lazyconinfo);
	if(status != 0) {
		std::cerr << "Error adding lazy constraing callback function." << std::endl;
		return (status);
	}

	//Free pointers
	delete[] obj_var_coefs;
	delete[] lower_bounds;
	delete[] upper_bounds;
	delete[] column_types;
	delete[] rhs;
	delete[] sense;
	delete[] rmatbeg;
	delete[] rmatind;
	delete[] rmatind2;
	delete[] rmatval;

	return 0;
}

int initial_solution_mip(Problem* prob, CPXENVptr env, CPXLPptr lp) {
	int status = 0;
	int edge_vars_quantity = prob->schedules*prob->N*(prob->N-1); //complete directed graph with split depot, for prob->schedules days
	int vertex_vars_quantity = prob->schedules*(prob->N+1); //+1 because of depot_dst
	int vars_quantity = edge_vars_quantity+vertex_vars_quantity;
	auto* rmatind = new int[vars_quantity];
	auto* rmatval = new double[vars_quantity];

	status = CPXXsetintparam(env, CPX_PARAM_ADVIND, 1);
	if(status != 0) {
		std::cerr << "Error setting initial feasible solution parameter." << std::endl;
		return (status);
	}

	int def_vars = 0;
	int next_cust_client = prob->depots+prob->every_day;
	for(int day = 0; day < prob->schedules; ++day) {
		rmatval[def_vars] = 1.0;
		rmatind[def_vars] = vertex_var_number(prob, day, 0); //depot_src //TODO: no esta claro que hacer si hay más de 1 depot
		def_vars++;

		double availably_capacity = prob->K;
		int cur_vertex = 0; //depot_src
		for(int i = prob->depots; i < prob->depots+prob->every_day; ++i) {
			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = vertex_var_number(prob, day, i);
			def_vars++;

			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = edge_var_number(prob, day, cur_vertex, i);
			def_vars++;

			cur_vertex = i;
			availably_capacity -= prob->demands[i];
		}

		for(;next_cust_client < prob->N && availably_capacity - prob->demands[next_cust_client] >= 0; ++next_cust_client) {
			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = vertex_var_number(prob, day, next_cust_client);
			def_vars++;

			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = edge_var_number(prob, day, cur_vertex, next_cust_client);
			def_vars++;

			availably_capacity -= prob->demands[next_cust_client];
			cur_vertex = next_cust_client;
		}

		//depot_dst
		rmatval[def_vars] = 1.0;
		rmatind[def_vars] = vertex_var_number(prob, day, prob->N);
		def_vars++;

		rmatval[def_vars] = 1.0;
		rmatind[def_vars] = edge_var_number(prob, day, cur_vertex, prob->N);
		def_vars++;
	}

	assert(def_vars
				 == prob->schedules*(prob->depots+prob->every_day+1) //every day + depot_dst vertexes
						+prob->schedules*(prob->depots+prob->every_day)  //every day + depot_dst edges
						+prob->cust_clients //cust clients + depot_dst vertexes
						+prob->cust_clients-1 // cust clients + depot_dst edges
						+prob->schedules-1); //every day + cust clients edges union
	CPXNNZ beg[1] = {0};
	int effortlevel[1] = {CPX_MIPSTART_AUTO};
	status = CPXXaddmipstarts(env, lp, 1, def_vars, beg, rmatind, rmatval, effortlevel, nullptr);
	if(status != 0) {
		std::cerr << "Error adding initial feasible solution." << std::endl;
		return (status);
	}

	delete[] rmatval;
	delete[] rmatind;

	return status;
}

int solveMIP(const Problem* prob,
						 CPXENVptr env,
						 CPXLPptr lp,
						 double* const objval,
						 std::vector<std::vector<int>>& schedules) {
	int status = 0;
	status = CPXXmipopt(env, lp);
	if(status != 0) { return status; }

	//Retrieve optimal solution
	status = CPXXgetobjval(env, lp, objval);
	if(status != 0) { return status; }

	int cur_numcols = CPXXgetnumcols(env, lp);
	auto* vars = new double[cur_numcols];
	status = CPXXgetx(env, lp, vars, 0, cur_numcols-1); // Get solution's variable values
	if(status != 0) { return status; }

	//Build solution according to instance vertex numbers
	schedules.clear();
	for(int day = 0; day < prob->schedules; ++day) {
		schedules.emplace_back();
		int cur_farm = 0;
		do {
			schedules[day].push_back(cur_farm);
			int j = 1;
			while(j == cur_farm || vars[edge_var_number(prob, day, cur_farm, j)] != 1.0) { ++j; }
			/*while(j == cur_farm
						|| (cur_farm < prob->N
								&& j < prob->N
								&& vars[edge_var_number(prob, day, cur_farm, j)] != 1.0)) { ++j; }*/
			cur_farm = j;
		} while(cur_farm != prob->N);
		schedules[day].push_back(0);
	}

	delete[] vars;
	return status;
}

/**********************/
int solve(Problem* prob, std::vector<std::vector<int>>& schedules, double& objval) {
	// CPlex environment
	int status = 0;
	CPXENVptr env = nullptr;
	CPXLPptr lp = nullptr;

	status = initialize_structures(env, lp);
	if(status != 0) {
		free_structures(env, lp);
		return status;
	}

	status = set_parameters(env);
	if(status != 0) {
		std::cout << "Error setting CPLEX parameters" << std::endl;
		free_structures(env, lp);
		return status;
	}

	std::cout << "Initializing MIP" << std::endl;
	status = initialize_mip(prob, env, lp);
	if(status != 0) {
		free_structures(env, lp);
		return status;
	}

	/*std::cout << "Providing initial integer solution" << std::endl;
	status = initial_solution_mip(prob, env, lp);
	if(status != 0) {
		free_structures(env, lp);
		return status;
	}*/

	//DEBUG.
	/*status = CPXXwriteprob(env, lp, "1e.lp", nullptr);
	if(status != 0)
	{
		free_structures(env, lp);
		return status;
	}*/

	std::cout << "Solving" << std::endl;
	status = solveMIP(prob, env, lp, &objval, schedules);
	if(status != 0) {
		free_structures(env, lp);
		return status;
	}
	free_structures(env, lp);

	return status;
}