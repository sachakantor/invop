#include<cassert>
#include<iostream>
#include<limits>
#include<unordered_map>
#include<vector>

#include<problem.hpp>
#include<ilcplex/cplexx.h>


/**********************/
Problem::Problem() = default;

/**********************/
mip_vars::mip_vars(const Problem* prob)
	: oil_q(prob->V+prob->NV),
		month_q(prob->M),
		refined_base_offset(0),
		storage_base_offset(refined_base_offset + oil_q*month_q),
		buy_base_offset(storage_base_offset + oil_q*month_q),
		oil_use_switch_base_offset(buy_base_offset + oil_q*month_q) {}

int mip_vars::var_num(mip_vars::var_type t, int oil, int month) const {
	assert(0 <= month && month < month_q); //TODO: se debería hacer verdadero manejo de errores acá
	assert(0 <= oil && oil <= oil_q);
	int base_offset = 0;
	switch(t) {
		case mip_vars::var_type::refined 		: base_offset = refined_base_offset; break;
		case mip_vars::var_type::storage 		: base_offset = storage_base_offset; break;
		case mip_vars::var_type::buy 				: base_offset = buy_base_offset; break;
		case mip_vars::var_type::oil_switch : base_offset = oil_use_switch_base_offset; break;
	}

	return base_offset+ oil*month_q + month;
}

int mip_vars::total() const {
	return oil_use_switch_base_offset + this->oil_use_switch_total();
}

int mip_vars::refined_total() const {
	return oil_q*month_q;
}

int mip_vars::storage_total() const {
	return oil_q*month_q;
}

int mip_vars::buy_total() const {
	return oil_q*month_q;
}

int mip_vars::oil_use_switch_total() const {
	return oil_q*month_q;
}

/**********************/
//lazy_constrain_info::lazy_constrain_info(CPXLPptr lp, Problem* prob) : lp(lp), prob(prob) {}

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

	status = CPXXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
	if(status != 0) { return status; }

	/****** Output ******/
	status = CPXXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status != 0) { return status; }

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPDISPLAY, 5);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1 );
	if(status != 0) { return(status); }*/
	/**********************/

	/****** Presolve ******/
	status = CPXXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if(status != 0) { return status; }

	/*// Ver bien por que sin esto no actualiza las duales.
	// Este es el workaround que encontro Agus (Crack).
	status = CPXXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 0);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_PROBE, -1);
	if(status != 0) { return(status); }*/
	/**********************/

	/** CPLEX heuristics **/
	/*status = CPXXsetintparam(env, CPX_PARAM_HEURFREQ, -1);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_RINSHEUR, -1);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FPHEUR, -1);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_LBHEUR, -1);
	if(status != 0) { return(status); }*/
	/**********************/

	/***** CPLEX cuts *****/
	/*status = CPXXsetintparam(env, CPX_PARAM_CUTPASS, 100);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_CLIQUES, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_COVERS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_DISJCUTS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FLOWCOVERS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FLOWPATHS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FRACCUTS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_GUBCOVERS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MCFCUTS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_IMPLBD, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MIRCUTS, -1 );
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_ZEROHALFCUTS, -1 );
	if(status != 0) { return(status); }*/
	/**********************/

	/** Execution decisions **/
	status = CPXXsetdblparam(env, CPX_PARAM_TILIM, 3600.0);
	if(status != 0) { return status; }

	/*status = CPXXsetintparam(env, CPX_PARAM_NODELIM, 1);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_PARALLELMODE, 1);
	if(status != 0) { return(status); }*/

	status = CPXXsetintparam(env, CPX_PARAM_THREADS, 8);
	if(status != 0) { return status; }

	//status = CPXXsetintparam(env, CPX_PARAM_NODESEL, CPX_NODESEL_DFS);
	status = CPXXsetintparam(env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTBOUND);
	if(status != 0) { return status; }

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_PRIMAL); CHECKSTATUS;
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_BARCROSSALG, -1); CHECKSTATUS;
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL ); CHECKSTATUS;
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_VARSEL, CPX_VARSEL_MININFEAS); CHECKSTATUS;
	if(status != 0) { return(status); }*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPORDIND, CPX_ON); CHECKSTATUS;
	if(status != 0) { return(status); }*/

	status = CPXXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	if(status != 0) { return status; }
	/**********************/

	/** Lazy Constraints **/
	/* Assure linear mappings between the presolved and original models */
	/*status = CPXXsetintparam(env, CPXPARAM_Preprocessing_Linear, CPX_OFF);
	if(status != 0) { return status; }*/

	/* Turn on traditional search for use with control callbacks */
	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL); //TODO: mepa que se solapa con otra
	if(status != 0) { return status; }*/

	/* Let MIP callbacks work on the original model */
	/*status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
	if(status != 0) { return status; }*/
	/**********************/

	/******* Logging ******/
	/*CPXFILEptr logfile = nullptr;
	logfile = CPXXfopen("BPGC.log", "a");
	CPXXsetlogfile(env, logfile);*/
	/**********************/

	return status;
}

int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp) {
	//Problem is minimization
	CPXXchgobjsen(env, lp, CPX_MIN);
	int status;

	//Model's variables definition
	mip_vars vars(prob);

	auto* obj_var_coefs = new double[vars.total()];
	auto* lower_bounds = new double[vars.total()];
	auto* upper_bounds = new double[vars.total()];
	auto* column_types = new char[vars.total()];

	int def_vars = 0;
	for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
		for(int month = 0; month < prob->M; ++month) {
			int refined_var_num = vars.var_num(mip_vars::var_type::refined, oil, month);
			obj_var_coefs[refined_var_num] = prob->oil_profit_per_month[oil][month];
			lower_bounds[refined_var_num] = 0.0;
			upper_bounds[refined_var_num] = std::numeric_limits<double>::max();
			column_types[refined_var_num] = 'I';
			def_vars++;

			int storage_var_num = vars.var_num(mip_vars::var_type::storage, oil, month);
			obj_var_coefs[storage_var_num] = -prob->oil_storage_cost_per_month[oil][month];
			lower_bounds[storage_var_num] = 0.0;
			upper_bounds[storage_var_num] = std::numeric_limits<double>::max();
			column_types[storage_var_num] = 'I';
			def_vars++;

			int buy_var_num = vars.var_num(mip_vars::var_type::buy, oil, month);
			obj_var_coefs[buy_var_num] = -prob->oil_buy_cost_per_month[oil][month];
			lower_bounds[buy_var_num] = 0.0;
			upper_bounds[buy_var_num] = std::numeric_limits<double>::max();
			column_types[buy_var_num] = 'I';
			def_vars++;

			int oil_switch_var_num = vars.var_num(mip_vars::var_type::oil_switch, oil, month);
			obj_var_coefs[oil_switch_var_num] = 0.0;
			lower_bounds[oil_switch_var_num] = 0.0;
			upper_bounds[oil_switch_var_num] = 1.0;
			column_types[oil_switch_var_num] = 'B';
			def_vars++;
		}
	}

	//Add variables defs/cols to lp
	assert(def_vars == vars.total());
	status = CPXXnewcols(env, lp, def_vars, obj_var_coefs, lower_bounds, upper_bounds, column_types, nullptr);
	if(status != 0) {
		std::cerr << "Variable definitions/Columns could not be defined." << std::endl;
		return status;
	}

	//Model's constraints
	auto* rhs = new double[1];
	auto* sense = new char[1];
	auto* rmatbeg = new CPXNNZ[1];
	auto* rmatind = new int[vars.total()];
	auto* rmatval = new double[vars.total()];

	rmatbeg[0] = 0; //Just one constraint per addition

	/* Refinery capacity
	 * Sum_{i,m}(r_im) <= max_refined_per_month_vegetal  forall i in veg_oils, m in months
	 * Sum_{i,m}(r_im) <= max_refined_per_month_non_veg  forall i in non_veg_oils, m in months
	 *
	 * Hardness
	 * min_hardness <= Sum_i(h_i r_im) / Sum_i(r_im) <= max_hardness  forall i in oils, m in months
	 *
	 * Max oil types used per month
	 * Sum_i(y_im) <= max_oils_used_per_month  forall i in oils, m in months
	 *
	 * Vegetable force first non-veg use
	 * Sum_i(y_im) <= I y_jm  forall i in veg_oils, m in months, k=first non-veg oil, I=|veg_oils|
	 */
	sense[0] = 'L';
	for(int month = 0; month < prob->M; ++month) {
		//Vegetable oils capacity
		rhs[0] = prob->max_refined_per_month_V;
		int def_vars = 0;
		for(int oil = 0; oil < prob->V; ++oil) {
			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;
		}
		assert(def_vars == prob->V);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		//Non-Vegetable oils capacity
		rhs[0] = prob->max_refined_per_month_NV;
		def_vars = 0;
		for(int oil = prob->V; oil < prob->V+prob->NV; ++oil) {
			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;
		}
		assert(def_vars == prob->NV);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		//Hardness
		rhs[0] = 0.0;
		sense[0] = 'G';
		def_vars = 0;
		for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
			rmatval[def_vars] = prob->max_hardness;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;

			rmatval[def_vars] = -prob->hardness[oil];
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;
		}
		assert(def_vars == 2*(prob->V+prob->NV));
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		def_vars = 0;
		for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
			rmatval[def_vars] = -prob->min_hardness;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;

			rmatval[def_vars] = prob->hardness[oil];
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;
		}
		assert(def_vars == 2*(prob->V+prob->NV));
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		//Max oil types used per month
		rhs[0] = prob->max_oils_per_month;
		sense[0] = 'L';
		def_vars = 0;
		for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::oil_switch, oil, month);
			def_vars++;
		}
		assert(def_vars == prob->V+prob->NV);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		//Vegetable force first non-veg use
		rhs[0] = 0.0;
		sense[0] = 'L';
		def_vars = 0;
		for(int veg_oil = 0; veg_oil < prob->V; ++veg_oil) {
			rmatval[def_vars] = 1.0;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::oil_switch, veg_oil, month);
			def_vars++;
		}
		rmatval[def_vars] = -static_cast<double>(prob->V);
		rmatind[def_vars] = vars.var_num(mip_vars::var_type::oil_switch, prob->V, month); //prob->V = first non-veg idx
		def_vars++;

		assert(def_vars == prob->V+1);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

	}

	/* Storage capacity
	 * s_im <= storage_capacity[i]  forall i in oils, m in months
	 *
	 * Oil_switchs/refined_oils variables' consistency
	 * r_im <= y_im M  forall i in oils, m in months, M=infinity
	 *
	 * Minimum tons per oil if refined
	 * r_im >= y_im 20  forall i in oils, m in months
	 */
	for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
		for(int month = 0; month < prob->M; ++month) {
			//Storage capacity
			rhs[0] = prob->storage_capacity[oil];
			sense[0] = 'L';
			rmatval[0] = 1.0;
			rmatind[0] = vars.var_num(mip_vars::var_type::storage, oil, month);

			status = CPXXaddrows(env, lp, 0, 1, 1, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added." << std::endl;
				return (status);
			}

			//Oil switch
			rhs[0] = 0.0;
			sense[0] = 'G';
			rmatval[0] = std::numeric_limits<double>::max();
			rmatind[0] = vars.var_num(mip_vars::var_type::oil_switch, oil, month);

			rmatval[1] = -1.0;
			rmatind[1] = vars.var_num(mip_vars::var_type::refined, oil, month);

			status = CPXXaddrows(env, lp, 0, 1, 2, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added." << std::endl;
				return (status);
			}

			//Min tons per oil if refined
			rhs[0] = 0.0;
			sense[0] = 'G';
			rmatval[0] = 1.0;
			rmatind[0] = vars.var_num(mip_vars::var_type::refined, oil, month);

			rmatval[1] = -prob->min_tons_oil_if_refined[oil];
			rmatind[1] = vars.var_num(mip_vars::var_type::oil_switch, oil, month);

			status = CPXXaddrows(env, lp, 0, 1, 2, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added." << std::endl;
				return (status);
			}
		}
	}

	/* Initial storaged oil
	 * s_i0 = initial_storaged_i  forall i in oils
	 *
	 * Final storaged oil
	 * s_iM >= initial_storaged_i  forall i in oils, M=last month index
	 *
	 * Storage flow consistency through months
	 * s_im = b_im + s_i(m-1) - r_im  forall i in oils, m=1..M
	 */
	for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
		rhs[0] = prob->initial_storage[oil];

		//Initial storage
		sense[0] = 'E';
		rmatval[0] = 1.0;
		rmatind[0] = vars.var_num(mip_vars::var_type::storage, oil, 0);
		status = CPXXaddrows(env, lp, 0, 1, 1, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		//Final storage
		sense[0] = 'G';
		rmatval[0] = 1.0;
		rmatind[0] = vars.var_num(mip_vars::var_type::storage, oil, prob->M-1);
		status = CPXXaddrows(env, lp, 0, 1, 1, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added." << std::endl;
			return (status);
		}

		//Storage Flow
		rhs[0] = 0.0;
		sense[0] = 'E';
		for(int month = 1; month < prob->M; ++month) {
			rmatval[0] = 1.0;
			rmatind[0] = vars.var_num(mip_vars::var_type::storage, oil, month);

			rmatval[1] = 1.0;
			rmatind[1] = vars.var_num(mip_vars::var_type::refined, oil, month);

			rmatval[2] = -1.0;
			rmatind[2] = vars.var_num(mip_vars::var_type::storage, oil, month-1);

			rmatval[3] = -1.0;
			rmatind[3] = vars.var_num(mip_vars::var_type::buy, oil, month);

			status = CPXXaddrows(env, lp, 0, 1, 4, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added." << std::endl;
				return (status);
			}
		}
	}

	// Este sector de codigo permite proveerle a CPLEX una solucion factible de entrada
	/*status = CPXsetintparam(env, CPX_PARAM_ADVIND, 1);
	if(status != 0)exit(-1);

	int beg[1] = {0};
	int nzcnt = prob->N+prob->bestObj;
	int* varindices = (int*) malloc(sizeof(int)*nzcnt);
	for(int i = 0; i < prob->N; i++)
	{
		varindices[i] = i*prob->bestObj+prob->colores[i];
	}
	for(int i = 0; i < prob->bestObj; i++)
	{
		varindices[i+prob->N] = i+prob->N*prob->bestObj;
	}


	double* values = (double*) malloc(sizeof(double)*nzcnt);
	for(int i = 0; i < nzcnt; ++i)values[i] = 1.0;
	int effortlevel[1] = {CPX_MIPSTART_AUTO};
	status = CPXaddmipstarts(env, lp, 1, nzcnt, beg, varindices, values, effortlevel, NULL);

	free(varindices);
	free(values);*/

	delete[] obj_var_coefs;
	delete[] lower_bounds;
	delete[] upper_bounds;
	delete[] column_types;
	delete[] rhs;
	delete[] sense;
	delete[] rmatbeg;
	delete[] rmatind;
	delete[] rmatval;

	return 0;
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
	/*schedules.clear();
	for(int day = 0; day < prob->schedules; ++day) {
		schedules.emplace_back();
		int cur_farm = 0;
		do {
			schedules[day].push_back(cur_farm);
			int j = 1;
			while(j == cur_farm || vars[edge_var_number(prob, day, cur_farm, j)] != 1.0) ++j;
			cur_farm = j;
		} while(cur_farm != prob->N);
		schedules[day].push_back(0);
	}*/

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

	// DEBUG.
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