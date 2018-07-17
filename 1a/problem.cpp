#include<cassert>
#include<iostream>
#include<limits>
#include<unordered_map>
#include<thread>
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
	assert(0 <= oil && oil < oil_q);
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

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_RINSHeur, 0);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_FPHeur, 0);
	if(status != 0) { return(status); }

	status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_LBHeur, CPX_ON);
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

	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_DFS);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_BESTBOUND);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_BESTEST);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_NodeSelect, CPX_NODESEL_BESTEST_ALT);
	//if(status != 0) { return status; }

	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Branch, -1);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Branch, 0);
	//status = CPXXsetintparam(env, CPXPARAM_MIP_Strategy_Branch, 1);
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

	/******* Logging ******/
	//CPXXsetlogfilename(env, "1a.log", "a");

	/**********************/

	return status;
}

int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp) {
	//Problem is maximization
	CPXXchgobjsen(env, lp, CPX_MAX);
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
			//lower_bounds[refined_var_num] = 0.0;
			//upper_bounds[refined_var_num] = std::numeric_limits<double>::max();
			column_types[refined_var_num] = 'I';
			def_vars++;

			int storage_var_num = vars.var_num(mip_vars::var_type::storage, oil, month);
			obj_var_coefs[storage_var_num] = -prob->oil_storage_cost_per_month[oil][month];
			//lower_bounds[storage_var_num] = 0.0;
			//upper_bounds[storage_var_num] = std::numeric_limits<double>::max();
			column_types[storage_var_num] = 'I';
			def_vars++;

			int buy_var_num = vars.var_num(mip_vars::var_type::buy, oil, month);
			obj_var_coefs[buy_var_num] = -prob->oil_buy_cost_per_month[oil][month];
			//lower_bounds[buy_var_num] = 0.0;
			//upper_bounds[buy_var_num] = std::numeric_limits<double>::max();
			column_types[buy_var_num] = 'I';
			def_vars++;

			int oil_switch_var_num = vars.var_num(mip_vars::var_type::oil_switch, oil, month);
			obj_var_coefs[oil_switch_var_num] = 0.0;
			//lower_bounds[oil_switch_var_num] = 0.0;
			//upper_bounds[oil_switch_var_num] = 1.0;
			column_types[oil_switch_var_num] = 'B';
			def_vars++;
		}
	}

	//Add variables defs/cols to lp
	assert(def_vars == vars.total());
	status = CPXXnewcols(env, lp, def_vars, obj_var_coefs, nullptr, nullptr, column_types, nullptr);
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
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
			return (status);
		}

		//Hardness
		rhs[0] = 0.0;
		sense[0] = 'G';
		def_vars = 0;
		for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
			rmatval[def_vars] = prob->max_hardness-prob->hardness[oil];
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;
		}
		assert(def_vars == prob->V+prob->NV);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
			return (status);
		}

		def_vars = 0;
		for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
			rmatval[def_vars] = prob->hardness[oil]-prob->min_hardness;
			rmatind[def_vars] = vars.var_num(mip_vars::var_type::refined, oil, month);
			def_vars++;
		}
		assert(def_vars == prob->V+prob->NV);
		status = CPXXaddrows(env, lp, 0, 1, def_vars, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
				std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
				return (status);
			}

			//Oil switch
			rhs[0] = 0.0;
			sense[0] = 'G';
			//rmatval[0] = std::numeric_limits<double>::max();
			rmatval[0] = std::max(prob->max_refined_per_month_NV, prob->max_refined_per_month_V);
			rmatind[0] = vars.var_num(mip_vars::var_type::oil_switch, oil, month);

			rmatval[1] = -1.0;
			rmatind[1] = vars.var_num(mip_vars::var_type::refined, oil, month);

			status = CPXXaddrows(env, lp, 0, 1, 2, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status != 0) {
				std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
				std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
			return (status);
		}

		//Final storage
		sense[0] = 'G';
		rmatval[0] = 1.0;
		rmatind[0] = vars.var_num(mip_vars::var_type::storage, oil, prob->M-1);
		status = CPXXaddrows(env, lp, 0, 1, 1, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status != 0) {
			std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
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
				std::cerr << "Row/Constraint could not be added. Line: " << __LINE__ << std::endl;
				return (status);
			}
		}
	}

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
						 std::vector<std::vector<month_oil_solution>>& solution) {
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

	//Build solution according to problem
	mip_vars mip_var_nums(prob);
	solution.clear();
	for(int month = 0; month < prob->M; ++month) {
		solution.emplace_back();
		for(int oil = 0; oil < prob->V+prob->NV; ++ oil) {
			month_oil_solution res;
			res.month = month;
			res.tons_buyed = vars[mip_var_nums.var_num(mip_vars::var_type::buy, oil, month)];
			res.tons_refined = vars[mip_var_nums.var_num(mip_vars::var_type::refined, oil, month)];
			res.tons_storaged = vars[mip_var_nums.var_num(mip_vars::var_type::storage, oil, month)];
			solution.back().push_back(res);
		}
	}

	delete[] vars;
	return status;
}

/**********************/
int solve(Problem* prob, double& objval, std::vector<std::vector<month_oil_solution>>& solution) {
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
	status = solveMIP(prob, env, lp, &objval, solution);
	if(status != 0) {
		free_structures(env, lp);
		return status;
	}
	free_structures(env, lp);

	return status;
}