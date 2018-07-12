#include<cassert>
#include<iostream>
#include<vector>

#include<problem.hpp>
#include<ilcplex/cplexx.h>


/**********************/
Problem::Problem() = default;
Problem::Problem(int n) : costs(n, std::vector<double>(n,0.0)), demands(n,0.0) {}

/**********************/
//lazy_constrain_info::lazy_constrain_info(CPXLPptr lp, Problem* prob) : lp(lp), prob(prob) {}

/**********************/
int initialize_structures(CPXENVptr& env, CPXLPptr& lp)
{

	int status = 0;
	env = CPXXopenCPLEX(&status);

	//std::cout<<env<<" "<<status<<std::endl;
	if(env == nullptr)
	{
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXXgeterrorstring(env, status, errmsg);
		std::cerr<<"Could not open CPLEX env\n"<<errmsg<<std::endl;
		return 1;
	}

	/* Create the problem. */
	std::string lpname = "1e.lp";
	lp = CPXXcreateprob(env, &status, lpname.c_str());

	if(lp == nullptr)
	{
		std::cerr<<"Failed to create LP"<<std::endl;
		return 1;
	}

	return 0;
}

void free_structures(CPXENVptr env, CPXLPptr lp)
{
	int status = 0;

	if(lp != nullptr)
	{
		status = CPXXfreeprob(env, &lp);
		if(status)
		{
			std::cerr<<"Problem feeing lp - status: "<<status<<std::endl;
		}
	}

	if(env != nullptr)
	{
		status = CPXXcloseCPLEX(&env);
		if(status)
		{
			std::cerr<<"Problem closing CPLEX environment"<<std::endl;
			char errmsg[CPXMESSAGEBUFSIZE];
			CPXXgeterrorstring(env, status, errmsg);
			std::cerr<<errmsg<<std::endl;
		}
	}

}

int set_parameters(CPXENVptr env)
{
	int status;

	status = CPXXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
	if(status) return (-1);

	/****** Output ******/
	status = CPXXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status) return (-1);

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPDISPLAY, 5);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPINTERVAL, 1 );
	if(status) return(-1);*/
	/**********************/

	/****** Presolve ******/
	status = CPXXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if(status) return (-1);

	/*// Ver bien por que sin esto no actualiza las duales.
	// Este es el workaround que encontro Agus (Crack).
	status = CPXXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 0);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_PROBE, -1);
	if(status) return(-1);*/
	/**********************/

	/** CPLEX heuristics **/
	/*status = CPXXsetintparam(env, CPX_PARAM_HEURFREQ, -1);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_RINSHEUR, -1);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FPHEUR, -1);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_LBHEUR, -1);
	if(status) return(-1);*/
	/**********************/

	/***** CPLEX cuts *****/
	/*status = CPXXsetintparam(env, CPX_PARAM_CUTPASS, 100);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_CLIQUES, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_COVERS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_DISJCUTS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FLOWCOVERS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FLOWPATHS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_FRACCUTS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_GUBCOVERS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MCFCUTS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_IMPLBD, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MIRCUTS, -1 );
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_ZEROHALFCUTS, -1 );
	if(status) return(-1);*/
	/**********************/

	/** Execution decisions **/
	status = CPXXsetdblparam(env, CPX_PARAM_TILIM, 3600.0);
	if(status) return (-1);

	/*status = CPXXsetintparam(env, CPX_PARAM_NODELIM, 1);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_PARALLELMODE, 1);
	if(status) return(-1);*/

	status = CPXXsetintparam(env, CPX_PARAM_THREADS, 8);
	if(status) return (-1);

	//status = CPXXsetintparam(env, CPX_PARAM_NODESEL, CPX_NODESEL_DFS);
	status = CPXXsetintparam(env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTBOUND);
	if(status) return (-1);

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_PRIMAL); CHECKSTATUS;
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_BARCROSSALG, -1); CHECKSTATUS;
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL ); CHECKSTATUS;
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_VARSEL, CPX_VARSEL_MININFEAS); CHECKSTATUS;
	if(status) return(-1);*/

	/*status = CPXXsetintparam(env, CPX_PARAM_MIPORDIND, CPX_ON); CHECKSTATUS;
	if(status) return(-1);*/

	status = CPXXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	if(status) return (-1);
	/**********************/

	/** Lazy Constraints **/
	/* Assure linear mappings between the presolved and original models */
	/*status = CPXXsetintparam (env, CPXPARAM_Preprocessing_Linear, 0);
	if(status) return (-1);*/

	/* Turn on traditional search for use with control callbacks */
	/*status = CPXXsetintparam (env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL); //mepa que se solapa con otra
	if(status) return (-1);*/

	/* Let MIP callbacks work on the original model */
	/*status = CPXXsetintparam (env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
	if(status) return (-1);*/
	/**********************/

	/******* Logging ******/
	/*CPXFILEptr logfile = nullptr;
	logfile = CPXXfopen("BPGC.log", "a");
	CPXXsetlogfile(env, logfile);*/
	/**********************/

	return status;
}

int edge_var_number(const Problem* prob, int schedule, int from, int to)
{
	assert(0 <= schedule && schedule < prob->schedules);
	assert(from != to);
	//This is a picewise function that maps the edge from 'from' to 'to', taking into consideration
	//that the graph doesn't have loop edges.
	int edge_number = schedule*prob->N*(prob->N-1) + (prob->N-1)*from + to; //schedule_offset + row_offset + col_offset
	if(to > from) edge_number--; //take into consideration loop edges not labeled/numbered.

	return edge_number;
}

int every_other_day_var_number(const Problem* prob, int schedule, int node)
{
	assert(0 <= schedule && schedule < prob->schedules);
	assert(node >= prob->depots+prob->every_day); //TODO: ojo con los casos bordes de los offsets aca
	return prob->schedules*prob->N*(prob->N-1) + prob->every_other_day*schedule + (node - prob->depots - prob->every_day); //edge_vars + schedule_offset + farm_node - (depots+every_day)
}

/*static int CPXPUBLIC subtour_constraint_generator(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	int status = 0;

	*useraction_p = CPX_CALLBACK_DEFAULT;
	auto* lazyconinfo = static_cast<lazy_constrain_info*>(cbhandle);

	double rhs;
	int cutnz = 0;
	int cur_numcols = CPXXgetnumcols(env, lazyconinfo->lp);
	auto* cutind = new int[cur_numcols];
	auto* cutval = new double[cur_numcols];
	auto* x = new double[cur_numcols];
	status = CPXXgetcallbacknodex(env, cbdata, wherefrom, x, 0, cur_numcols-1);
	if(status)
	{
		std::cerr << "Failed to get node solution." << std::endl;
		return (-1);
	}

	//Search for subtours in integer solution
	//std::vector<int> visited_among_subtours(lazyconinfo->prob->schedules,0);
	for(int day = 0; day < lazyconinfo->prob->schedules; ++day)
	{
		int visited_all_subtours = 0;
		for(int i = 0; i < lazyconinfo->prob->N; ++i)
		{
			if
		}

		int visited_depot_subtour;
		int cur_farm = 0;
		do
		{
			visited_depot_subtour++;
			int j = 0; while(j == cur_farm || x[edge_var_number(lazyconinfo->prob, day, cur_farm, j)] != 1.0) ++j;
			cur_farm = j;
		} while(cur_farm != 0);
	}



	status = CPXXcutcallbackadd(env, cbdata, wherefrom, cutnz, rhs, 'L', cutind, cutval, CPX_USECUT_PURGE);
	if(status) {
		std::cerr << "Failed to add cut." << std::endl;
		return (-1);

	} else {
		*//* Tell CPLEX that cuts have been created *//*
		*useraction_p = CPX_CALLBACK_SET;
	}

	return (status);
}*/

int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp)
{
	//Problem is minimization
	CPXXchgobjsen(env, lp, CPX_MIN);
	int status;

	//Model's variables definition
	std::cout << "-> Defining columns/variables" << std::endl;

	int edge_vars_quantity = prob->schedules * prob->N * (prob->N-1); //complete directed graph, for prob->schedules schedules
	int every_other_day_vars_quantity = prob->schedules * prob->every_other_day;
	int vars_quantity = edge_vars_quantity + every_other_day_vars_quantity;

	auto* obj_var_coefs = new double[vars_quantity];
	auto* lower_bounds = new double[vars_quantity];
	auto* upper_bounds = new double[vars_quantity];
	auto* column_types = new char[vars_quantity];

	int def_vars = 0;
	for(int i = 0; i < prob->N; ++i)
	{
		for(int day = 0; day < prob->schedules; ++day)
		{
			for(int j = 0; j < prob->N; ++j)
			{
				if(i != j)
				{
					obj_var_coefs[edge_var_number(prob, day, i, j)] = prob->costs[i][j];
					lower_bounds[edge_var_number(prob, day, i, j)] = 0.0;
					upper_bounds[edge_var_number(prob, day, i, j)] = 1.0;
					column_types[edge_var_number(prob, day, i, j)] = 'B';
					def_vars++;
				}
			}

			if(i >= prob->depots+prob->every_day)
			{
				obj_var_coefs[every_other_day_var_number(prob, day, i)] = 0.0;
				lower_bounds[every_other_day_var_number(prob, day, i)] = 0.0;
				upper_bounds[every_other_day_var_number(prob, day, i)] = 1.0;
				column_types[every_other_day_var_number(prob, day, i)] = 'B';
				def_vars++;
			}
		}
	}

	assert(def_vars == vars_quantity);
	status = CPXXnewcols(env, lp, vars_quantity, obj_var_coefs, lower_bounds, upper_bounds, column_types, nullptr);
	if(status) return (-1);

	//Model's restrictions
	std::cout << "-> Adding restrictions" << std::endl;
	auto* rhs = new double[1];
	auto* sense = new char[1];
	auto* rmatbeg = new CPXNNZ[1];
	auto* rmatind = new int[vars_quantity];
	auto* rmatind2 = new int[vars_quantity];
	auto* rmatval = new double[vars_quantity];

	//Every day client is visited exactly one time in each schedule
	// Sum_i ek_ij = 1, where k is the schedule number
	// Sum_i ek_ji = 1
	std::cout << "--> Every day farm is visited every day" << std::endl;
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = 1.0;
	sense[0] = 'E';
	for(int j = prob->depots; j < prob->depots + prob->every_day; ++j)
	{
		for(int day = 0; day < prob->schedules; ++day)
		{
			for(int i = 0, k = 0; i < prob->N; ++i)
			{
				if(i != j)
				{
					rmatval[k] = 1.0; //Variable coef
					rmatind[k] = edge_var_number(prob, day, i, j); //Variable number
					rmatind2[k] = edge_var_number(prob, day, j, i); //Variable number
					k++;
				}
			}

			//Add restrictions one by one
			status = CPXXaddrows(env, lp, 0, 1, prob->N-1, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status) return (-1);

			status = CPXXaddrows(env, lp, 0, 1, prob->N-1, rhs, sense, rmatbeg, rmatind2, rmatval, nullptr, nullptr);
			if(status) return (-1);
		}
	}

	//Every other day client is visited exactly one time in only one schedule
	// Sum_i e1_ij + e2_ij = 1
	// Sum_i e1_ji + e2_ji = 1
	std::cout << "--> Every other day farm is visited in only of the two day schedules" << std::endl;
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = 1.0;
	sense[0] = 'E';
	for(int j = prob->depots+prob->every_day; j < prob->N; ++j)
	{
		for(int day = 0, k = 0; day < prob->schedules; ++day)
		{
			for(int i = 0; i < prob->N; ++i)
			{
				if(i != j)
				{
					rmatval[k] = 1.0; //Variable coef
					rmatind[k] = edge_var_number(prob, day, i, j); //Variable number
					rmatind2[k] = edge_var_number(prob, day, j, i); //Variable number
					k++;
				}
			}
		}
		//Add restrictions one by one
		status = CPXXaddrows(env, lp, 0, 1, prob->schedules*(prob->N-1), rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status) return (-1);

		status = CPXXaddrows(env, lp, 0, 1, prob->schedules*(prob->N-1), rhs, sense, rmatbeg, rmatind2, rmatval, nullptr, nullptr);
		if(status) return (-1);
	}

	// Sum_i e1_ij - Sum_i e1_ji = 0
	std::cout << "---> Same schedule day arrival and departure for every other day farm" << std::endl;
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = 0.;
	sense[0] = 'E';
	for(int j = prob->depots+prob->every_day; j < prob->N; ++j)
	{
		for(int day = 0; day < prob->schedules; ++day)
		{
			for(int i = 0, k = 0; i < prob->N; ++i)
			{
				if(i != j)
				{
					rmatval[k] = 1.0; //Variable coef
					rmatind[k] = edge_var_number(prob, day, i, j); //Variable number

					rmatval[k+1] = -1.0; //Variable coef
					rmatind[k+1] = edge_var_number(prob, day, j, i); //Variable number
					k+=2;
				}
			}
			//Add restrictions one by one
			status = CPXXaddrows(env, lp, 0, 1, (prob->N-1)*2, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status) return (-1);
		}
	}

	//Vehicle capacity constraint
	// Sum_i,j ek_ij dj <= 80, where k = schedule number
	std::cout << "--> Vehicle capacity is not violated" << std::endl;
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = 80.;
	sense[0] = 'L';
	for(int day = 0; day < prob->schedules; ++day)
	{
		for(int i = 0, k = 0; i < prob->N; ++i)
		{
			for(int j = 0; j < prob->N; ++j)
			{
				if(i != j)
				{
					rmatval[k] = prob->demands[j]; //Variable coef
					rmatind[k] = edge_var_number(prob, day, i, j); //Variable number
					k++;
				}
			}
		}
		//Add restrictions one by one
		status = CPXXaddrows(env, lp, 0, 1, prob->N*(prob->N-1), rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status) return (-1);
	}

	//Depot belongs to each schedule
	// Sum_j ek_1j = 1, where k = schedule number
	// Sum_j ek_j1 = 1
	std::cout << "--> The deposit is part of every schedule" << std::endl;
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = 1.0;
	sense[0] = 'E';
	for(int day = 0; day < prob->schedules; ++day)
	{
		for(int j = 1; j < prob->N; ++j)
		{
			rmatval[j-1] = 1.0;
			rmatind[j-1] = edge_var_number(prob, day, 0, j);
			rmatind2[j-1] = edge_var_number(prob, day, j, 0);
		}
		//Add restrictions one by one
		status = CPXXaddrows(env, lp, 0, 1, prob->N-1, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status) return (-1);

		status = CPXXaddrows(env, lp, 0, 1, prob->N-1, rhs, sense, rmatbeg, rmatind2, rmatval, nullptr, nullptr);
		if(status) return (-1);
	}

	//Subtours restrictions
	std::cout << "--> No subtours in each schedule" << std::endl;

	//Consistency between edges selected and selected day for every other day clients.
	//sum_{i,j}(e1_ij) = 1+ every_day + sum_k(1-u_k) for k in every_other_day
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = static_cast<double>(1 + prob->every_day); //Every schedule uses exactly 1 depot
	sense[0] = 'E';
	for(int day = 0; day < prob->schedules; ++day)
	{
		int added = 0;
		for(int i = 0; i < prob->N; ++i)
		{
			for(int j = 0; j < prob->N; ++j)
			{
				if(i != j)
				{
					rmatval[added] = 1.0; //Variable coef
					rmatind[added] = edge_var_number(prob, day, i, j); //Variable number
					added++;
				}
			}

			if(i >= prob->depots+prob->every_day)
			{
				rmatval[added] = -1.0;
				rmatind[added] = every_other_day_var_number(prob, day, i);
				added++;
			}
		}

		//Add restrictions one by one
		assert(prob->N*(prob->N-1) + (prob->N - prob->depots - prob->every_day) == added);
		status = CPXXaddrows(env, lp, 0, 1, added, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
		if(status) return (-1);
	}

	//If U_ik says that i is visited during day k, then schedule must arrive till i during day k
	//U_j1 = Sum_{i} e1_ij
	rmatbeg[0] = 0; //Just one restriction
	rhs[0] = 0.0;
	sense[0] = 'E';
	for(int day = 0; day < prob->schedules; ++day)
	{
		for(int j = prob->depots+prob->every_day; j < prob->N; ++j)
		{
			int added = 0;
			for(int i = 0; i<prob->N; ++i)
			{
				if(j != i)
				{
					rmatval[added] = -1.0;
					rmatind[added] = edge_var_number(prob, day, i, j);
					added++;
				}
			}
			rmatval[added] = 1.0;
			rmatind[added] = every_other_day_var_number(prob, day, j);
			added++;

			//Add restrictions one by one
			assert((prob->N-1) + 1 == added);
			status = CPXXaddrows(env, lp, 0, 1, added, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status) return (-1);
		}
	}


	/* Set up to use MIP lazyconstraint callback. */
	/*auto* lazyconinfo = new lazy_constrain_info(lp, prob);
	status = CPXXsetlazyconstraintcallbackfunc(env, subtour_constraint_generator, lazyconinfo);
	if(status) return (-1);*/

	std::cout << "--> No more restrictions" << std::endl;

	// Este sector de codigo permite proveerle a CPLEX una solucion factible de entrada
	/*status = CPXsetintparam(env, CPX_PARAM_ADVIND, 1);
	if(status)exit(-1);


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
	delete[] rmatind2;
	delete[] rmatval;

	return 0;
}

/**********************/
double solve(Problem* prob, std::vector<std::vector<int>>& schedules)
{
	// CPlex environment
	int status = 0;
	CPXENVptr env = nullptr;
	CPXLPptr lp = nullptr;

	status = initialize_structures(env, lp);
	//std::cout << env << std::endl;
	if(status != 0)
	{
		free_structures(env, lp);
		return status;
	}

	status = set_parameters(env);
	if(status)
	{
		std::cout<<"Error seteando parametros"<<std::endl;
		free_structures(env, lp);
		return status;
	}

	std::cout << "Initializing MIP" << std::endl;
	status = initialize_mip(prob, env, lp);
	if(status)
	{
		free_structures(env, lp);
		return status;
	}

	// DEBUG.
	/*status = CPXXwriteprob(env, lp, "1e.lp", nullptr);
	if(status)
	{
		free_structures(env, lp);
		return status;
	}*/

	//Solve MIP //TODO: ver de pasar a su propia funcion la parte donde se optimiza y construye el resultado
	std::cout << "Solving" << std::endl;
	status = CPXXmipopt(env,lp);
	if(status) return(-1);

	//Retrieve optimal solution
	double objval;
	//double* y = nullptr;
	status = CPXXgetobjval (env, lp, &objval); // Get optimum objective function value
	if(status) return(-1);

	//int cur_numrows = CPXXgetnumrows(env,lp);
	int cur_numcols = CPXXgetnumcols(env,lp);
	auto* vars = new double[cur_numcols];
	status = CPXXgetx (env, lp, vars, 0, cur_numcols-1); // Pedimos los valores de las variables en la solucion
	if(status) return(-1);

	//Armamos el orden de cada schedule respecto de la nomeclatura de las variables
	schedules.clear();
	for(int day = 0; day < prob->schedules; ++day)
	{
		schedules.emplace_back();
		int cur_farm = 0;
		do
		{
			schedules[day].push_back(cur_farm);
			int j = 0; while(j == cur_farm || vars[edge_var_number(prob, day, cur_farm, j)] != 1.0) ++j;
			cur_farm = j;
		} while(cur_farm != 0);
		schedules[day].push_back(cur_farm);
	}

	delete[] vars;
	free_structures(env, lp);

	return objval;
}