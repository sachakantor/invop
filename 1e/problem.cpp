#include<cassert>
#include<iostream>
#include<vector>

#include<problem.hpp>
#include<ilcplex/cplexx.h>

Problem::Problem() = default;

Problem::Problem(int n) : costs(n, std::vector<double>(n,0.0)), demands(n,0.0) {}

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

	/******* Logging ******/
	/*CPXFILEptr logfile = nullptr;
	logfile = CPXXfopen("BPGC.log", "a");
	CPXXsetlogfile(env, logfile);*/
	/**********************/

	return status;
}

int edge_var_number(const Problem* prob, int schedule, int from, int to)
{
	assert(from != to);
	//This is a picewise function that maps the edge from 'from' to 'to', taking into consideration
	//that the graph doesn't have loop edges.
	int edge_number = schedule*prob->N*(prob->N-1) + (prob->N-1)*from + to; //schedule_offset + row_offset + col_offset
	if(to > from) edge_number--; //take into consideration loop edges not labeled/numbered.

	return edge_number;
}

int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp)
{
	//Problem is minimization
	CPXXchgobjsen(env, lp, CPX_MIN);
	int status;

	std::cout << "Definiendo columnas/variables" << std::endl;
	//Model's variables definition
	int var_number = prob->schedules * prob->N * (prob->N-1); //complete directed graph, for prob->schedules schedules
	auto* obj_var_coefs = new double[var_number];
	auto* lower_bounds = new double[var_number];
	auto* upper_bounds = new double[var_number];
	auto* column_types = new char[var_number];

	for(int i = 0; i < prob->N; ++i)
	{
		for(int j = 0; j < prob->N; ++j)
		{
			if(i != j)
			{
				for(int day = 0; day < prob->schedules; ++day)
				{
					obj_var_coefs[edge_var_number(prob, day, i, j)] = prob->costs[i][j];
					lower_bounds[edge_var_number(prob, day, i, j)] = 0.0;
					upper_bounds[edge_var_number(prob, day, i, j)] = 1.0;
					column_types[edge_var_number(prob, day, i, j)] = 'B';
				}
			}
		}
	}

	status = CPXXnewcols(env, lp, var_number, obj_var_coefs, lower_bounds, upper_bounds, column_types, nullptr);
	if(status) return (-1);

	std::cout << "Agregando restricciones" << std::endl;

	//Model's restrictions
	auto* rhs = new double[1];
	auto* sense = new char[1];
	auto* rmatbeg = new CPXNNZ[1];
	auto* rmatind = new int[var_number];
	auto* rmatind2 = new int[var_number];
	auto* rmatval = new double[var_number];

	std::cout << "-> Cada cliente que debe ser visitado todos los dias lo es." << std::endl;
	//Every day client is visited exactly one time in each schedule
	// Sum_i ek_ij = 1, where k is the schedule number
	// Sum_i ek_ji = 1
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

	std::cout << "-> Cada cliente que debe ser visitado UN día, lo es." << std::endl;
	//Every other day client is visited exactly one time in only one schedule
	// Sum_i e1_ij + e2_ij = 1
	// Sum_i e1_ji + e2_ji = 1
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

	std::cout << "--> Se entra y sale el mismo día" << std::endl;
	// Sum_i e1_ij - Sum_i e1_ji = 0
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

	std::cout << "-> Capacidad del vehículo" << std::endl;
	//Vehicle capacity constraint
	// Sum_i,j ek_ij dj <= 80, where k = schedule number
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

	std::cout << "-> Deposito está en todos los schedules" << std::endl;
	//Depot belongs to each schedule
	// Sum_j ek_1j = 1, where k = schedule number
	// Sum_j ek_j1 = 1
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

	std::cout << "-> Fin restricciones" << std::endl;
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

int solve(Problem* prob)
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

	std::cout << "Inicializando MIP" << std::endl;

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

	std::cout << "Resolver" << std::endl;
	// Variables para traer los resultados
	double objval;
	//double* y = nullptr;

	//Solve MIP
	status = CPXXmipopt(env,lp);
	if(status) return(-1);

	status = CPXXgetobjval (env, lp, &objval); // Get optimum objective function value
	if(status) return(-1);

	// Pido memoria para traer la solucion
	//int cur_numrows = CPXXgetnumrows(env,lp);
	//int cur_numcols = CPXXgetnumcols(env,lp);

	//double* x = new double[cur_numcols];

	/*// Pedimos los valores de las variables en la solucion
	status = CPXXgetx (env, lp, x, 0, cur_numcols-1);
	if(status) return(-1);

	// Armamos la asignacion de colores en nuestra estructura con la informacion que nos trajimos de CPLEX
	if(prob->bestObj>=objval){
		for(int i=0; i<prob->N; ++i){
			for(int j=0; j<prob->bestObj; ++j){
				if(x[i*prob->bestObj+j]>0.001){
					prob->colores[i]=j;
				}
			}
		}
		prob->bestObj = objval;
	}*/

	free_structures(env, lp);

	return status;
}