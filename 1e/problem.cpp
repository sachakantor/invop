#include<cassert>
#include<iostream>
#include<vector>

#include<problem.hpp>
#include<ilcplex/cplexx.h>


/**********************/
Problem::Problem() = default;
Problem::Problem(int n) : costs(n, std::vector<double>(n,0.0)), demands(n,0.0) {}

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
	status = CPXXsetintparam (env, CPXPARAM_Preprocessing_Linear, 0);
	if(status) return (-1);

	/* Turn on traditional search for use with control callbacks */
	status = CPXXsetintparam (env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL); //mepa que se solapa con otra
	if(status) return (-1);

	/* Let MIP callbacks work on the original model */
	status = CPXXsetintparam (env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
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

static int CPXPUBLIC subtour_constraint_generator(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	int status = 0;

	Problem* prob = static_cast<Problem*>(cbhandle);

	/*int      numcols  = cutinfo->numcols;
	int      numcuts  = cutinfo->num;
	double   *x       = cutinfo->x;
	int      *beg     = cutinfo->beg;
	int      *ind     = cutinfo->ind;
	double   *val     = cutinfo->val;
	double   *rhs     = cutinfo->rhs;
	int      *cutind  = nullptr;
	double   *cutval  = nullptr;
	double   cutvio;
	int      addedcuts = 0;
	int      i, j, k, cutnz;*/

	*useraction_p = CPX_CALLBACK_DEFAULT;

	status = CPXXgetcallbacknodex(env, cbdata, wherefrom, x, 0, numcols-1);
	if(status)
	{
		std::cerr << "Failed to get node solution." << std::endl;
		return (-1);
	}

	for (i = 0; i < numcuts; i++) {
		cutvio = -rhs[i];
		k = beg[i];
		cutnz = beg[i+1] - k;
		cutind = ind + k;
		cutval = val + k;
		for (j = 0; j < cutnz; j++) {
			cutvio += x[cutind[j]] * cutval[j];
		}

		/* Use a cut violation tolerance of 0.01 */
		if ( cutvio > 0.01 ) {
			status = CPXXcutcallbackadd (env, cbdata, wherefrom, cutnz, rhs[i], 'L', cutind, cutval, CPX_USECUT_PURGE);
			if(status) {
				fprintf (stderr, "Failed to add cut.");
				return (-1);
			}
			addedcuts++;
		}
	}

	/* Tell CPLEX that cuts have been created */
	if(addedcuts > 0) *useraction_p = CPX_CALLBACK_SET;

	return (status);
} /* END mycutcallback */

static int makelazyconstraint(CPXENVptr env, CPXLPptr lp, cutinfo* lazyconinfo)
{
	int status = 0;

	int beg[] = {0, 5};
	double val[] = {1, 1, 1, 1, 1};
	char* varname[] = { "W11", "W12", "W13", "W14", "W15" };
	double rhs[] = {3};

	int    *cutbeg = nullptr;
	int    *cutind = nullptr;
	double *cutval = nullptr;
	double *cutrhs = nullptr;

	int i, varind;
	int nz   = 5;
	int cuts = 1;

	int cur_numcols = CPXXgetnumcols(env, lp);

	lazyconinfo->lp = lp;
	lazyconinfo->numcols = cur_numcols;

	lazyconinfo->x = new double[cur_numcols];
	if (lazyconinfo->x == nullptr) {
		std::cerr << "No memory for solution values." << std::endl;
		return (-1);
	}

	cutbeg = (int *)    malloc ((cuts+1) * sizeof (int));
	cutind = (int *)    malloc (nz * sizeof (int));
	cutval = (double *) malloc (nz * sizeof (double));
	cutrhs = (double *) malloc (cuts * sizeof (double));

	if ( cutbeg == NULL ||
			 cutind == NULL ||
			 cutval == NULL ||
			 cutrhs == NULL   ) {
		fprintf (stderr, "No memory.\n");
		status = CPXERR_NO_MEMORY;
		goto TERMINATE;
	}

	for (i = 0; i < nz; i++) {
		status = CPXXgetcolindex(env, lp, varname[i], &varind);
		if ( status )  {
			fprintf (stderr,
							 "Failed to get index from variable name.\n");
			goto TERMINATE;
		}
		cutind[i] = varind;
		cutval[i] = val[i];
	}

	for (i = 0; i < cuts; i++) {
		cutbeg[i] = beg[i];
		cutrhs[i] = rhs[i];
	}
	cutbeg[cuts] = beg[cuts];

	lazyconinfo->num      = cuts;
	lazyconinfo->beg      = cutbeg;
	lazyconinfo->ind      = cutind;
	lazyconinfo->val      = cutval;
	lazyconinfo->rhs      = cutrhs;

	TERMINATE:

	if ( status ) {
		free_and_null ((char **) &cutbeg);
		free_and_null ((char **) &cutind);
		free_and_null ((char **) &cutval);
		free_and_null ((char **) &cutrhs);
	}

	return (status);
} /* END makelazyconstraint */

/* Init information on the node objval for the user cut callback */
/*static void initnodeobjvalinfo(CPXENVptr env, CPXLPptr lp, cutinfo* cutinfo)
{
	cutinfo->nodeid = -1;
	cutinfo->nodeobjval = 0.0;
	cutinfo->objsen = CPXgetobjsen (env, lp);
	if ( cutinfo->objsen == CPX_MIN )
		cutinfo->objsen = 1;
	else
		cutinfo->objsen = -1;
}*/ /* END initnodeobjvalinfo */

int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp)
{
	//Problem is minimization
	CPXXchgobjsen(env, lp, CPX_MIN);
	int status;

	//Model's variables definition
	std::cout << "-> Defining columns/variables" << std::endl;
	int vars_quantity = prob->schedules * prob->N * (prob->N-1); //complete directed graph, for prob->schedules schedules
	auto* obj_var_coefs = new double[vars_quantity];
	auto* lower_bounds = new double[vars_quantity];
	auto* upper_bounds = new double[vars_quantity];
	auto* column_types = new char[vars_quantity];

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
	/*rmatbeg[0] = 0; //Just one restriction
	sense[0] = 'L';
	for(int day = 0; day < prob->schedules; ++day)
	{
		std::cout << "Cantidad de variables" << vars_quantity/prob->schedules << std::endl;
		for(unsigned long long subset = 1, subset_size = 0; subset < (static_cast<unsigned long long>(1)<<(vars_quantity/prob->schedules))-1; ++subset, subset_size=0)
		{
			//std::cout << "Considerando subset: " << subset << std::endl;
			for(int var_number = 0; var_number < vars_quantity/prob->schedules; ++var_number)
			{
				if((subset & (static_cast<unsigned int>(1) << var_number)) > 0)
				{
					rmatval[subset_size] = 1.0;
					rmatind[subset_size] = day*prob->N*(prob->N-1) + var_number; //Estaría bueno que el comportamiento de day en la numeración de la variable esté encapsulado
					subset_size++;
				}
			}

			rhs[0] = subset_size - 1.0;

			//Add restrictions one by one
			status = CPXXaddrows(env, lp, 0, 1, subset_size, rhs, sense, rmatbeg, rmatind, rmatval, nullptr, nullptr);
			if(status) return (-1);
		}
	}*/

	/* Create a lazy constraint to alter the optimum */
	//cutinfo lazyconinfo; //TODO: esto ha de estar dentro de Problem y ser liberado luego de usarse
	//status = makelazyconstraint(env, lp, &lazyconinfo);
	//if(status) return (-1);

	/* Set up to use MIP lazyconstraint callback. The callback funtion registered is the same, but the data will be different. */
	//status = CPXXsetlazyconstraintcallbackfunc(env, subtour_constraint_generator, &lazyconinfo);
	status = CPXXsetlazyconstraintcallbackfunc(env, subtour_constraint_generator, prob);
	if(status) return (-1);
	//status = myoptimize (env, lp, -39.0);

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