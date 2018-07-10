#ifndef INVOP_1E_PROBLEM_HPP
#define INVOP_1E_PROBLEM_HPP

#include<vector>

#include<ilcplex/cplexx.h>


struct Problem
{
	int N;  				//number of clients (including depot)
	int K;  				//vehicle capacity
	int schedules;	//number of schedules
	int every_day;	//max client number who must be visited every day
	std::vector<std::vector<int>> costs;
	std::vector<int> demands;

	Problem();
};

int initialize_structures(CPXENVptr& env, CPXLPptr& lp);
void free_structures(CPXENVptr env, CPXLPptr lp);
int set_parameters(CPXENVptr env);
int edge_var_number(const Problem* prob, int schedule, int from, int to);
int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp);
int solve(Problem* prob);

#endif //INVOP_PROBLEM_HPP
