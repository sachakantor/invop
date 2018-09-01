#ifndef INVOP_1A_PROBLEM_HPP
#define INVOP_1A_PROBLEM_HPP

#include<vector>

#include<ilcplex/cplexx.h>


struct Problem {
	int V; //number of vegetable oils
	double max_refined_per_month_V; //max tons of vegetable oils refined per month
	int NV; //number of non-vegetable oils
	double max_refined_per_month_NV; //max tons of vegetable oils refined per month
	int M; //number of months
	double max_hardness, min_hardness;
	int max_oils_per_month;

	//TODO: cambiar de vectores sincronizados a un solo vector de structs acordes
	std::vector<double> hardness;
	std::vector<std::vector<double>> oil_buy_cost_per_month;
	std::vector<std::vector<double>> oil_storage_cost_per_month;
	std::vector<std::vector<double>> oil_profit_per_month;
	std::vector<double> storage_capacity;
	std::vector<double> min_tons_oil_if_refined;
	std::vector<double> initial_storage;

	Problem();
	//explicit Problem(int n);
};

struct month_oil_solution {
	double tons_refined;
	double tons_storaged;
	double tons_buyed;
	int month;
};

class mip_vars {
public:
	explicit mip_vars(const Problem* prob);
	enum class var_type { refined, storage, buy, oil_switch};
	int var_num(var_type t, int oil, int month) const;
	int total() const;
	int refined_total() const;
	int storage_total() const;
	int buy_total() const;
	int oil_use_switch_total() const;

private:
	const int oil_q;
	const int month_q;
	const int refined_base_offset;
	const int storage_base_offset;
	const int buy_base_offset;
	const int oil_use_switch_base_offset;
};

int initialize_structures(CPXENVptr& env, CPXLPptr& lp);
void free_structures(CPXENVptr env, CPXLPptr lp);
int set_parameters(CPXENVptr env);
int edge_var_number(const Problem* prob, int day, int from, int to);
int vertex_var_number(const Problem* prob, int day, int vertex);
int initialize_mip(Problem* prob, CPXENVptr env, CPXLPptr lp);
int solveMIP(const Problem* prob,
						 CPXENVptr env,
						 CPXLPptr lp,
						 double* const objval,
						 std::vector<std::vector<month_oil_solution>>& solution);
int solve(Problem* prob, double& objval, std::vector<std::vector<month_oil_solution>>& solution);

#endif //INVOP_1A_PROBLEM_HPP
