#include<cassert>
#include<iostream>
#include<fstream>

#include<problem.hpp>
#include<ilcplex/cplexx.h>


int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cerr << "El formato de ejecucion debe recibir al menos 1 parámetro:" << std::endl;
		std::cerr << argv[0] << " input1.in [ input2.in input3.in  ... ]" << std::endl;
		return 1;
	}

	for(int k = 1; k < argc; ++k)
	{
		std::ifstream input_file(argv[k]);
		if(input_file.fail())
		{
			std::cerr << "No se pudo abrir el archivo \"" << argv[k] << "\"" << std::endl;
			break;
		};

		//Load instance
		int veg_q, non_veg_q, months_q, max_oils_per_month;
		double max_veg_tons_per_month, max_non_veg_tons_per_month, min_hardness, max_hardness;
		input_file
			>> veg_q
			>> max_veg_tons_per_month
			>> non_veg_q
			>> max_non_veg_tons_per_month
			>> months_q
			>> min_hardness
			>> max_hardness
			>> max_oils_per_month;

		auto* prob = new Problem();
		prob->V = veg_q;
		prob->max_refined_per_month_V = max_veg_tons_per_month;
		prob->NV = non_veg_q;
		prob->max_refined_per_month_NV = max_non_veg_tons_per_month;
		prob->M = months_q;
		prob->min_hardness = min_hardness;
		prob->max_hardness = max_hardness;
		prob->max_oils_per_month = max_oils_per_month;

		for(int oil = 0; oil < veg_q+non_veg_q; ++oil) {
			double hardness;
			input_file >> hardness;
			prob->hardness.push_back(hardness);

			double storage_capacity;
			input_file >> storage_capacity;
			prob->storage_capacity.push_back(storage_capacity);

			double min_tons_if_refined;
			input_file >> min_tons_if_refined;
			prob->min_tons_oil_if_refined.push_back(min_tons_if_refined);

			double initial_storage;
			input_file >> initial_storage;
			prob->initial_storage.push_back(initial_storage);

			prob->oil_buy_cost_per_month.emplace_back();
			for(int m = 0; m < months_q; ++m) {
				double cost;
				input_file >> cost;
				prob->oil_buy_cost_per_month.back().push_back(cost);
			}

			prob->oil_storage_cost_per_month.emplace_back();
			for(int m = 0; m < months_q; ++m) {
				double cost;
				input_file >> cost;
				prob->oil_storage_cost_per_month.back().push_back(cost);
			}

			prob->oil_profit_per_month.emplace_back();
			for(int m = 0; m < months_q; ++m) {
				double profit;
				input_file >> profit;
				prob->oil_profit_per_month.back().push_back(profit);
			}
		}

		//Solve
		std::vector<std::vector<int>> schedules;
		double cost;
		int status;
		status = solve(prob, schedules, cost);
		if(status != 0)
		{
			return status;
		}

		//Print solution
		/*std::cout << "Costo solucion optima: " << cost << std::endl;
		for(int day = 0; day < prob->schedules; ++day)
		{
			std::cout << "Ruta día " << day << ": ";
			for(auto i : schedules[day]) std::cout << i+1 << " ";
			std::cout << std::endl;
		}*/

		//End
		delete prob;
	}

	return 0;
}