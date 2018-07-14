#include<algorithm>
#include<cassert>
#include<iomanip>
#include<iostream>
#include<fstream>
#include<string>

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
		std::vector<std::vector<month_oil_solution>> solution;
		double cost;
		int status;
		status = solve(prob, cost, solution);
		if(status != 0)
		{
			return status;
		}

		//Print solution
		std::cout << "Costo solucion optima: " << cost << std::endl;

		//Header //TODO: Todo lo que viene es demasiado casero. Habría que hacerlo bien.
		std::cout << std::left << std::setw(3) << "Mes";
		std::string underline = "===";
		for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
			char oil_type = 'V';
			int oil_num = oil+1;
			if(oil >= prob->V) {
				oil_type = 'O';
				oil_num -= prob->V;
			}
			std::cout << std::setw(9) << " || Ref " << std::setw(1) << oil_type << std::setw(2) << oil_num
				<< std::setw(7) << " | Alm " << std::setw(1) << oil_type << std::setw(2) << oil_num
				<< std::setw(7) << " | Com " << std::setw(1) << oil_type << std::setw(2) << oil_num;

			underline += "=##==========|=========|========";
		}
		std::cout << std::endl << underline << std::endl;

		//Values
		std::replace(underline.begin(),underline.end(), '=', '-');
		std::replace(underline.begin(),underline.end(), '#', '+');
		for(int month = 0; month < prob->M; ++month) {
			std::cout << std::left << std::setw(3) << month+1;
			for(int oil = 0; oil < prob->V+prob->NV; ++oil) {
				std::cout << std::setw(4) << " || " << std::right << std::setw(8) << solution[month][oil].tons_refined
					<< std::setw(3) << " | " << std::right << std::setw(7) << solution[month][oil].tons_storaged
					<< std::setw(3) << " | " << std::right << std::setw(7) << solution[month][oil].tons_buyed;
			}
			std::cout << std::endl << underline << std::endl;
		}

		//End
		delete prob;
	}

	return 0;
}