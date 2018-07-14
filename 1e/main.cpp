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
		int n, depots, every_day_clients, every_other_day_clients;
		double vehicle_capacity;
		input_file >> n >> depots >> every_day_clients >> every_other_day_clients >> vehicle_capacity;

		assert(n == depots+every_day_clients+every_other_day_clients);
		//TODO: agergar assert de factibilidad de solucion SIN multitrip

		auto* prob = new Problem(n);
		prob->N = n;
		prob->depots = depots;
		prob->every_day = every_day_clients;
		prob->cust_clients = every_other_day_clients;
		prob->cust_clients_freq = 1;
		prob->K = vehicle_capacity;
		prob->schedules = 2; //TODO: Esto está hardcodeado, debería ser pasado a las instancias de entrada
		for(int i=0; i<depots; ++i) prob->demands[i] = 0.0;
		for(int i=depots; i<n; ++i) input_file >> prob->demands[i];
		for(int i=0; i<n; ++i)
		{
			for(int j=0; j<n; ++j)
			{
				input_file >> prob->costs[i][j];
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
		std::cout << "Costo solucion optima: " << cost << std::endl;
		for(int day = 0; day < prob->schedules; ++day)
		{
			std::cout << "Ruta día " << day << ": ";
			for(auto i : schedules[day]) std::cout << i+1 << " ";
			std::cout << std::endl;
		}

		//End
		delete prob;
	}

	return 0;
}