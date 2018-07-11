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
		int n, depots, every_day_clients;
		double vehicle_capacity;
		input_file >> n >> depots >> every_day_clients >> vehicle_capacity;

		auto* prob = new Problem(n);
		prob->N = n;
		prob->depots = depots;
		prob->every_day = every_day_clients;
		prob->K = vehicle_capacity;
		prob->schedules = 2; //Esto está hardcodeado, debería ser pasado a las instancias de entrada
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
		solve(prob);

		delete prob;
	}

	return 0;
}