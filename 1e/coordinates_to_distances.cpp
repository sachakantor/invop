#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>

#include<coordinates.hpp>

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cerr << "El formato de ejecucion debe 1 archivo de entrada y otro de salida:" << std::endl;
		std::cerr << argv[0] << " coordinates.in distances.out" << std::endl;
		return 1;
	}

	std::ifstream input_file(argv[1]);
	if(input_file.fail())
	{
		std::cerr << "No se pudo abrir el archivo \"" << argv[1] << "\"" << std::endl;
		return -1;
	};

	std::ofstream output_file(argv[2]);
	if(output_file.fail())
	{
		std::cerr << "No se pudo abrir el archivo \"" << argv[2] << "\"" << std::endl;
		return -1;
	};

	//Read input
	int n, depots, every_day;
	std::vector<coordinate> coordinates;
	std::vector<double> demands;

	input_file >> n >> depots >> every_day;

	for(int i = 0; i < depots; ++i)
	{
		double x, y;
		input_file >> x >> y;
		coordinates.emplace_back(x, y);
		demands.emplace_back(0.0);
	}

	for(int i = depots; i < n; ++i)
	{
		double x, y, d;
		input_file >> x >> y >> d;
		coordinates.emplace_back(x, y);
		demands.emplace_back(d);
	}

	//Calculate euclidean distances
	std::vector<std::vector<double>> distances(n, std::vector<double>(n, 0));
	for(int i = 0; i < n; ++i)
	{
		for(int j = i+1; j < n; ++j)
		{
			distances[i][j] = coordinates[i].distance(coordinates[j]);
			distances[j][i] = distances[i][j];
		}
	}

	//save output
	output_file << n << " " << depots << " " << every_day << std::endl;
	output_file << demands[depots];
	for(int i=depots+1; i<n; ++i)
		output_file << " " << demands[i];
	output_file << std::endl;

	for(int i=0; i<n; ++i)
	{
		for(auto distance : distances[i]) output_file << distance << " ";
		output_file << std::endl;
	}

	return 0;
}