#include<iostream>
#include<vector>

#include<problem.hpp>
#include<ilcplex/cplexx.h>


int main(int argc, char* argv[])
{
	std::cout << "Hello world" << std::endl;
	return solve(new Problem());
}