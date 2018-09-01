#include<cmath>

#include<coordinates.hpp>

coordinate::coordinate() = default;

coordinate::coordinate(double x, double y) : x(x), y(y) {}

double coordinate::distance(const coordinate& c)
{
	return std::sqrt(pow(x-c.x, 2) + pow(y-c.y, 2));
}