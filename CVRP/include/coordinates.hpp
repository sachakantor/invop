#ifndef COORDINATES_HPP
#define COORDINATES_HPP

struct coordinate
{
	double x, y;

	coordinate();
	coordinate(double x, double y);

	double distance(const coordinate&);
};


#endif //COORDINATES_HPP
