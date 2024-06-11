
#include "Volterra.hpp"

#include <cmath>
#include <iostream>

double Source( double t ) { return std::exp( t ); };
double Kernel( double t, double s ) { return -std::exp( t - s ); };


int main( int, char ** )
{
	TrapezoidSolver tsolve;
	std::vector<double> data;
	unsigned int N = 10;
	data = tsolve.Solve( Source, Kernel, 1.0, N );
	for( unsigned int i=0; i <= N; ++i ) 
	{
		std::cout << i*1.0/N << '\t' << 1.0 - data[i] << std::endl;
	}
	return 0;
}
