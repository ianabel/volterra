
#include "Volterra.hpp"

#include <cmath>
#include <iostream>

double Source( double t ) { return std::exp( t ); };
double Kernel( double t, double s ) { return -std::exp( t - s ); };

double Source2( double t ) { return 1; };
double Kernel2( double t, double s ) { return std::sin(t-s); };

int main( int, char ** )
{
	TrapezoidSolver<double> tsolve;
	std::vector<double> data;
	unsigned int N = 10;

	data = tsolve.Solve( Source, Kernel, 1.0, N );
	for( unsigned int i=0; i <= N; ++i ) 
	{
		std::cout << i*1.0/N << '\t' << 1.0 - data[i] << std::endl;
	}
	std::cout << std::endl;

	data = tsolve.Solve( Source2, Kernel2, 1.0, N );
	for( unsigned int i=0; i <= N; ++i ) 
	{
		double ti = i * 1.0 / N;
		std::cout << ti << '\t' << 1.0 + (ti*ti/2.0) - data[i] << std::endl;
	}

	std::cout << std::endl;
	SimpsonSolver<double> ssolve;
	data = ssolve.Solve( Source, Kernel, 1.0, N );
	for( unsigned int i=0; i <= N; ++i ) 
	{
		double ti = i * 1.0 / N;
		std::cout << ti << '\t' << 1.0 - data[i] << std::endl;
	}

	std::cout << std::endl;
	data = ssolve.Solve( Source2, Kernel2, 1.0, N );
	for( unsigned int i=0; i <= N; ++i ) 
	{
		double ti = i * 1.0 / N;
		std::cout << ti << '\t' << 1.0 + (ti*ti/2.0) - data[i] << std::endl;
	}
	return 0;
}
