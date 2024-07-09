
#include "Volterra.hpp"
#include "GK.hpp"

#include <cmath>
#include <iostream>
#include <fstream>


int main( int, char ** )
{
	unsigned int N = 1024.0;
	double tEnd = 30.0;
	SimpsonSolver< std::complex<double> > ssolve;
	LinearGK lin_gk( 1.0, 50.0, 2.0, 0.0, 1.0, 1.0 );

	auto Source = lin_gk.getSourceFn();
	auto Kernel = lin_gk.getKernelFn();
	auto data = ssolve.Solve( Source, Kernel, tEnd, N );

	std::cout << "# t\tPhi^2" << std::endl;
	
	for( unsigned int i=0; i <= N; ++i )
	{
		double ti = i * tEnd / N;
		std::cout << ti << '\t' << std::norm(data[i]) << std::endl;
	}
	return 0;
}
