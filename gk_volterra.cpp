
#include "Volterra.hpp"
#include "GK.hpp"

#include <cmath>
#include <iostream>
#include <fstream>

using std::sqrt;

int main( int, char ** )
{
	unsigned int N = 2048;
	double tEnd = 200.0 * sqrt(2.0);
	SimpsonSolver< std::complex<double> > ssolve;
	SlabLinearGK lin_gk( 1.0, 0.05, 5.0/sqrt(2), 2.0, 0.5*sqrt(2) );

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
