
#include "Volterra.hpp"
#include "GK.hpp"

#include <cmath>
#include <iostream>
#include <fstream>


int main( int, char ** )
{
	LinearGK lin_gk( 1.0, 50.0, 2.0, 0.0, 1.0, 1.0 );

	std::cout << "# lambda\t Gamma0(lambda,lambda)" << std::endl;
	
	std::vector<double> lambda{ 1.0, 4.0, 64.0, 256.0, 1024.0, 4096.0, 8192.0 };
	for( auto l : lambda )
	{
		std::cout << l << '\t' << lin_gk.Gamma0( l, l ) << std::endl;
	}

	std::cout << "# lambda\tlambda'\t Gamma0(lambda,lambda')" << std::endl;
	std::vector< std::pair< double, double > > lambdas{ {746,8}, {750,10} };
	for( auto pa : lambdas )
	{
		std::cout << pa.first << '\t' << pa.second << '\t' << lin_gk.Gamma0( pa.first, pa.second ) << std::endl;
	}
	return 0;
}
