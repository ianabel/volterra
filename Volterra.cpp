
#include "Volterra.hpp"

#include <utility>

TrapezoidSolver::TrapezoidSolver() {};

double TrapezoidSolver::Fn( unsigned int n )
{
	if( n == 0 )
		return g(0);

	double tn = n*h;
	double x = g(tn) + (h/2.0) * ( K( tn, 0 )*F[0] );

	for( unsigned int i = 1; i < n; ++i )
	{
		double ti = i*h;
		x += h*K(tn,ti)*F[i];
	}

	// x = Fn - (h/2)*K(tn,tn)*Fn
	
	return x / (1.0 - (h/2.0)*K(tn,tn));
}


Vector&& TrapezoidSolver::Solve ( SourceFunction f, KernelFunction Kfn, double t, unsigned int N ) 
{
	h = t/static_cast<double>(N);
	K = Kfn;
	g = f;
	F.resize( N + 1 );

	for( unsigned int i = 0; i < N + 1; ++i )
	{
		F[i] = Fn( i );
	}

	return std::move( F );
}
