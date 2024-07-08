
#include "Volterra.hpp"

#include <utility>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/Dense>


template<typename T> T TrapezoidSolver<T>::Fn( unsigned int n )
{
	if( n == 0 )
		return g(0);

	double tn = n*h;
	T x = g(tn) + (h/2.0) * ( K( tn, 0 )*F[0] );

	for( unsigned int i = 1; i < n; ++i )
	{
		double ti = i*h;
		x += h*K(tn,ti)*F[i];
	}

	// x = Fn - (h/2)*K(tn,tn)*Fn
	
	return x / (1.0 - (h/2.0)*K(tn,tn));
}

template<typename T> typename VolterraSolver<T>::Vector&& TrapezoidSolver<T>::Solve( typename VolterraSolver<T>::SourceFunction f, typename VolterraSolver<T>::KernelFunction Kfn, double t, unsigned int N )
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

double w_ni( unsigned int n, unsigned int i )
{
	if( n % 2 == 0 ) {
		if( i == 0 || i == n ) {
			return 1./3.;
		} else if ( i % 2 == 0 ) {
			return 2./3.;
		} else {
			return 4./3.;
		}
	} else {
		if( i == 0  ) {
			if( n == 3 ) {
				return 3./8.;
			} else {
				return 1./3.;
			}
		} else if ( i == n ) {
			return 3./8.;
		} else if ( i + 2 >= n ) {
			return 9./8.;
		} else if ( i + 3 == n ) {
			return 17./24.;
		} else if ( i % 2 == 0 ) {
			return 2./3.;
		} else {
			return 4./3.;
		}
	}
}

template<typename T> T SimpsonSolver<T>::Fn( unsigned int n )
{
	if( n <= 2 )
		throw std::logic_error("Should have pre-computed F[0],F[1],F[2]");

	double tn = n*h;
	T x = g(tn);

	for( unsigned int i = 0; i < n; ++i )
	{
		double ti = i*h;
		x += h*w_ni(n,i)*K(tn,ti)*F[i];
	}

	// x = Fn - (h/2)*w_nn*K(tn,tn)*Fn
	
	return x / (1.0 - h*w_ni(n,n)*K(tn,tn));
}


template<typename T> typename VolterraSolver<T>::Vector&& SimpsonSolver<T>::Solve( typename VolterraSolver<T>::SourceFunction f, typename VolterraSolver<T>::KernelFunction Kfn, double t, unsigned int N )
{
	h = t/static_cast<double>(N);
	K = Kfn;
	g = f;
	F.resize( N + 1 );

	if( N == 0 )
		throw std::logic_error("Cannot have N = 0!");
	if( N == 1 )
		throw std::logic_error("Higher order methods need more than one step");

	/*
		Use the method in Section 7.6 of Analytical and Numerical Methods for Volterra Equations. 1985, 95-128
		to solve for F[1] F[2]
	*/

	Eigen::Matrix< T, 2, 1> F12,G;
	Eigen::Matrix< T, 2, 2> A;

	G(0) = g( h ) + (h/6.0) * ( K( h, 0.0 ) + (3./2.)*K( h, h/2.0 ) ) * g( 0.0 );
	G(1) = g( 2.0 * h ) + (h/3.0)*K( 2.0*h, 0.0 )*g( 0.0 );
	A(0,0) = (1./2.)*( 3.0*K( h, h/2. ) + K(h,h) );
	A(0,1) = -(1./4.)*K( h, h/2. );
	A(1,0) = 4. * K( 2. * h, h );
	A(1,1) = K(2. * h, 2. * h);

	F12 = ( Eigen::Matrix<T,2,2>::Identity() - (h/3.)*A ).inverse() * G;

	F[0] = g( 0.0 );
	F[1] = F12(0);
	F[2] = F12(1);

	if( N == 2 )
		return std::move(F);

	for( unsigned int i = 3; i < N + 1; ++i )
	{
		F[i] = Fn( i );
	}

	return std::move( F );
}

template class TrapezoidSolver<double>;
template class SimpsonSolver<double>;
template class TrapezoidSolver<std::complex<double>>;
template class SimpsonSolver<std::complex<double>>;
