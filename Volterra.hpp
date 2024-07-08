#ifndef VOLTERRA_HPP
#define VOLTERRA_HPP

#include <functional>
#include <vector>

/*
 * Solves
 *
 * u(t) = f(t) + Int[0,t] { K(t,t') u(t') dt }
 *
 */


template<typename T> class VolterraSolver {
	public:
		using Vector = std::vector<T>;
		using SourceFunction = std::function<T(double)>;
		using KernelFunction = std::function<T(double,double)>;

		virtual Vector&& Solve( SourceFunction f, KernelFunction K, double t, unsigned int N ) = 0;
};

template<typename T> class TrapezoidSolver : public VolterraSolver<T> {
	private:
		typename VolterraSolver<T>::Vector F;
		typename VolterraSolver<T>::SourceFunction g;
		typename VolterraSolver<T>::KernelFunction K;
		double h;
		T Fn( unsigned int n );
	public:
		TrapezoidSolver() {};
		typename VolterraSolver<T>::Vector&& Solve( typename VolterraSolver<T>::SourceFunction, typename VolterraSolver<T>::KernelFunction, double, unsigned int );
};

template<typename T> class SimpsonSolver : public VolterraSolver<T> {
	private:
		typename VolterraSolver<T>::Vector F;
		typename VolterraSolver<T>::SourceFunction g;
		typename VolterraSolver<T>::KernelFunction K;
		unsigned int N;
		double h;
		T Fn( unsigned int n );
	public:
		SimpsonSolver() {};
		typename VolterraSolver<T>::Vector&& Solve( typename VolterraSolver<T>::SourceFunction, typename VolterraSolver<T>::KernelFunction, double, unsigned int );
};

#endif
