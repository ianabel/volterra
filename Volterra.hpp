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

using Vector = std::vector<double>;
using SourceFunction = std::function<double(double)>;
using KernelFunction = std::function<double(double,double)>;

class VolterraSolver {
	public:

	virtual Vector&& Solve( SourceFunction f, KernelFunction K, double t, unsigned int N ) = 0;
};

class TrapezoidSolver : public VolterraSolver {
	private:
		Vector F;
		SourceFunction g;
		KernelFunction K;
		double h;
		double Fn( unsigned int n );
	public:
		TrapezoidSolver();
		Vector&& Solve( SourceFunction, KernelFunction, double, unsigned int );
};

#endif
