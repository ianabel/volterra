
#include <complex>
#include <cmath>
#include <limits>
#include <boost/math/special_functions/factorials.hpp>

using std::exp;
using std::sqrt;
using std::log;
using std::cyl_bessel_i;

class LinearGK
{
	private:
		double TauOverZ,qOverEpsilon;
		double OmegaS,OmegaStar;
		double eta;
		double ky;
	public:

		LinearGK( double tau, double qeps, double om_s, double om_star, double eta_i, double ky_rho )
			: TauOverZ( tau ), qOverEpsilon( qeps ), OmegaS( om_s ), OmegaStar( om_star ), eta( eta_i ), ky( ky_rho )
		{
		}

		static double ScaledBesselI( double x ) {
			if ( x > 22.0 ) { /* Asymptotic series expansion */
				double ExpI0 = 1.0;	// Factor out the 1/sqrt(2*pi*x)
				double Ak = 1.0;
				double z = 1.0/x;
				double log_z = log(z);
				const unsigned int k_max = 46;
				for( unsigned int k = 0; k < k_max; ++k )
				{
					// A_{k+1}
					double Ak1 = Ak * z * (0.5 + k) * (0.5 + k) * (0.5) / ( 1.0 + k );
					ExpI0 += Ak1;
					Ak = Ak1;
				}
				return ExpI0/sqrt( 2.0 * M_PI * x );
			} else {
				return exp(-x) * cyl_bessel_i(0,x);
			}
		}

		double Gamma0( double L1, double L2 ) const {
			double x = (L1 + L2)/2.0;
			double s = sqrt( L1 * L2 ); 

			// Gamma0 = exp(-x)*I_0(s) = exp( s - x ) * [ exp(-s) * I_0(s) ]
			double E = exp( s - x );
			double ScaledI = ScaledBesselI( s );

			return E * ScaledI;	

		};

		double Lambda( double L1, double L2 ) const {
			double x = std::sqrt( L1 * L2 );
			return 1.0 - (L1 + L2)/2.0 + x * std::cyl_bessel_i( 1.0, x ) / std::cyl_bessel_i( 0.0, x );
		}

		double lambda( double t ) const { 
			return (ky*ky + OmegaS * OmegaS * t * t ) / 2.0;
		};

		std::complex<double> K( double t, double tPrime ) const {
			double dt2 = (t - tPrime)/2.0;
			double l1 = lambda(t);
			double l2 = lambda(tPrime);
			std::complex<double> I(0.0,1.0);
			return std::exp( -1.0 * (dt2*dt2) ) * (
					( qOverEpsilon * OmegaS - 1.0 ) * dt2
					- I * OmegaStar * ( 1.0 + eta * ( Lambda( l1, l2 ) - 1 - dt2*dt2 ) )
					) * Gamma0( l1, l2 );
		};

		double g( double t ) const {
			double l = lambda(t);
			return 1.0 + TauOverZ - Gamma0( l, l );
		};

		std::complex<double> Kernel( double t, double tPrime ) const {
			return K(t,tPrime)/g(t);
		};

		std::complex<double> Source( double t ) const {
			return g(0)/g(t);
		};

		std::function< std::complex<double> (double) > getSourceFn() {
			return std::bind(&LinearGK::Source,this,std::placeholders::_1);
		}

		std::function< std::complex<double> (double,double) > getKernelFn() {
			return std::bind(&LinearGK::Kernel,this,std::placeholders::_1,std::placeholders::_2);
		}
};

/*
class SlabLinearGK
{
	private:
		double TauOverZ,qOverEpsilon;
		double OmegaS,OmegaStar;
		double eta;
		double ky;
	public:

		SlabLinearGK( double tau, double qeps, double om_s, double om_star, double eta_i, double ky_rho )
			: TauOverZ( tau ), qOverEpsilon( qeps ), OmegaS( om_s ), OmegaStar( om_star ), eta( eta_i ), ky( ky_rho )
		{
		}

		double Gamma0( double L1, double L2 ) const { return std::exp( -(L1 + L2)/2.0 ) * std::cyl_bessel_i( 0.0, std::sqrt( L1 * L2 ) ); }; 
		double Lambda( double L1, double L2 ) const { 
			double x = std::sqrt( L1 * L2 );
			return 1.0 - (L1 + L2)/2.0 + x * std::cyl_bessel_i( 1.0, x ) / std::cyl_bessel_i( 0.0, x );
		}
		double lambda( double t ) const { return (ky*ky + OmegaS * OmegaS * t * t ) / 2.0; };

		std::complex<double> K( double t, double tPrime ) const {
			double dt2 = (t - tPrime)/2.0;
			double l1 = lambda(t);
			double l2 = lambda(tPrime);
			std::complex<double> I(0.0,1.0);
			return std::exp( -(dt2*dt2) ) * ( ( qOverEpsilon * OmegaS - 1.0 ) * dt2 - I * OmegaStar * ( 1.0 + eta * ( Lambda( l1, l2 ) - 1 - dt2*dt2 ) ) ) * Gamma0( l1, l2 );
		};

		double g( double t ) const {
			double l = lambda(t);
			return 1.0 + TauOverZ - Gamma0( l, l );
		};

		std::complex<double> Kernel( double t, double tPrime ) const {
			return K(t,tPrime)/g(t);
		};

		std::complex<double> Source( double t ) const {
			return g(0)/g(t);
		};

		std::function< std::complex<double> (double) > getSourceFn() {
			return std::bind(&SlabLinearGK::Source,this,std::placeholders::_1);
		}

		std::function< std::complex<double> (double,double) > getKernelFn() {
			return std::bind(&SlabLinearGK::Kernel,this,std::placeholders::_1,std::placeholders::_2);
		}
};
*/
