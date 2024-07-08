
#include <complex>
#include <cmath>

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

		double Gamma0( double L1, double L2 ) const {
			double x = (L1 + L2)/2.0;
			double y = (L1 - L2)/2.0;
			return std::exp( -x ) * std::cyl_bessel_i( 0.0, std::sqrt( x * x - y * y ) ); 
		};

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
			return std::exp( -(dt2*dt2) ) * (
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
