
#include <cmath>

class LinearGK
{
	private:
		double TauOverZ,qOverEpsilon;
		double OmegaS,OmegaStar;
		double eta;
		double ky;
	public:
		double Gamma0( double L1, double L2 ) const { return std::exp( -(L1 + L2)/2.0 ) * std::cyl_bessel_i( 0.0, std::sqrt( L1 * L2 ) ); }; 
		double Lambda( double L1, double L2 ) const { 
			double x = std::sqrt( L1 * L2 );
			return 1.0 - (L1 + L2)/2.0 + x * std::cyl_bessel_i( 1.0, x ) / std::cyl_bessel_i( 0.0, x );
		}
		double lambda( double t ) const { return (ky*ky + OmegaS * OmegaS * t * t ) / 2.0; };

		double K( double t, double tPrime ) const {
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

		double Kernel( double t, double tPrime ) const {
			return K(t,tPrime)/g(t);
		};

		double Source( double t ) const { 
			return g(0)/g(t);
		};

}

