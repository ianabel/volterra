
#include <cmath>

double Gamma0( double L1, double L2 )
{
	return std::exp( -(L1 + L2)/2.0 ) * std::cyl_bessel_i( 0.0, std::sqrt( L1 * L2 ) );
}

double Lambda( double L1, double L2 )
{

}
