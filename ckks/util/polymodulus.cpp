
#include "polymodulus.h"
#include "numeric.h"


namespace cpet
{
	PolyModulus::PolyModulus() :PolyModulus(2) {}

	PolyModulus::PolyModulus(uint64_t poly_modulus_degree)
	{
		if (!is_power_of_two(poly_modulus_degree) || poly_modulus_degree < 2)
		{
			throw std::invalid_argument("The polynomial modulus degree must be a power of two and greater than or equal to 2.");
		}

		poly_modulus_degree_ = poly_modulus_degree;
	}

	uint64_t PolyModulus::degree() const
	{
		return poly_modulus_degree_;
	}
}