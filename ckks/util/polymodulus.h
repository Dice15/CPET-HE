#pragma once

#include <cstdint>
#include <stdexcept>


namespace cpet
{
	class PolyModulus
	{
	public:
		PolyModulus();

		PolyModulus(uint64_t poly_modulus_degree);

		uint64_t degree() const;

	private:
		uint64_t poly_modulus_degree_;
	};
}