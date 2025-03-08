#pragma once

#include <cstdint>
#include <stdexcept>


namespace cpet
{
	class Modulus
	{
	public:
		Modulus();

		Modulus(uint64_t modulus);

		uint64_t value() const;

	private:
		uint64_t modulus_;
	};
}