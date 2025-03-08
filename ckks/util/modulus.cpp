
#include "modulus.h"
#include "numeric.h"


namespace cpet
{
    Modulus::Modulus() :Modulus(2) {}

    Modulus::Modulus(uint64_t modulus)
    {
        if (is_over_60_bit(modulus) || !is_prime(modulus) || modulus < 2)
        {
            throw std::invalid_argument("The modulus must be a prime number greater than or equal to 2 and less than or equal to 60 bits.");
        }

        modulus_ = modulus;
    }

	uint64_t Modulus::value() const
    {
        return modulus_;
    }
}