
#include "arithmod.h"


namespace cpet
{
    // TODO: 빠른 모듈러 연산이 가능한 barrett reduction 구현.
    uint64_t mod(uint64_t value, uint64_t modulus)
    {
        return value < modulus ? value : value % modulus;
    }

    uint64_t add_mod(uint64_t value1, uint64_t value2, uint64_t modulus)
    {
        return mod(value1 + value2, modulus);
    }

    uint64_t sub_mod(uint64_t value1, uint64_t value2, uint64_t modulus)
    {
        return mod(value1 + negate_mod(value2, modulus), modulus);
    }

    uint64_t mul_mod(uint64_t value1, uint64_t value2, uint64_t modulus)
    {
        if (value1 == 0ULL || value2 == 0ULL)
        {
            return 0ULL;
        }

        uint64_t high = 0ULL;
        uint64_t low = _umul128(value1, value2, &high);

        if (high == 0ULL)
        {
            return mod(low, modulus);
        }

        uint64_t remainder = 0ULL;
        _udiv128(high, low, modulus, &remainder);

        return remainder;
    }

    uint64_t pow_mod(uint64_t base, uint64_t exponent, uint64_t modulus)
    {
        if (exponent == 0ULL)
        {
            return 1ULL;
        }

        if (exponent == 1ULL)
        {
            return mod(base, modulus);
        }

        base = mod(base, modulus);

        uint64_t result = 1ULL;

        // Binary exponentiation: 3^13 = 3^(1101)_2 = 3^8 * 3^4 * 3^1
        while (exponent > 0ULL)
        {
            if (exponent & 1ULL)
            {
                result = mul_mod(result, base, modulus);
            }

            base = mul_mod(base, base, modulus);

            exponent >>= 1ULL;
        }

        return result;
    }

    uint64_t negate_mod(uint64_t value, uint64_t modulus)
    {
        if (value == 0ULL)
        {
            return 0ULL;
        }

        uint64_t result = mod(value, modulus);

        return result == 0ULL ? 0ULL : modulus - result;
    }

    // The modulus must be prime.
    // However, this is not a problem since all moduli in the RNS CKKS scheme are prime.
    uint64_t inverse_mod(uint64_t value, uint64_t modulus)
    {
        uint64_t result = 1ULL;
        uint64_t exp = modulus - 2ULL;

        // a^-1 ≡ a^(p - 2) mod p, where p is prime.
        while (exp > 0)
        {
            if (exp & 1)
            {
                result = mul_mod(result, value, modulus);
            }

            value = mul_mod(value, value, modulus);
            exp >>= 1;
        }

        return result;
    }
}