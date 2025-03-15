#pragma once

#include <cstdint>
#include <intrin.h>


namespace cpet
{
    uint64_t mod(uint64_t value, uint64_t modulus);

    uint64_t mod_switch(uint64_t value, uint64_t modulus1, uint64_t modulus2);

    uint64_t add_mod(uint64_t value1, uint64_t value2, uint64_t modulus);

    uint64_t sub_mod(uint64_t value1, uint64_t value2, uint64_t modulus);

    uint64_t mul_mod(uint64_t value1, uint64_t value2, uint64_t modulus);

    uint64_t pow_mod(uint64_t base, uint64_t exponent, uint64_t modulus);

    uint64_t negate_mod(uint64_t value, uint64_t modulus);

    uint64_t inverse_mod(uint64_t value, uint64_t modulus);
}