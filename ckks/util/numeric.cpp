
#include "numeric.h"
#include "arithmod.h"


namespace cpet
{
    bool is_over_60_bit(uint64_t value)
    {
        return (value & 0xF000000000000000ULL) != 0;
    }

    bool is_power_of_two(uint64_t value)
    {
        return value > 0 && (value & (value - 1)) == 0;
    }

    bool is_prime(uint64_t value, size_t num_rounds)
    {
        if (value < 2)
        {
            return false;
        }
        if (2 == value)
        {
            return true;
        }
        if (0 == (value & 0x1))
        {
            return false;
        }
        if (3 == value)
        {
            return true;
        }
        if (0 == (value % 3))
        {
            return false;
        }
        if (5 == value)
        {
            return true;
        }
        if (0 == (value % 5))
        {
            return false;
        }
        if (7 == value)
        {
            return true;
        }
        if (0 == (value % 7))
        {
            return false;
        }
        if (11 == value)
        {
            return true;
        }
        if (0 == (value % 11))
        {
            return false;
        }
        if (13 == value)
        {
            return true;
        }
        if (0 == (value % 13))
        {
            return false;
        }

        // Miller-Rabin test.
        // Find r and odd d that satisfy value = 2^r * d + 1.
        uint64_t d = value - 1;
        uint64_t r = 0;
        while (0 == (d & 0x1))
        {
            d >>= 1;
            r++;
        }
        if (r == 0)
        {
            return false;
        }

        // 1) Pick a = 2, check a^(value - 1).
        // 2) Pick a randomly from [3, value - 1], check a^(value - 1).
        // 3) Repeat 2) for another num_rounds - 2 times.
        std::random_device rand;
        std::uniform_int_distribution<uint64_t> dist(3, value - 1);
        uint64_t modulus = value;

        for (size_t i = 0; i < num_rounds; i++)
        {
            uint64_t a = i ? dist(rand) : 2;
            uint64_t x = pow_mod(a, d, modulus);
            if (x == 1 || x == value - 1)
            {
                continue;
            }
            uint64_t count = 0;
            do
            {
                x = mul_mod(x, x, modulus);
                count++;
            } while (x != value - 1 && count < r - 1);
            if (x != value - 1)
            {
                return false;
            }
        }
        return true;
    }

    void try_prime(uint64_t factor, uint64_t bit_size, uint64_t& destination)
    {
        if (bit_size > 60 || bit_size < 2)
        {
            throw std::invalid_argument("The prime number must be greater than or equal to 2 and less than or equal to 60 bits.");
        }

        // (2^bit_size - 1) / factor * factor + 1 부터 factor를 빼면서 소수를 찾음.
        uint64_t candidate = ((static_cast<uint64_t>(1) << bit_size) - 1) / factor * factor + 1;
        uint64_t lower_bound = static_cast<uint64_t>(1) << (bit_size - 1);
        bool found = false;

        while (!found && candidate > lower_bound)
        {
            if (is_prime(candidate))
            {
                destination = candidate;
                found = true;
            }

            candidate -= factor;
        }

        if (found == false)
        {
            throw std::logic_error("failed to find enough qualifying prime");
        }
    }

    void try_primes(const uint64_t factor, const uint64_t bit_size, uint64_t count, std::vector<uint64_t>& destination)
    {
        if (bit_size > 60 || bit_size < 2)
        {
            throw std::invalid_argument("The prime number must be greater than or equal to 2 and less than or equal to 60 bits.");
        }

        // (2^bit_size - 1) / factor * factor + 1 부터 factor를 빼면서 소수를 찾음.
        uint64_t candidate = ((static_cast<uint64_t>(1) << bit_size) - 1) / factor * factor + 1;
        uint64_t lower_bound = static_cast<uint64_t>(1) << (bit_size - 1);
      
        destination.clear();

        while (count >0 && candidate > lower_bound)
        {
            if (is_prime(candidate))
            {
                destination.push_back(candidate);
                count--;
            }

            candidate -= factor;
        }

        if (count > 0)
        {
            throw std::logic_error("failed to find enough qualifying primes");
        }
    }
}