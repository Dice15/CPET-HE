#pragma once

#include <cstdint>
#include <random>


namespace cpet
{
	bool is_over_60_bit(uint64_t value);

	bool is_power_of_two(uint64_t poly_modulus_degree);

	bool is_prime(uint64_t value, size_t num_rounds = 40);

	void try_prime(const uint64_t factor, const uint64_t bit_size, uint64_t& destination);
}