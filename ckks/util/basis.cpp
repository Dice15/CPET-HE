
#include "basis.h"
#include "numeric.h"
#include "arithmod.h"
#include <stdexcept>
#include <unordered_map>


namespace cpet
{
	// Basis Constants
	Basis::constant::constant(
		uint64_t poly_modulus_degree,
		const std::vector<uint64_t>& basis_p_bit_sizes,
		const std::vector<uint64_t>& basis_q_bit_sizes)
	{
		basis_p_size_ = basis_p_bit_sizes.size();
		basis_q_size_ = basis_q_bit_sizes.size();
		basis_d_size_ = basis_p_size_ + basis_q_size_;


		// Create basis p and q
		std::unordered_map<uint64_t, uint64_t> modulus_bit_sizes_count;
		std::unordered_map<uint64_t, std::vector<uint64_t>> modulus_map;

		for (auto modulus_bit_size : basis_p_bit_sizes)
		{
			if (is_over_60_bit(modulus_bit_size))
			{
				throw std::invalid_argument("The modulus must be a prime number greater than or equal to 2 and less than or equal to 60 bits.");
			}

			modulus_bit_sizes_count[modulus_bit_size]++;
		}

		for (auto modulus_bit_size : basis_q_bit_sizes)
		{
			if (is_over_60_bit(modulus_bit_size))
			{
				throw std::invalid_argument("The modulus must be a prime number greater than or equal to 2 and less than or equal to 60 bits.");
			}

			modulus_bit_sizes_count[modulus_bit_size]++;
		}

		for (const auto& [modulus_bit_size, count] : modulus_bit_sizes_count)
		{
			try_primes(poly_modulus_degree * 2, modulus_bit_size, count, modulus_map[modulus_bit_size]);
		}

		for (auto& [bit_size, moduli] : modulus_map)
		{
			std::sort(moduli.begin(), moduli.end(), std::greater<uint64_t>());
		}

		basis_p_.reserve(basis_p_size_);
		basis_q_.reserve(basis_q_size_);
		basis_d_.reserve(basis_p_size_ + basis_q_size_);

		for (auto basis_p_bit_size : basis_p_bit_sizes)
		{
			basis_p_.push_back(modulus_map[basis_p_bit_size].back());
			basis_d_.push_back(modulus_map[basis_p_bit_size].back());
			modulus_map[basis_p_bit_size].pop_back();
		}

		for (auto basis_q_bit_size : basis_q_bit_sizes)
		{
			basis_q_.push_back(modulus_map[basis_q_bit_size].back());
			basis_d_.push_back(modulus_map[basis_q_bit_size].back());
			modulus_map[basis_q_bit_size].pop_back();
		}


		// Compute the constants for basis p
		inv_P_q_l_.resize(basis_q_size_);
		inv_p_hats_p_.resize(basis_p_size_);
		p_hats_q_l_.resize(basis_p_size_, std::vector<uint64_t>(basis_q_size_));

		for (uint64_t i = 0; i < basis_p_size_; i++)
		{
			uint64_t p_hat_p_i = 1;

			for (uint64_t ii = 0; ii < basis_p_size_; ii++)
			{
				if (i == ii)
				{
					continue;
				}

				p_hat_p_i = mul_mod(p_hat_p_i, basis_p_[ii], basis_p_[i]);
			}

			inv_p_hats_p_[i] = inverse_mod(p_hat_p_i, basis_p_[i]);
		}

		for (uint64_t j = 0; j < basis_q_size_; j++)
		{
			uint64_t P_q_j = 1;

			for (uint64_t i = 0; i < basis_p_size_; i++)
			{
				P_q_j = mul_mod(P_q_j, basis_p_[i], basis_q_[j]);

				uint64_t p_hat_q_j = 1;

				for (uint64_t ii = 0; ii < basis_p_size_; ii++)
				{
					if (i == ii)
					{
						continue;
					}

					p_hat_q_j = mul_mod(p_hat_q_j, basis_p_[ii], basis_q_[j]);
				}

				p_hats_q_l_[i][j] = p_hat_q_j;
			}

			inv_P_q_l_[j] = inverse_mod(P_q_j, basis_q_[j]);
		}


		// Compute the constants for basis q
		inv_q_l_hats_q_l_.resize(basis_q_size_);
		q_l_hats_p_.resize(basis_q_size_);

		for (uint64_t l = 0; l < basis_q_size_; l++)
		{
			inv_q_l_hats_q_l_[l].resize(l + 1);

			for (uint64_t j = 0; j <= l; j++)
			{
				uint64_t q_l_hat_q_j = 1;

				for (uint64_t jj = 0; jj <= l; jj++)
				{
					if (j == jj)
					{
						continue;
					}

					q_l_hat_q_j = mul_mod(q_l_hat_q_j, basis_q_[jj], basis_q_[j]);
				}

				inv_q_l_hats_q_l_[l][j] = inverse_mod(q_l_hat_q_j, basis_q_[j]);
			}
		}

		for (uint64_t l = 0; l < basis_q_size_; l++)
		{
			q_l_hats_p_[l].resize(l + 1, std::vector<uint64_t>(basis_p_size_));

			for (uint64_t i = 0; i < basis_p_size_; i++)
			{
				for (uint64_t j = 0; j <= l; j++)
				{
					uint64_t q_l_hat_p_i = 1;

					for (uint64_t jj = 0; jj <= l; jj++)
					{
						if (j == jj)
						{
							continue;
						}

						q_l_hat_p_i = mul_mod(q_l_hat_p_i, basis_q_[jj], basis_p_[i]);
					}

					q_l_hats_p_[l][j][i] = q_l_hat_p_i;
				}
			}
		}
	}


	// Basis
	Basis::Basis(std::shared_ptr<const constant> const_data) :
		const_data_(const_data),
		pop_count_of_basis_q_(0),
		curr_basis_(basis_type::basis_q) {}

	bool Basis::operator==(const Basis& other) const
	{
		return const_data_ == other.const_data_
			&& pop_count_of_basis_q_ == other.pop_count_of_basis_q_
			&& curr_basis_ == other.curr_basis_;
	}

	bool Basis::operator!=(const Basis& other) const
	{
		return const_data_ != other.const_data_
			|| pop_count_of_basis_q_ != other.pop_count_of_basis_q_
			|| curr_basis_ != other.curr_basis_;
	}

	uint64_t Basis::operator[](uint64_t index) const
	{
		if (index < begin() || end() <= index)
		{
			throw std::out_of_range("Index out of range.");
		}

		return const_data_->basis_d_[index];
	}

	uint64_t Basis::begin() const
	{
		return curr_basis_ == Basis::basis_type::basis_d ? 0 : const_data_->basis_p_size_;
	}

	uint64_t Basis::end() const
	{
		return const_data_->basis_d_size_ - pop_count_of_basis_q_;
	}

	uint64_t Basis::size() const
	{
		return end() - begin();
	}

	void Basis::convert_basis(Basis::basis_type basis)
	{
		curr_basis_ = basis;
	}

	void Basis::pop_basis_q()
	{
		if (const_data_->basis_q_.size() < 2)
		{
			throw std::out_of_range("Basis q must have at least one elements.");
		}

		++pop_count_of_basis_q_;
	}

	uint64_t Basis::basis_p_size() const
	{
		return const_data_->basis_p_.size();
	}

	uint64_t Basis::basis_q_size() const
	{
		return const_data_->basis_q_.size() - pop_count_of_basis_q_;
	}

	uint64_t Basis::basis_d_size() const
	{
		return const_data_->basis_d_.size() - pop_count_of_basis_q_;
	}

	const std::vector<uint64_t>& Basis::basis_p() const
	{
		return const_data_->basis_p_;
	}

	const std::vector<uint64_t>& Basis::basis_q() const
	{
		return const_data_->basis_q_;
	}

	const std::vector<uint64_t>& Basis::basis_d() const
	{
		return const_data_->basis_d_;
	}

	const std::vector<uint64_t>& Basis::inv_P_q() const
	{
		return const_data_->inv_P_q_l_;
	}

	const std::vector<uint64_t>& Basis::inv_p_hats_p() const
	{
		return const_data_->inv_p_hats_p_;
	}

	const std::vector<std::vector<uint64_t>>& Basis::p_hats_q() const
	{
		return const_data_->p_hats_q_l_;
	}

	const std::vector<uint64_t>& Basis::inv_q_hats_q() const
	{
		return const_data_->inv_q_l_hats_q_l_[basis_q_size() - 1ULL];
	}

	const std::vector<std::vector<uint64_t>>& Basis::q_hats_p() const
	{
		return const_data_->q_l_hats_p_[basis_q_size() - 1ULL];
	}
}