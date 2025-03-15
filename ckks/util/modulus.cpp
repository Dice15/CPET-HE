
#include "modulus.h"
#include "numeric.h"
#include "arithmod.h"
#include <stdexcept>
#include <unordered_map>
#include <iostream>


namespace cpet
{
	// ModulusChain Constant
	Basis::ConstData::ConstData() :
		basis_p_(std::vector<uint64_t>()),
		basis_q_(std::vector<uint64_t>()),
		inv_q_l_hats_q_(std::vector<std::vector<uint64_t>>()),
		q_l_hats_p_(std::vector<std::vector<std::vector<uint64_t>>>()) {}


	Basis::ConstData::ConstData(
		uint64_t poly_modulus_degree,
		const std::vector<uint64_t>& basis_p_bit_sizes,
		const std::vector<uint64_t>& basis_q_bit_sizes)
	{
		uint64_t basis_p_size = basis_p_bit_sizes.size();
		uint64_t basis_q_size = basis_q_bit_sizes.size();


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

		basis_p_.reserve(basis_p_size);
		basis_q_.reserve(basis_q_size);

		for (auto basis_q_bit_size : basis_q_bit_sizes)
		{
			basis_q_.push_back(modulus_map[basis_q_bit_size].back());
			modulus_map[basis_q_bit_size].pop_back();
		}

		for (auto basis_p_bit_size : basis_p_bit_sizes)
		{
			basis_p_.push_back(modulus_map[basis_p_bit_size].back());
			modulus_map[basis_p_bit_size].pop_back();
		}

		// Compute the constants for p
		inv_P_q_.resize(basis_q_size);
		inv_p_hats_p_.resize(basis_p_size);
		p_hats_q_.resize(basis_p_size, std::vector<uint64_t>(basis_q_size));

		for (uint64_t i = 0; i < basis_p_size; i++)
		{
			uint64_t p_hat_p_i = 1;

			for (uint64_t ii = 0; ii < basis_p_size; ii++)
			{
				if (i == ii)
				{
					continue;
				}

				p_hat_p_i = mul_mod(p_hat_p_i, basis_p_[ii], basis_p_[i]);
			}

			inv_p_hats_p_[i] = inverse_mod(p_hat_p_i, basis_p_[i]);
		}

		for (uint64_t j = 0; j < basis_q_size; j++)
		{
			uint64_t P_q_j = 1;

			for (uint64_t i = 0; i < basis_p_size; i++)
			{
				P_q_j = mul_mod(P_q_j, basis_p_[i], basis_q_[j]);

				uint64_t p_hat_q_j = 1;

				for (uint64_t ii = 0; ii < basis_p_size; ii++)
				{
					if (i == ii)
					{
						continue;
					}

					p_hat_q_j = mul_mod(p_hat_q_j, basis_p_[ii], basis_q_[j]);
				}

				p_hats_q_[i][j] = p_hat_q_j;
			}

			inv_P_q_[j] = inverse_mod(P_q_j, basis_q_[j]);
		}


		// Compute the constants for q
		inv_q_l_hats_q_.resize(basis_q_size);
		q_l_hats_p_.resize(basis_q_size);

		for (uint64_t l = 0; l < basis_q_size; l++)
		{
			inv_q_l_hats_q_[l].resize(l + 1);

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

				inv_q_l_hats_q_[l][j] = inverse_mod(q_l_hat_q_j, basis_q_[j]);
			}
		}

		for (uint64_t l = 0; l < basis_q_size; l++)
		{
			q_l_hats_p_[l].resize(l + 1, std::vector<uint64_t>(basis_p_size));

			for (uint64_t i = 0; i < basis_p_size; i++)
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

	uint64_t Basis::ConstData::basis_p_size() const
	{
		return basis_p_.size();
	}

	uint64_t Basis::ConstData::basis_q_size() const
	{
		return basis_q_.size();
	}

	std::vector<uint64_t>::const_iterator Basis::ConstData::basis_p() const
	{
		return basis_p_.begin();
	}

	std::vector<uint64_t>::const_iterator Basis::ConstData::basis_q() const
	{
		return basis_q_.begin();
	}

	std::vector<uint64_t>::const_iterator Basis::ConstData::inv_P_q() const
	{
		return inv_P_q_.begin();
	}

	std::vector<uint64_t>::const_iterator Basis::ConstData::inv_p_hats_p() const
	{
		return inv_p_hats_p_.begin();
	}

	std::vector<std::vector<uint64_t>>::const_iterator Basis::ConstData::p_hats_q() const
	{
		return p_hats_q_.begin();
	}

	std::vector<std::vector<uint64_t>>::const_iterator Basis::ConstData::inv_q_l_hats_q() const
	{
		return inv_q_l_hats_q_.begin();
	}

	std::vector<std::vector<std::vector<uint64_t>>>::const_iterator Basis::ConstData::q_l_hats_p() const
	{
		return q_l_hats_p_.begin();
	}


	// Basis iterator
	Basis::basis_iterator::basis_iterator(
		std::shared_ptr<const ConstData> const_data,
		basis_type basis,
		uint64_t pop_modulus_of_basis_q_count,
		iterator_type iter)
	{
		const_data_ = const_data;
		begin_index_ = (basis == basis_type::basis_d ? 0 : const_data->basis_p_size());
		end_index_ = const_data->basis_p_size() + const_data->basis_q_size() - pop_modulus_of_basis_q_count;
		curr_index_ = (iter == iterator_type::begin ? begin_index_ : end_index_);
	}

	uint64_t Basis::basis_iterator::operator-(basis_iterator other) const
	{
		return end_index_ - begin_index_;
	}

	Basis::basis_iterator Basis::basis_iterator::operator+(uint64_t offset) const
	{
		if (this->curr_index_ + offset > this->end_index_)
		{
			throw std::out_of_range("Iterator addition out of range");
		}

		basis_iterator result = *this;

		result.curr_index_ += offset;

		return result;
	}

	Basis::basis_iterator& Basis::basis_iterator::operator++()
	{
		if (curr_index_ >= end_index_)
		{
			throw std::out_of_range("Iterator cannot be incremented beyond the valid range");
		}

		curr_index_++;
		return *this;
	}

	Basis::basis_iterator Basis::basis_iterator::operator++(int) 
	{ 
		basis_iterator temp = *this;
		++(*this);
		return temp;
	}

	Basis::basis_iterator& Basis::basis_iterator::operator--()
	{ 
		if (curr_index_ <= begin_index_)
		{
			throw std::out_of_range("Iterator cannot be decremented below the valid range");
		}

		curr_index_--;
		return *this; 
	}

	Basis::basis_iterator Basis::basis_iterator::operator--(int)
	{
		basis_iterator temp = *this;
		--(*this);
		return temp;
	}

	bool Basis::basis_iterator::operator==(const basis_iterator& other) const
	{
		return const_data_ == other.const_data_
			&& begin_index_ == other.begin_index_
			&& end_index_ == other.end_index_
			&& curr_index_ == other.curr_index_;
	}

	bool Basis::basis_iterator::operator!=(const basis_iterator& other) const
	{
		return !(const_data_ == other.const_data_
			&& begin_index_ == other.begin_index_
			&& end_index_ == other.end_index_
			&& curr_index_ == other.curr_index_);
	}

	const uint64_t& Basis::basis_iterator::operator*() const
	{
		if (curr_index_ < begin_index_ || curr_index_ >= end_index_)
		{
			throw std::out_of_range("Dereferencing iterator out of valid range");
		}

		return curr_index_ < const_data_->basis_p_size()
			? const_data_->basis_p()[curr_index_]
			: const_data_->basis_q()[curr_index_ - const_data_->basis_p_size()];
	}


	// Modulus chain
	Basis::Basis() :
		const_data_(std::make_shared<const Basis::ConstData>()),
		pop_modulus_of_basis_q_count_(0),
		curr_basis_(basis_type::basis_q) {}

	Basis::Basis(std::shared_ptr<const ConstData> const_data) :
		const_data_(const_data),
		pop_modulus_of_basis_q_count_(0),
		curr_basis_(basis_type::basis_q) {}

	Basis::Basis(const Basis& basis) :
		const_data_(basis.const_data_),
		pop_modulus_of_basis_q_count_(basis.pop_modulus_of_basis_q_count_),
		curr_basis_(basis.curr_basis_) {}

	bool Basis::operator==(const Basis& other) const
	{
		return const_data_ == other.const_data_
			&& pop_modulus_of_basis_q_count_ == other.pop_modulus_of_basis_q_count_
			&& curr_basis_ == other.curr_basis_;
	}

	bool Basis::operator!=(const Basis& other) const
	{
		return !(const_data_ == other.const_data_
			&& pop_modulus_of_basis_q_count_ == other.pop_modulus_of_basis_q_count_
			&& curr_basis_ == other.curr_basis_);
	}

	uint64_t Basis::size() const
	{
		return end() - begin();
	}

	Basis::basis_iterator Basis::begin() const
	{
		return Basis::basis_iterator(const_data_, curr_basis_, pop_modulus_of_basis_q_count_, iterator_type::begin);
	}

	Basis::basis_iterator Basis::end() const
	{
		return Basis::basis_iterator(const_data_, curr_basis_, pop_modulus_of_basis_q_count_, iterator_type::end);
	}

	void Basis::convert_basis(Basis::basis_type basis)
	{
		curr_basis_ = basis;
	}

	Basis::basis_type Basis::current_basis() const
	{
		return curr_basis_;
	}

	void Basis::pop_modulus_of_basis_q()
	{
		if (const_data_->basis_q_size() < 2)
		{
			throw std::out_of_range("Basis q must have at least one elements.");
		}

		pop_modulus_of_basis_q_count_++;
	}

	uint64_t Basis::basis_p_size() const
	{
		return const_data_->basis_p_size();
	}

	uint64_t Basis::basis_q_size() const
	{
		return const_data_->basis_q_size() - pop_modulus_of_basis_q_count_;
	}

	std::vector<uint64_t>::const_iterator Basis::basis_p() const
	{
		return const_data_->basis_p();
	}

	std::vector<uint64_t>::const_iterator Basis::basis_q() const
	{
		return const_data_->basis_q();
	}

	std::vector<uint64_t>::const_iterator Basis::inv_P_q() const
	{
		return const_data_->inv_P_q();
	}

	std::vector<uint64_t>::const_iterator Basis::inv_p_hats_p() const
	{
		return const_data_->inv_p_hats_p();
	}

	std::vector<std::vector<uint64_t>>::const_iterator Basis::p_hats_q() const
	{
		return const_data_->p_hats_q();
	}

	std::vector<uint64_t>::const_iterator Basis::inv_q_hats_q() const
	{
		return const_data_->inv_q_l_hats_q()[basis_q_size() - 1].begin();
	}

	std::vector<std::vector<uint64_t>>::const_iterator Basis::q_hats_p() const
	{
		return const_data_->q_l_hats_p()[basis_q_size() - 1].begin();
	}
}