
#include "basis.h"
#include "numeric.h"
#include "arithmod.h"
#include <stdexcept>
#include <unordered_map>


namespace cpet
{
	Basis::Basis() :
		constant_(nullptr),
		basis_type_(BasisType::basis_q),
		drop_q_count_(0)
	{}

	Basis::Basis(
		uint64_t poly_modulus_degree,
		const std::vector<uint64_t>& moduli_p_bit_sizes,
		const std::vector<uint64_t>& moduli_q_bit_sizes
	) :
		basis_type_(BasisType::basis_q),
		drop_q_count_(0)
	{
		// Check count of basis.
		if (moduli_p_bit_sizes.empty())
		{
			throw std::invalid_argument("The basis p cannot be empty.");
		}

		if (moduli_q_bit_sizes.empty())
		{
			throw std::invalid_argument("The basis q cannot be empty.");
		}


		// Create constant
		Constant constant;

		constant.p_begin_ = 0;
		constant.p_end_ = moduli_p_bit_sizes.size();

		constant.q_begin_ = constant.p_end_;
		constant.q_end_ = constant.q_begin_ + moduli_q_bit_sizes.size();


		// Create basis p and q
		std::unordered_map<uint64_t, uint64_t> modulus_bit_sizes_count;
		std::unordered_map<uint64_t, std::vector<uint64_t>> modulus_map;

		for (uint64_t modulus_bit_size : moduli_p_bit_sizes)
		{
			if (is_over_60_bit(modulus_bit_size))
			{
				throw std::invalid_argument("The modulus must be a prime number greater than or equal to 2 and less than or equal to 60 bits.");
			}

			modulus_bit_sizes_count[modulus_bit_size]++;
		}

		for (uint64_t modulus_bit_size : moduli_q_bit_sizes)
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

		std::vector<uint64_t> moduli_p;
		std::vector<uint64_t> moduli_q;

		moduli_p.reserve(moduli_p_bit_sizes.size());
		moduli_q.reserve(moduli_q_bit_sizes.size());
		constant.moduli_.reserve(moduli_p_bit_sizes.size() + moduli_q_bit_sizes.size());

		for (auto basis_p_bit_size : moduli_p_bit_sizes)
		{
			moduli_p.push_back(modulus_map[basis_p_bit_size].back());		
			modulus_map[basis_p_bit_size].pop_back();
			constant.moduli_.push_back(moduli_p.back());
		}

		for (auto basis_q_bit_size : moduli_q_bit_sizes)
		{
			moduli_q.push_back(modulus_map[basis_q_bit_size].back());
			modulus_map[basis_q_bit_size].pop_back();
			constant.moduli_.push_back(moduli_q.back());
		}


		// Compute the constants for fast basis convert (p <-> q)
		constant.P_mod_ql_.resize(moduli_q_bit_sizes.size());
		constant.inv_P_mod_ql_.resize(moduli_q_bit_sizes.size());
		constant.inv_p_hats_mod_p.resize(moduli_p_bit_sizes.size());
		constant.p_hats_mod_ql_.resize(moduli_p_bit_sizes.size(), std::vector<uint64_t>(moduli_q_bit_sizes.size()));

		for (uint64_t i = 0; i < moduli_p_bit_sizes.size(); i++)
		{
			uint64_t p_hat_p_i = 1;

			for (uint64_t ii = 0; ii < moduli_p_bit_sizes.size(); ii++)
			{
				if (i == ii)
				{
					continue;
				}

				p_hat_p_i = mul_mod(p_hat_p_i, moduli_p[ii], moduli_p[i]);
			}

			constant.inv_p_hats_mod_p[i] = inverse_mod(p_hat_p_i, moduli_p[i]);
		}

		for (uint64_t j = 0; j < moduli_q_bit_sizes.size(); j++)
		{
			uint64_t P_q_j = 1;

			for (uint64_t i = 0; i < moduli_p_bit_sizes.size(); i++)
			{
				P_q_j = mul_mod(P_q_j, moduli_p[i], moduli_q[j]);

				uint64_t p_hat_q_j = 1;

				for (uint64_t ii = 0; ii < moduli_p_bit_sizes.size(); ii++)
				{
					if (i == ii)
					{
						continue;
					}

					p_hat_q_j = mul_mod(p_hat_q_j, moduli_p[ii], moduli_q[j]);
				}

				constant.p_hats_mod_ql_[i][j] = p_hat_q_j;
			}
	
			constant.P_mod_ql_[j] = P_q_j;
			constant.inv_P_mod_ql_[j] = inverse_mod(P_q_j, moduli_q[j]);
		}


		// Compute the constants fast basis convert (p <-> q)
		constant.inv_ql_hats_mod_ql_.resize(moduli_q_bit_sizes.size());
		constant.ql_hats_mod_p.resize(moduli_q_bit_sizes.size());

		for (uint64_t l = 0; l < moduli_q_bit_sizes.size(); l++)
		{
			constant.inv_ql_hats_mod_ql_[l].resize(l + 1);

			for (uint64_t j = 0; j <= l; j++)
			{
				uint64_t q_l_hat_q_j = 1;

				for (uint64_t jj = 0; jj <= l; jj++)
				{
					if (j == jj)
					{
						continue;
					}

					q_l_hat_q_j = mul_mod(q_l_hat_q_j, moduli_q[jj], moduli_q[j]);
				}

				constant.inv_ql_hats_mod_ql_[l][j] = inverse_mod(q_l_hat_q_j, moduli_q[j]);
			}
		}

		for (uint64_t l = 0; l < moduli_q_bit_sizes.size(); l++)
		{
			constant.ql_hats_mod_p[l].resize(l + 1, std::vector<uint64_t>(moduli_p_bit_sizes.size()));

			for (uint64_t i = 0; i < moduli_p_bit_sizes.size(); i++)
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

						q_l_hat_p_i = mul_mod(q_l_hat_p_i, moduli_q[jj], moduli_p[i]);
					}

					constant.ql_hats_mod_p[l][j][i] = q_l_hat_p_i;
				}
			}
		}


		// init constant
		constant_ = std::make_shared<const Constant>(constant);
	}

	bool Basis::operator==(const Basis& other) const
	{
		return constant_ == other.constant_
			&& basis_type_ == other.basis_type_
			&& drop_q_count_ == other.drop_q_count_;
	}

	bool Basis::operator!=(const Basis& other) const
	{
		return constant_ != other.constant_
			|| basis_type_ != other.basis_type_
			|| drop_q_count_ != other.drop_q_count_;
	}

	Basis::BasisType Basis::get_type() const
	{
		return basis_type_;
	}

	void Basis::set_type(BasisType basis_type)
	{
		basis_type_ = basis_type;
	}

	const std::vector<uint64_t>& Basis::get_moduli() const
	{
		return constant_->moduli_;
	}

	uint64_t Basis::begin() const
	{
		switch (basis_type_)
		{
		case BasisType::basis_q: 
		{
			return constant_->q_begin_;
		}
		case BasisType::basis_pq:
		{
			return constant_->p_begin_;
		}
		default:
			throw std::out_of_range("Invaild basis type.");
		}
	}

	uint64_t Basis::end() const
	{
		switch (basis_type_)
		{
		case BasisType::basis_q:
		{
			return constant_->q_end_ - drop_q_count_;
		}
		case BasisType::basis_pq:
		{
			return constant_->q_end_ - drop_q_count_;
		}
		default:
			throw std::out_of_range("Invaild basis type.");
		}
	}

	uint64_t Basis::size() const
	{
		return end() - begin();
	}

	uint64_t Basis::capacity() const
	{
		return constant_->q_end_;
	}

	uint64_t Basis::at(uint64_t index) const
	{
		if (index < begin() || end() <= index)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->moduli_[index];
	}

	void Basis::drop_q()
	{
		if (static_cast<int64_t>(constant_->q_end_) - ++drop_q_count_ == constant_->q_begin_)
		{
			throw std::out_of_range("Basis q must have at least one elements.");
		}
	}

	const std::vector<uint64_t>& Basis::P_mod_ql() const
	{
		return constant_->P_mod_ql_;
	}

	const std::vector<uint64_t>& Basis::inv_P_mod_ql() const
	{
		return constant_->inv_P_mod_ql_;
	}

	const std::vector<uint64_t>& Basis::inv_p_hats_mod_p() const
	{
		return constant_->inv_p_hats_mod_p;
	}

	const std::vector<std::vector<uint64_t>>& Basis::p_hats_mod_ql() const
	{
		return constant_->p_hats_mod_ql_;
	}

	const std::vector<uint64_t>& Basis::inv_ql_hats_mod_ql() const
	{
		return constant_->inv_ql_hats_mod_ql_[constant_->q_end_ - constant_->q_begin_ - 1ULL];
	}

	const std::vector<std::vector<uint64_t>>& Basis::ql_hats_mod_p() const
	{
		return constant_->ql_hats_mod_p[constant_->q_end_ - constant_->q_begin_ - 1ULL];
	}
}