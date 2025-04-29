
#include "basis.h"
#include "numeric.h"
#include "arithmod.h"
#include <stdexcept>
#include <unordered_map>


namespace cpet
{
	Basis::Basis() :
		constant_(nullptr),
		basis_type_(Basis::Type::Q),
		drop_q_count_(0)
	{}

	Basis::Basis(
		uint64_t poly_modulus_degree,
		Basis::Type default_type,
		const std::vector<uint64_t>& moduli_p_bit_sizes,
		const std::vector<uint64_t>& moduli_q_bit_sizes
	) :
		basis_type_(default_type),
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
		constant.P_mod_p_and_q_.resize(moduli_q_bit_sizes.size());
		constant.inv_P_mod_q_.resize(moduli_q_bit_sizes.size());
		constant.inv_p_hats_mod_p.resize(moduli_p_bit_sizes.size());
		constant.p_hats_mod_q_.resize(moduli_p_bit_sizes.size(), std::vector<uint64_t>(moduli_q_bit_sizes.size()));

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

				constant.p_hats_mod_q_[i][j] = p_hat_q_j;
			}
	
			constant.P_mod_p_and_q_[j] = P_q_j;
			constant.inv_P_mod_q_[j] = inverse_mod(P_q_j, moduli_q[j]);
		}


		// Compute the constants fast basis convert (p <-> q)
		constant.inv_q_hats_mod_q_by_lev_.resize(moduli_q_bit_sizes.size());
		constant.q_hats_mod_p_by_lev.resize(moduli_q_bit_sizes.size());

		for (uint64_t l = 0; l < moduli_q_bit_sizes.size(); l++)
		{
			constant.inv_q_hats_mod_q_by_lev_[l].resize(l + 1);

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

				constant.inv_q_hats_mod_q_by_lev_[l][j] = inverse_mod(q_l_hat_q_j, moduli_q[j]);
			}
		}

		for (uint64_t l = 0; l < moduli_q_bit_sizes.size(); l++)
		{
			constant.q_hats_mod_p_by_lev[l].resize(l + 1, std::vector<uint64_t>(moduli_p_bit_sizes.size()));

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

					constant.q_hats_mod_p_by_lev[l][j][i] = q_l_hat_p_i;
				}
			}
		}


		// init constant
		constant_ = std::make_shared<const Constant>(constant);
	}

	bool Basis::operator==(const Basis& other) const
	{
		if (constant_ != other.constant_ || basis_type_ != other.basis_type_)
		{
			throw std::out_of_range("Constant and type of basis is mismatched.");
		}

		return drop_q_count_ == other.drop_q_count_;
	}

	bool Basis::operator!=(const Basis& other) const
	{
		if (constant_ != other.constant_ || basis_type_ != other.basis_type_)
		{
			throw std::out_of_range("Constant and type of basis is mismatched.");
		}

		return drop_q_count_ != other.drop_q_count_;
	}

	bool Basis::operator<(const Basis& other) const
	{
		if (constant_ != other.constant_ || basis_type_ != other.basis_type_)
		{
			throw std::out_of_range("Constant and type of basis is mismatched.");
		}

		return drop_q_count_ > other.drop_q_count_;
	}

	bool Basis::operator>(const Basis& other) const
	{
		if (constant_ != other.constant_ || basis_type_ != other.basis_type_)
		{
			throw std::out_of_range("Constant and type of basis is mismatched.");
		}

		return drop_q_count_ < other.drop_q_count_;
	}

	const std::vector<uint64_t>& Basis::get_moduli() const
	{
		return constant_->moduli_;
	}

	uint64_t Basis::begin() const
	{
		switch (basis_type_)
		{
		case Basis::Type::Q:
		{
			return constant_->q_begin_;
		}
		case Basis::Type::PQ:
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
		case Basis::Type::Q:
		{
			return constant_->q_end_ - drop_q_count_;
		}
		case Basis::Type::PQ:
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

	uint64_t Basis::level() const
	{
		return constant_->q_end_ - constant_->q_begin_ - drop_q_count_;
	}

	void Basis::drop_basis()
	{
		if (static_cast<int64_t>(constant_->q_end_) - ++drop_q_count_ == constant_->q_begin_)
		{
			throw std::out_of_range("Basis q must have at least one elements.");
		}
	}

	Basis::Type Basis::get_basis_type() const
	{
		return basis_type_;
	}

	void Basis::convert_basis(Basis::Type basis_type)
	{
		basis_type_ = basis_type;
	}

	uint64_t Basis::capacity() const
	{
		return constant_->q_end_;
	}

	uint64_t Basis::at(uint64_t basis_idx) const
	{
		if (constant_->q_end_ <= basis_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->moduli_[basis_idx];
	}

	uint64_t Basis::P_mod_p_and_q(uint64_t basis_idx) const
	{
		if (constant_->q_end_ <= basis_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return basis_idx < constant_->q_begin_ ? 0 : constant_->P_mod_p_and_q_[basis_idx - constant_->q_begin_];
	}

	uint64_t Basis::inv_P_mod_q(uint64_t basis_q_idx) const
	{
		if (basis_q_idx < constant_->q_begin_ || constant_->q_end_ <= basis_q_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->inv_P_mod_q_[basis_q_idx - constant_->q_begin_];
	}

	uint64_t Basis::inv_p_hats_mod_p(uint64_t basis_p_idx) const
	{
		if (basis_p_idx < constant_->p_begin_ || constant_->p_end_ <= basis_p_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->inv_p_hats_mod_p[basis_p_idx];
	}

	uint64_t Basis::p_hats_mod_q(uint64_t basis_p_idx, uint64_t basis_q_idx) const
	{
		if (basis_p_idx < constant_->p_begin_ || constant_->p_end_ <= basis_p_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		if (basis_q_idx < constant_->q_begin_ || constant_->q_end_ <= basis_q_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->p_hats_mod_q_[basis_p_idx][basis_q_idx - constant_->q_begin_];
	}

	uint64_t Basis::inv_q_hats_mod_q(uint64_t basis_q_idx) const
	{
		if (basis_q_idx < constant_->q_begin_ || constant_->q_end_ <= basis_q_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->inv_q_hats_mod_q_by_lev_[level() - 1][basis_q_idx - constant_->q_begin_];
	}

	uint64_t Basis::q_hats_mod_p(uint64_t basis_q_idx, uint64_t basis_p_idx) const
	{
		if (basis_q_idx < constant_->q_begin_ || constant_->q_end_ <= basis_q_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		if (basis_p_idx < constant_->p_begin_ || constant_->p_end_ <= basis_p_idx)
		{
			throw std::out_of_range("Index out of range.");
		}

		return constant_->q_hats_mod_p_by_lev[level() - 1][basis_q_idx - constant_->q_begin_][basis_p_idx];
	}
}