
#include "rnscycloring.h"
#include "arithmod.h"
#include <stdexcept>
#include <iostream>

namespace cpet
{
	RnsCycloRing::RnsCycloRing() :
		poly_modulus_degree_(0),
		basis_(nullptr),
		rns_coeffs_({}),
		ntt_handler_(nullptr),
		ntt_form_(false)
	{}

	RnsCycloRing::RnsCycloRing(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		poly_modulus_degree_(poly_modulus_degree),
		basis_(&basis),
		rns_coeffs_(std::vector<std::vector<uint64_t>>(
			basis.capacity(),
			std::vector<uint64_t>(poly_modulus_degree_, 0))
		),
		ntt_handler_(ntt_handler),
		ntt_form_(false)
	{}

	RnsCycloRing::RnsCycloRing(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		int64_t value
	) :
		poly_modulus_degree_(poly_modulus_degree),
		basis_(&basis),
		rns_coeffs_(std::vector<std::vector<uint64_t>>(
			basis.capacity(),
			std::vector<uint64_t>(poly_modulus_degree_, 0))
		),
		ntt_handler_(ntt_handler),
		ntt_form_(false)
	{
		for (uint64_t rns_idx = basis.begin(); rns_idx < basis.end(); ++rns_idx)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree; ++coeff_idx)
			{
				if (coeff_idx >= rns_coeffs_[rns_idx].size())
				{
					throw std::out_of_range("Index out of range.");
				}

				if (value < 0)
				{
					rns_coeffs_[rns_idx][coeff_idx] = negate_mod(std::llabs(value), basis_->at(rns_idx));
				}
				else
				{
					rns_coeffs_[rns_idx][coeff_idx] = mod(value, basis_->at(rns_idx));
				}
			}
		}
	}

	void RnsCycloRing::operator()(const RnsCycloRing& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_)
		{
			throw std::out_of_range("Poly modulus degree is mismatched.");
		}

		if (basis_ != other.basis_ || ntt_handler_ != other.ntt_handler_)
		{
			throw std::out_of_range("Basis is mismatched.");
		}

		for (uint64_t rns_idx = basis_->begin(); rns_idx < basis_->end(); ++rns_idx)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
			{
				rns_coeffs_[rns_idx][coeff_idx] = other.rns_coeffs_[rns_idx][coeff_idx];
			}
		}

		ntt_form_ = other.ntt_form_;
	}

	int64_t RnsCycloRing::get_coeff(uint64_t coeff_idx) const
	{
		uint64_t rns_idx = basis_->begin();

		if (coeff_idx >= rns_coeffs_[rns_idx].size())
		{
			throw std::out_of_range("Index out of range.");
		}
		
		if (rns_coeffs_[rns_idx][coeff_idx] < (basis_->at(rns_idx) >> 1))
		{
			return rns_coeffs_[rns_idx][coeff_idx];
		}
		else
		{
			return -static_cast<int64_t>(negate_mod(rns_coeffs_[rns_idx][coeff_idx], basis_->at(rns_idx)));
		}
	}

	void RnsCycloRing::set_coeff(uint64_t coeff_idx, int64_t value)
	{
		for (uint64_t rns_idx = basis_->begin(); rns_idx < basis_->end(); ++rns_idx)
		{
			if (coeff_idx >= rns_coeffs_[rns_idx].size())
			{
				throw std::out_of_range("Index out of range.");
			}

			if (value < 0)
			{
				rns_coeffs_[rns_idx][coeff_idx] = negate_mod(std::llabs(value), basis_->at(rns_idx));
			}
			else
			{
				rns_coeffs_[rns_idx][coeff_idx] = mod(value, basis_->at(rns_idx));
			}
		}
	}

	uint64_t RnsCycloRing::get_rns_coeff(uint64_t rns_idx, uint64_t coeff_idx) const
	{
		if (rns_idx >= rns_coeffs_.size())
		{
			throw std::out_of_range("Index out of range.");
		}

		if (coeff_idx >= rns_coeffs_[rns_idx].size())
		{
			throw std::out_of_range("Index out of range.");
		}

		return rns_coeffs_[rns_idx][coeff_idx];
	}

	void RnsCycloRing::set_rns_coeff(uint64_t rns_idx, uint64_t coeff_idx, uint64_t value)
	{
		if (rns_idx >= rns_coeffs_.size())
		{
			throw std::out_of_range("Index out of range.");
		}

		if (coeff_idx >= rns_coeffs_[rns_idx].size())
		{
			throw std::out_of_range("Index out of range.");
		}

		rns_coeffs_[rns_idx][coeff_idx] = mod(value, basis_->at(rns_idx));
	}

	uint64_t RnsCycloRing::poly_modulus_degree() const
	{
		return poly_modulus_degree_;
	}

	const Basis& RnsCycloRing::get_basis() const
	{
		return *basis_;
	}

	void RnsCycloRing::set_basis(const Basis* const basis)
	{
		basis_ = basis;

		for (uint64_t b = 0; b < basis_->begin(); ++b)
		{
			rns_coeffs_[b] = {};
		}

		for (uint64_t b = basis_->end(); b < rns_coeffs_.size(); ++b)
		{
			rns_coeffs_[b] = {};
		}
	}

	/*void RnsCycloRing::fast_basis_conversion(
		uint64_t basis_a_size,
		uint64_t basis_b_size,
		std::vector<uint64_t>::const_iterator basis_a,
		std::vector<uint64_t>::const_iterator basis_b,
		std::vector<uint64_t>::const_iterator inv_a_hats_a,
		std::vector<std::vector<uint64_t>>::const_iterator a_hats_b,
		uint64_t congruence_size,
		const std::vector<std::vector<uint64_t>>& congruences_a,
		std::vector<std::vector<uint64_t>>& destination) const
	{
		std::vector<std::vector<uint64_t>> result(basis_b_size, std::vector<uint64_t>(congruence_size, 0));

		for (uint64_t i = 0; i < basis_b_size; i++)
		{
			for (uint64_t j = 0; j < basis_a_size; j++)
			{
				for (uint64_t c = 0; c < congruence_size; c++)
				{
					uint64_t inner = mul_mod(congruences_a[j][c], inv_a_hats_a[j], basis_a[j]);
					uint64_t outer = mul_mod(inner, a_hats_b[j][i], basis_b[i]);
					result[i][c] = add_mod(result[i][c], outer, basis_b[i]);
				}
			}
		}

		for (uint64_t i = 0; i < basis_b_size; i++)
		{
			for (uint64_t j = 0; j < basis_a_size; j++)
			{
				uint64_t inv_a_hat_a = 1;
				uint64_t a_hat_b = 1;

				for (uint64_t jj = 0; jj < basis_a_size; jj++)
				{
					if (j == jj) continue;
					inv_a_hat_a = mul_mod(inv_a_hat_a, basis_a[jj], basis_a[j]);
					a_hat_b = mul_mod(a_hat_b, basis_a[jj], basis_b[i]);
				}

				inv_a_hat_a = inverse_mod(inv_a_hat_a, basis_a[j]);

				if (inv_a_hats_a[j] != inv_a_hat_a)
				{
					throw std::invalid_argument("");
				}

				if (a_hats_b[j][i] != a_hat_b)
				{
					throw std::invalid_argument("");
				}
			}
		}

		destination = result;
	}

	void RnsCycloRing::approximate_modulus_raising(
		uint64_t basis_p_size,
		uint64_t basis_q_size,
		std::vector<uint64_t>::const_iterator basis_p,
		std::vector<uint64_t>::const_iterator basis_q,
		std::vector<uint64_t>::const_iterator inv_q_hats_q,
		std::vector<std::vector<uint64_t>>::const_iterator q_hats_p,
		uint64_t congruence_size,
		const std::vector<std::vector<uint64_t>>& congruences_q,
		std::vector<std::vector<uint64_t>>& destination) const
	{
		std::vector<std::vector<uint64_t>> congruences_d;

		fast_basis_conversion(
			basis_q_size, basis_p_size, basis_q, basis_p, inv_q_hats_q, q_hats_p, congruence_size, congruences_q, congruences_d);

		congruences_d.resize(basis_p_size + basis_q_size, std::vector<uint64_t>(congruence_size, 0));

		destination = congruences_d;
	}

	void RnsCycloRing::approximate_modulus_reduction(
		uint64_t basis_p_size,
		uint64_t basis_q_size,
		std::vector<uint64_t>::const_iterator basis_p,
		std::vector<uint64_t>::const_iterator basis_q,
		std::vector<uint64_t>::const_iterator inv_P_q,
		std::vector<uint64_t>::const_iterator inv_p_hats_p,
		std::vector<std::vector<uint64_t>>::const_iterator p_hats_q,
		uint64_t congruence_size,
		const std::vector<std::vector<uint64_t>>& congruences_d,
		std::vector<std::vector<uint64_t>>& destination) const
	{
		std::vector<std::vector<uint64_t>> congruences_q;

		fast_basis_conversion(
			basis_p_size, basis_q_size, basis_p, basis_q, inv_p_hats_p, p_hats_q, congruence_size, congruences_d, congruences_q);


		for (uint64_t j = 0; j < basis_q_size; j++)
		{
			for (uint64_t c = 0; c < congruence_size; c++)
			{
				uint64_t diff = sub_mod(congruences_d[basis_p_size + j][c], congruences_q[j][c], basis_q[j]);
				congruences_q[j][c] = mul_mod(inv_P_q[j], diff, basis_q[j]);
			}
		}

		destination = congruences_q;
	}

	void RnsCycloRing::convert_basis(Basis::basis_type basis_type)
	{
		if (basis_->bas == basis_type)
		{
			throw std::invalid_argument("Already same basis.");
		}

		switch (basis_type)
		{
		case Basis::basis_type::basis_d:
		{
			approximate_modulus_raising(
				basis_.basis_p_size(),
				basis_.basis_q_size(),
				basis_.basis_p_begin(),
				basis_.basis_q(),
				basis_.inv_q_hats_q(),
				basis_.q_hats_p(),
				poly_modulus_.degree(),
				congruences_,
				congruences_);
			break;
		}
		case cpet::Basis::basis_type::basis_q:
		{
			approximate_modulus_reduction(
				basis_.basis_p_size(),
				basis_.basis_q_size(),
				basis_.basis_p_begin(),
				basis_.basis_q(),
				basis_.inv_P_q(),
				basis_.inv_p_hats_p(),
				basis_.p_hats_q(),
				poly_modulus_.degree(),
				congruences_,
				congruences_);
			break;
		}
		default:
			throw std::invalid_argument("Unsupported basis.");
			break;
		}

		basis_.convert_basis(basis_type);
	}*/

	void RnsCycloRing::set_ntt_form()
	{
		if (ntt_form_)
		{
			throw std::out_of_range("This ring is already ntt form.");
		}

		ntt_handler_->ntt_negacyclic(rns_coeffs_, basis_->begin(), basis_->end());
		ntt_form_ = true;
	}

	void RnsCycloRing::set_normal_form()
	{
		if (!ntt_form_)
		{
			throw std::out_of_range("This ring is already normal form.");
		}

		ntt_handler_->inverse_ntt_negacyclic(rns_coeffs_, basis_->begin(), basis_->end());
		ntt_form_ = false;
	}

	void RnsCycloRing::add_inplace(const RnsCycloRing& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_)
		{
			throw std::invalid_argument("Rings must have matching poly modulus degree and basis.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		for (uint64_t b = basis_->begin(); b < basis_->end(); ++b)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[b][i] = add_mod(rns_coeffs_[b][i], other.rns_coeffs_[b][i], basis_->at(b));
			}
		}
	}

	void RnsCycloRing::sub_inplace(const RnsCycloRing& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_)
		{
			throw std::invalid_argument("Rings must have matching poly modulus degree and basis.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		for (uint64_t b = basis_->begin(); b < basis_->end(); ++b)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[b][i] = sub_mod(rns_coeffs_[b][i], other.rns_coeffs_[b][i], basis_->at(b));
			}
		}
	}

	void RnsCycloRing::mul_inplace(const RnsCycloRing& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_)
		{
			throw std::invalid_argument("Rings must have matching degree and basis.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		for (uint64_t b = basis_->begin(); b < basis_->end(); ++b)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[b][i] = mul_mod(rns_coeffs_[b][i], other.rns_coeffs_[b][i], basis_->at(b));
			}
		}
	}

	void RnsCycloRing::negate_inplace()
	{
		for (uint64_t b = basis_->begin(); b < basis_->end(); ++b)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[b][i] = negate_mod(rns_coeffs_[b][i], basis_->at(b));
			}
		}
	}

/*	void RnsCycloRing::modulus_reduction()
	{
		basis_->level_down();
		congruences_[basis_.end()].clear();
	}

	void RnsCycloRing::rescale()
	{
		set_normal_form();

		uint64_t last_index = basis_.end() - 1ULL;
		uint64_t last_modulus = basis_[last_index];

		for (uint64_t b = basis_.begin(); b < last_index; ++b)
		{
			uint64_t inv_modulus = inverse_mod(last_modulus, basis_[b]);

			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				congruences_[b][i] = mul_mod(inv_modulus, sub_mod(congruences_[b][i], congruences_[last_index][i], basis_[b]), basis_[b]);
			}
		}

		scale_ /= static_cast<double_t>(last_modulus);

		modulus_reduction();

		set_ntt_form();
	}*/
}