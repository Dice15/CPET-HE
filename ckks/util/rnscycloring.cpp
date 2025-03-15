
#include "rnscycloring.h"
#include "arithmod.h"
#include "ntt.h"
#include <stdexcept>


namespace cpet
{
	RnsCycloRing::RnsCycloRing() : RnsCycloRing(0, Basis()) {}

	RnsCycloRing::RnsCycloRing(const PolyModulus& poly_modulus, const Basis& basis, uint64_t value)
	{
		poly_modulus_ = poly_modulus;
		basis_ = basis;
		congruences_.assign(basis_.size(), std::vector<uint64_t>(poly_modulus_.degree()));

		auto congruence = congruences_.begin();

		for (auto modulus = basis_.begin(); modulus != basis_.end(); ++modulus, ++congruence)
		{
			for (auto& coeff : *congruence)
			{
				coeff = mod(value, *modulus);
			}
		}

		ntt_form_ = false;
	}

	RnsCycloRing::RnsCycloRing(RnsCycloRing& other)
	{
		poly_modulus_ = other.poly_modulus_;
		basis_ = other.basis_;
		congruences_ = other.congruences_;
		ntt_form_ = other.ntt_form_;
	}

	uint64_t RnsCycloRing::operator()(uint64_t congruence_index, uint64_t coeff_index) const
	{
		if (congruence_index >= congruences_.size() || coeff_index >= congruences_[congruence_index].size())
		{
			throw std::out_of_range("Index out of range.");
		}

		return congruences_[congruence_index][coeff_index];
	}

	void RnsCycloRing::operator()(uint64_t congruence_index, uint64_t coeff_index, uint64_t value)
	{
		if (congruence_index >= congruences_.size() || coeff_index >= congruences_[congruence_index].size())
		{
			throw std::out_of_range("Index out of range.");
		}

		congruences_[congruence_index][coeff_index] = mod(value, *(basis_.begin() + congruence_index));
	}

	void RnsCycloRing::assign(const PolyModulus& poly_modulus, const Basis& basis, uint64_t value)
	{
		poly_modulus_ = poly_modulus;
		basis_ = basis;
		congruences_.assign(basis_.size(), std::vector<uint64_t>(poly_modulus_.degree()));

		auto congruence = congruences_.begin();

		for (auto modulus = basis_.begin(); modulus != basis_.end(); ++modulus, ++congruence)
		{
			for (auto& coeff : *congruence)
			{
				coeff = mod(value, *modulus);
			}
		}

		ntt_form_ = false;
	}

	const PolyModulus& RnsCycloRing::poly_modulus() const
	{
		return poly_modulus_;
	}

	uint64_t RnsCycloRing::basis_size() const
	{
		return basis_.size();
	}

	Basis::basis_iterator RnsCycloRing::basis_begin() const
	{
		return basis_.begin();
	}

	Basis::basis_iterator RnsCycloRing::basis_end() const
	{
		return basis_.end();
	}

	void RnsCycloRing::fast_basis_conversion(
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
		if (basis_.current_basis() == basis_type)
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
				basis_.basis_p(),
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
				basis_.basis_p(),
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
	}

	void RnsCycloRing::set_ntt_form()
	{
		if (ntt_form_)
		{
			throw std::out_of_range("This ring is already ntt form.");
		}

		ntt_negacyclic(*this);

		ntt_form_ = true;
	}

	void RnsCycloRing::set_normal_form()
	{
		if (!ntt_form_)
		{
			throw std::out_of_range("This ring is already normal form.");
		}

		inverse_ntt_negacyclic(*this);

		ntt_form_ = false;
	}

	void RnsCycloRing::add(const RnsCycloRing& other, RnsCycloRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || basis_ != other.basis_)
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.basis_ = basis_;
		destination.congruences_.resize(congruences_.size(), std::vector<uint64_t>(poly_modulus_.degree()));
		destination.ntt_form_ = ntt_form_;

		auto modulus = basis_.begin();

		for (uint64_t m = 0; m < basis_.size(); ++m, ++modulus)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); ++i)
			{
				destination.congruences_[m][i] = add_mod(congruences_[m][i], other.congruences_[m][i], *modulus);
			}
		}
	}

	void RnsCycloRing::sub(const RnsCycloRing& other, RnsCycloRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || basis_ != other.basis_)
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.basis_ = basis_;
		destination.congruences_.resize(congruences_.size(), std::vector<uint64_t>(poly_modulus_.degree()));
		destination.ntt_form_ = ntt_form_;

		auto modulus = basis_.begin();

		for (uint64_t m = 0; m < basis_.size(); ++m, ++modulus)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); ++i)
			{
				destination.congruences_[m][i] = sub_mod(congruences_[m][i], other.congruences_[m][i], *modulus);
			}
		}
	}

	void RnsCycloRing::mul(const RnsCycloRing& other, RnsCycloRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || basis_ != other.basis_)
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.basis_ = basis_;
		destination.congruences_.resize(congruences_.size(), std::vector<uint64_t>(poly_modulus_.degree()));
		destination.ntt_form_ = ntt_form_;

		auto modulus = basis_.begin();

		for (uint64_t m = 0; m < basis_.size(); ++m, ++modulus)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); ++i)
			{
				destination.congruences_[m][i] = mul_mod(congruences_[m][i], other.congruences_[m][i], *modulus);
			}
		}
	}

	void RnsCycloRing::negate(RnsCycloRing& destination) const
	{
		destination.poly_modulus_ = poly_modulus_;
		destination.basis_ = basis_;
		destination.congruences_.resize(congruences_.size());
		destination.ntt_form_ = ntt_form_;

		auto modulus = basis_.begin();

		for (uint64_t m = 0; m < basis_.size(); ++m, ++modulus)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); ++i)
			{
				destination.congruences_[m][i] = negate_mod(congruences_[m][i], *modulus);
			}
		}
	}
}