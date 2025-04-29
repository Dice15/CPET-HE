
#include "rnscycloring.h"
#include "arithmod.h"
#include <stdexcept>
#include <iostream>

namespace cpet
{
	RnsCycloRing::RnsCycloRing() :
		scale_(1.0),
		poly_modulus_degree_(0),
		basis_(Basis()),
		rns_coeffs_({}),
		ntt_handler_(nullptr),
		form_(RnsCycloRing::Form::coeff)
	{}

	RnsCycloRing::RnsCycloRing(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		RnsCycloRing::Form default_form,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		scale_(scale),
		poly_modulus_degree_(poly_modulus_degree),
		basis_(basis),
		rns_coeffs_(std::vector<std::vector<uint64_t>>(basis.capacity())),
		form_(default_form),
		ntt_handler_(ntt_handler)
	{
		for (uint64_t basis_idx = basis.begin(); basis_idx < basis.end(); ++basis_idx)
		{
			rns_coeffs_[basis_idx].resize(poly_modulus_degree);
		}
	}

	void RnsCycloRing::operator()(const RnsCycloRing& other)
	{
		if (scale_ != other.scale_ || poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
			{
				rns_coeffs_[basis_idx][coeff_idx] = other.rns_coeffs_[basis_idx][coeff_idx];
			}
		}
	}

	uint64_t RnsCycloRing::get(uint64_t coeff_idx) const
	{
		if (coeff_idx >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		// TODO: CRT를 한 결과를 리턴하도록 수정해야함.
		return rns_coeffs_[basis_.begin()][coeff_idx];
	}

	void RnsCycloRing::set(uint64_t coeff_idx, uint64_t value)
	{
		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			if (coeff_idx >= poly_modulus_degree_)
			{
				throw std::out_of_range("Index out of range.");
			}

			rns_coeffs_[basis_idx][coeff_idx] = mod(value, basis_.at(basis_idx));
		}
	}

	uint64_t RnsCycloRing::get(uint64_t basis_idx, uint64_t coeff_idx) const
	{
		if (basis_idx < basis_.begin() || basis_idx >= basis_.end())
		{
			throw std::out_of_range("Index out of range.");
		}

		if (coeff_idx >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		return rns_coeffs_[basis_idx][coeff_idx];
	}

	void RnsCycloRing::set(uint64_t basis_idx, uint64_t coeff_idx, uint64_t value)
	{
		if (basis_idx < basis_.begin() || basis_idx >= basis_.end())
		{
			throw std::out_of_range("Index out of range.");
		}

		if (coeff_idx >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		rns_coeffs_[basis_idx][coeff_idx] = mod(value, basis_.at(basis_idx));
	}

	double_t RnsCycloRing::scale() const
	{
		return scale_;
	}

	uint64_t RnsCycloRing::poly_modulus_degree() const
	{
		return poly_modulus_degree_;
	}

	const Basis& RnsCycloRing::basis() const
	{
		return basis_;
	}

	RnsCycloRing::Form RnsCycloRing::form() const
	{
		return form_;
	}

	void RnsCycloRing::modulus_reduction()
	{
		if (basis_.level() < 2)
		{
			throw std::out_of_range("레벨 2미만에서는 modulus reduction을 진행할 수 없습니다.");
		}

		basis_.drop_basis();
		rns_coeffs_[basis_.end()] = {};
	}

	void RnsCycloRing::rescale()
	{
		if (basis_.level() < 2)
		{
			throw std::out_of_range("레벨 2미만에서는 modulus reduction을 진행할 수 없습니다.");
		}

		slot_to_coeff();

		uint64_t last_basis_idx = basis_.end() - 1ULL;

		for (uint64_t basis_idx = basis_.begin(); basis_idx < last_basis_idx; ++basis_idx)
		{
			uint64_t inv_modulus = inverse_mod(basis_.at(last_basis_idx), basis_.at(basis_idx));

			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
			{
				uint64_t residue = sub_mod(
					rns_coeffs_[basis_idx][coeff_idx], rns_coeffs_[last_basis_idx][coeff_idx], basis_.at(basis_idx)
				);

				rns_coeffs_[basis_idx][coeff_idx] = mul_mod(
					inv_modulus, residue, basis_.at(basis_idx)
				);
			}
		}

		scale_ /= basis_.at(basis_.end() - 1);
		basis_.drop_basis();
		rns_coeffs_[basis_.end()] = {};

		coeff_to_slot();
	}

	void RnsCycloRing::convert_scale_force(double_t scale)
	{
		scale_ = scale;
	}

	void RnsCycloRing::convert_basis_force(Basis::Type type)
	{
		if (basis_.get_basis_type() == type)
		{
			throw std::out_of_range("이미 같은 기저입니다.");
		}

		basis_.convert_basis(type);

		for (uint64_t basis_idx = 0; basis_idx < basis_.begin(); ++basis_idx)
		{
			rns_coeffs_[basis_idx] = {};
		}

		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			rns_coeffs_[basis_idx].resize(poly_modulus_degree_);
		}

		for (uint64_t basis_idx = basis_.end(); basis_idx < rns_coeffs_.size(); ++basis_idx)
		{
			rns_coeffs_[basis_idx] = {};
		}
	}

	void RnsCycloRing::convert_basis_approximate(Basis::Type type)
	{
		if (basis_.get_basis_type() == type)
		{
			throw std::out_of_range("이미 같은 기저입니다.");
		}

		slot_to_coeff();

		switch (type)
		{
		case Basis::Type::Q:
		{	
			basis_reduction_approximate();		
			break;
		}
		case Basis::Type::PQ:
		{
			basis_raising_approximate();
			break;
		}
		default:
			throw std::out_of_range("Invaild basis type.");
		}

		coeff_to_slot();
	}

	void RnsCycloRing::coeff_to_slot()
	{
		if (form_ == RnsCycloRing::Form::slot)
		{
			throw std::out_of_range("This ring is already slot form.");
		}

		ntt_handler_->ntt_negacyclic(rns_coeffs_, basis_.begin(), basis_.end());

		form_ = RnsCycloRing::Form::slot;
	}

	void RnsCycloRing::slot_to_coeff()
	{
		if (form_ == RnsCycloRing::Form::coeff)
		{
			throw std::out_of_range("This ring is already coeff form.");
		}

		ntt_handler_->inverse_ntt_negacyclic(rns_coeffs_, basis_.begin(), basis_.end());

		form_ = RnsCycloRing::Form::coeff;
	}

	void RnsCycloRing::add_inplace(const RnsCycloRing& other)
	{
		if (scale_ != other.scale_ || poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (form_ != RnsCycloRing::Form::slot || other.form_ != RnsCycloRing::Form::slot)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in slot form.");
		}

		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[basis_idx][i] = add_mod(rns_coeffs_[basis_idx][i], other.rns_coeffs_[basis_idx][i], basis_.at(basis_idx));
			}
		}
	}

	void RnsCycloRing::sub_inplace(const RnsCycloRing& other)
	{
		if (scale_ != other.scale_ || poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (form_ != RnsCycloRing::Form::slot || other.form_ != RnsCycloRing::Form::slot)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in slot form.");
		}

		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[basis_idx][i] = sub_mod(rns_coeffs_[basis_idx][i], other.rns_coeffs_[basis_idx][i], basis_.at(basis_idx));
			}
		}
	}

	void RnsCycloRing::mul_inplace(const RnsCycloRing& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (form_ != RnsCycloRing::Form::slot || other.form_ != RnsCycloRing::Form::slot)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in slot form.");
		}

		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[basis_idx][i] = mul_mod(rns_coeffs_[basis_idx][i], other.rns_coeffs_[basis_idx][i], basis_.at(basis_idx));
			}
		}

		scale_ *= other.scale_;
	}

	void RnsCycloRing::mul_inplace(uint64_t value)
	{
		if (form_ != RnsCycloRing::Form::slot)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in slot form.");
		}

		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
			{
				rns_coeffs_[basis_idx][i] = mul_mod(rns_coeffs_[basis_idx][i], value, basis_.at(basis_idx));
			}
		}
	}

	void RnsCycloRing::negate_inplace()
	{
		for (uint64_t basis_idx = basis_.begin(); basis_idx < basis_.end(); ++basis_idx)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
			{
				rns_coeffs_[basis_idx][coeff_idx] = negate_mod(rns_coeffs_[basis_idx][coeff_idx], basis_.at(basis_idx));
			}
		}
	}

	void RnsCycloRing::basis_raising_approximate()
	{
		if (basis_.get_basis_type() == Basis::Type::PQ)
		{
			throw std::out_of_range("이미 상위 기저입니다.");
		}

		uint64_t basis_p_begin = 0;
		uint64_t basis_p_end = basis_.begin();
		uint64_t basis_q_begin = basis_.begin();
		uint64_t basis_q_end = basis_.end();

		basis_.convert_basis(Basis::Type::PQ);

		for (uint64_t basis_p_idx = basis_p_begin; basis_p_idx < basis_p_end; ++basis_p_idx)
		{
			rns_coeffs_[basis_p_idx].resize(poly_modulus_degree_);
		}

		for (uint64_t basis_p_idx = basis_p_begin; basis_p_idx < basis_p_end; ++basis_p_idx)
		{
			for (uint64_t basis_q_idx = basis_q_begin; basis_q_idx < basis_q_end; ++basis_q_idx)
			{
				for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
				{
					if (basis_q_idx == basis_q_begin)
					{
						rns_coeffs_[basis_p_idx][coeff_idx] = 0;
					}

					uint64_t residue = mul_mod(
						rns_coeffs_[basis_q_idx][coeff_idx], basis_.inv_q_hats_mod_q(basis_q_idx), basis_.at(basis_q_idx)
					);

					residue = mod(
						residue, basis_.at(basis_p_idx)
					);

					residue = mul_mod(
						residue, basis_.q_hats_mod_p(basis_q_idx, basis_p_idx), basis_.at(basis_p_idx)
					);

					rns_coeffs_[basis_p_idx][coeff_idx] = add_mod(
						rns_coeffs_[basis_p_idx][coeff_idx], residue, basis_.at(basis_p_idx)
					);
				}
			}
		}
	}

	void RnsCycloRing::basis_reduction_approximate()
	{
		if (basis_.get_basis_type() == Basis::Type::Q)
		{
			throw std::out_of_range("이미 하위 기저입니다.");
		}

		basis_.convert_basis(Basis::Type::Q);

		uint64_t basis_p_begin = 0;
		uint64_t basis_p_end = basis_.begin();
		uint64_t basis_q_begin = basis_.begin();
		uint64_t basis_q_end = basis_.end();

		for (uint64_t basis_q_idx = basis_q_begin; basis_q_idx < basis_q_end; ++basis_q_idx)
		{
			for (uint64_t basis_p_idx = basis_p_begin; basis_p_idx < basis_p_begin; ++basis_p_idx)
			{
				for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
				{
					uint64_t residue = mul_mod(
						rns_coeffs_[basis_p_idx][coeff_idx], basis_.inv_p_hats_mod_p(basis_p_idx), basis_.at(basis_p_idx)
					);

					residue = mod(
						residue, basis_.at(basis_q_idx)
					);

					residue = mul_mod(
						residue, basis_.p_hats_mod_q(basis_p_idx, basis_q_idx), basis_.at(basis_q_idx)
					);

					rns_coeffs_[basis_q_idx][coeff_idx] = sub_mod(
						rns_coeffs_[basis_q_idx][coeff_idx], residue, basis_.at(basis_q_idx)
					);
				}
			}
		}

		for (uint64_t basis_p_idx = basis_p_begin; basis_p_idx < basis_p_end; ++basis_p_idx)
		{
			rns_coeffs_[basis_p_idx] = {};
		}

		for (uint64_t basis_q_idx = basis_q_begin; basis_q_idx < basis_q_end; ++basis_q_idx)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
			{
				rns_coeffs_[basis_q_idx][coeff_idx] = mul_mod(
					rns_coeffs_[basis_q_idx][coeff_idx], basis_.inv_P_mod_q(basis_q_idx), basis_.at(basis_q_idx)
				);
			}
		}
	}
}