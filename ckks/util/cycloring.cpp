
#include "cycloring.h"
#include <stdexcept>


namespace cpet
{
	CycloRing::CycloRing() :
		poly_modulus_degree_(0),
		coeffs_({}),
		fft_handler_(nullptr),
		form_(CycloRing::Form::slot)
	{}

	CycloRing::CycloRing(
		uint64_t poly_modulus_degree,
		const std::shared_ptr<const FFT>& fft_handler,
		CycloRing::Form form
	) :
		poly_modulus_degree_(poly_modulus_degree),
		coeffs_(std::vector<std::complex<double_t>>(poly_modulus_degree)),
		fft_handler_(fft_handler),
		form_(form)
	{}

	void CycloRing::operator()(const CycloRing& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || form_ != other.form_ || fft_handler_ != other.fft_handler_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree_; ++coeff_idx)
		{
			coeffs_[coeff_idx] = other.coeffs_[coeff_idx];
		}
	}

	const std::complex<double_t>& CycloRing::get(uint64_t coeff_idx) const
	{
		if (coeff_idx >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		return coeffs_[coeff_idx];
	}

	void CycloRing::set(uint64_t coeff_idx, const std::complex<double_t>& value)
	{
		if (coeff_idx >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		coeffs_[coeff_idx] = value;
	}

	uint64_t CycloRing::poly_modulus_degree() const
	{
		return poly_modulus_degree_;
	}

	CycloRing::Form CycloRing::form() const
	{
		return form_;
	}

	void CycloRing::slot_to_coeff(double_t scale)
	{
		if (form_ == CycloRing::Form::coeff)
		{
			throw std::out_of_range("This ring is already coeff form.");
		}

		fft_handler_->inverse_variant_canonical_embedding(coeffs_, scale);

		form_ = CycloRing::Form::coeff;
	}

	void CycloRing::coeff_to_slot(double_t scale)
	{
		if (form_ == CycloRing::Form::slot)
		{
			throw std::out_of_range("This ring is already slot form.");
		}

		fft_handler_->variant_canonical_embedding(coeffs_, scale);

		form_ = CycloRing::Form::slot;
	}
}