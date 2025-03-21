
#include "cycloring.h"
#include <stdexcept>


namespace cpet
{
	CycloRing::CycloRing() : CycloRing(PolyModulus(), std::make_shared<const FFT>()) {}

	CycloRing::CycloRing(
		const PolyModulus& poly_modulus, 
		const std::shared_ptr<const FFT>& fft_handler, 
		const std::complex<double_t>& value)
	{
		poly_modulus_degree_ = poly_modulus.degree();
		slot_count_ = poly_modulus_degree_ >> 1ULL;
		coeffs_.assign(poly_modulus_degree_, value);
		fft_handler_ = fft_handler;
		ifft_form_ = false;
	}

	const std::complex<double_t>& CycloRing::operator()(uint64_t index) const
	{
		if (index >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		return coeffs_[index];
	}

	void CycloRing::operator()(uint64_t index, const std::complex<double_t>& value)
	{
		if (index >= poly_modulus_degree_)
		{
			throw std::out_of_range("Index out of range.");
		}

		coeffs_[index] = value;
	}

	void CycloRing::assign(
		const PolyModulus& poly_modulus,
		const std::shared_ptr<const FFT>& fft_handler,
		const std::complex<double_t>& value)
	{
		poly_modulus_degree_ = poly_modulus.degree();
		slot_count_ = poly_modulus_degree_ >> 1ULL;
		coeffs_.assign(poly_modulus_degree_, value);
		fft_handler_ = fft_handler;
		ifft_form_ = false;
	}

	uint64_t CycloRing::poly_modulus_degree() const
	{
		return poly_modulus_degree_;
	}

	uint64_t CycloRing::slot_count() const
	{
		return slot_count_;
	}

	void CycloRing::set_ifft_form()
	{
		if (ifft_form_)
		{
			throw std::out_of_range("This ring is already ifft form.");
		}

		fft_handler_->inverse_variant_canonical_embedding(coeffs_);
		ifft_form_ = true;
	}

	void CycloRing::set_normal_form()
	{
		if (!ifft_form_)
		{
			throw std::out_of_range("This ring is already noraml form.");
		}

		fft_handler_->variant_canonical_embedding(coeffs_);
		ifft_form_ = false;
	}
}