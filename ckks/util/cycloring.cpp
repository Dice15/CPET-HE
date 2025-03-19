
#include "cycloring.h"
#include <stdexcept>


namespace cpet
{
	CycloRing::CycloRing() : CycloRing(PolyModulus(), std::make_shared<FFT>()) {}

	CycloRing::CycloRing(const PolyModulus& poly_modulus, const std::shared_ptr<FFT>& fft_handler, const std::complex<double_t>& value)
	{
		poly_modulus_ = poly_modulus;
		slot_count_ = poly_modulus.degree() >> 1ULL;
		coeffs_.assign(poly_modulus_.degree(), value);
		fft_handler_ = fft_handler;
		ifft_form_ = false;
	}

	CycloRing::CycloRing(const CycloRing& other)
	{
		poly_modulus_ = other.poly_modulus_;
		slot_count_ = other.slot_count_;
		coeffs_ = other.coeffs_;
		fft_handler_ = other.fft_handler_;
		ifft_form_ = other.ifft_form_;
	}

	const std::complex<double_t>& CycloRing::operator()(uint64_t index) const
	{
		if (index >= poly_modulus_.degree())
		{
			throw std::out_of_range("Index out of range.");
		}

		return coeffs_[index];
	}

	void CycloRing::operator()(uint64_t index, const std::complex<double_t>& value)
	{
		if (index >= poly_modulus_.degree())
		{
			throw std::out_of_range("Index out of range.");
		}

		coeffs_[index] = value;
	}

	void CycloRing::assign(PolyModulus poly_modulus, const std::complex<double_t>& value)
	{
		poly_modulus_ = poly_modulus;
		slot_count_ = poly_modulus_.degree() / 2;
		coeffs_.assign(poly_modulus_.degree(), value);
		ifft_form_ = false;
	}

	const PolyModulus& CycloRing::poly_modulus() const
	{
		return poly_modulus_;
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