
#include "cyclomodring.h"
#include "arithmod.h"
#include "ntt.h"
#include <stdexcept>


namespace cpet
{
	CycloModRing::CycloModRing() : CycloModRing(0, Modulus()) {}

	CycloModRing::CycloModRing(PolyModulus poly_modulus, Modulus modulus, uint64_t value)
	{
		poly_modulus_ = poly_modulus;
		modulus_ = modulus;
		coeffs_.assign(poly_modulus_.degree(), mod(value, modulus_.value()));
		ntt_form_ = false;
	}

	CycloModRing::CycloModRing(CycloModRing& other)
	{
		poly_modulus_ = other.poly_modulus_;
		modulus_ = other.modulus_;
		coeffs_ = other.coeffs_;
		ntt_form_ = other.ntt_form_;
	}

	uint64_t CycloModRing::operator()(uint64_t index) const
	{
		if (index >= coeffs_.size())
		{
			throw std::out_of_range("Index out of range.");
		}

		return coeffs_[index];
	}

	void CycloModRing::operator()(uint64_t index, uint64_t value)
	{
		if (index >= coeffs_.size())
		{
			throw std::out_of_range("Index out of range.");
		}

		coeffs_[index] = mod(value, modulus_.value());
	}

	void CycloModRing::assign(PolyModulus poly_modulus, Modulus modulus, uint64_t value)
	{
		poly_modulus_ = poly_modulus;
		modulus_ = modulus;
		coeffs_.assign(poly_modulus_.degree(), mod(value, modulus_.value()));
		ntt_form_ = false;
	}

	const PolyModulus& CycloModRing::poly_modulus() const
	{
		return poly_modulus_;
	}

	const Modulus& CycloModRing::modulus() const
	{
		return modulus_;
	}

	void CycloModRing::set_ntt_form()
	{
		if (ntt_form_)
		{
			throw std::out_of_range("This ring is already ntt form.");
		}

		ntt_negacyclic(*this);

		ntt_form_ = true;
	}

	void CycloModRing::set_normal_form()
	{
		if (!ntt_form_)
		{
			throw std::out_of_range("This ring is already normal form.");
		}

		inverse_ntt_negacyclic(*this);

		ntt_form_ = false;
	}

	void CycloModRing::add(const CycloModRing& other, CycloModRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || modulus_.value() != other.modulus_.value())
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.modulus_ = modulus_;
		destination.coeffs_.resize(coeffs_.size());
		destination.ntt_form_ = ntt_form_;

		for (uint64_t i = 0; i < destination.coeffs_.size(); i++)
		{
			destination.coeffs_[i] = add_mod(coeffs_[i], other.coeffs_[i], modulus_.value());
		}
	}

	void CycloModRing::sub(const CycloModRing& other, CycloModRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || modulus_.value() != other.modulus_.value())
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.modulus_ = modulus_;
		destination.coeffs_.resize(coeffs_.size());
		destination.ntt_form_ = ntt_form_;

		for (uint64_t i = 0; i < poly_modulus_.degree(); i++)
		{
			destination.coeffs_[i] = sub_mod(coeffs_[i], other.coeffs_[i], modulus_.value());
		}
	}

	void CycloModRing::mul(const CycloModRing& other, CycloModRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || modulus_.value() != other.modulus_.value())
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.modulus_ = modulus_;
		destination.coeffs_.resize(coeffs_.size());
		destination.ntt_form_ = ntt_form_;

		for (uint64_t i = 0; i < poly_modulus_.degree(); i++)
		{
			destination.coeffs_[i] = mul_mod(coeffs_[i], other.coeffs_[i], modulus_.value());
		}
	}

	void CycloModRing::negate(CycloModRing& destination) const
	{
		destination.poly_modulus_ = poly_modulus_;
		destination.modulus_ = modulus_;
		destination.coeffs_.resize(coeffs_.size());
		destination.ntt_form_ = ntt_form_;

		for (uint64_t i = 0; i < poly_modulus_.degree(); i++)
		{
			destination.coeffs_[i] = negate_mod(coeffs_[i], modulus_.value());
		}
	}
}