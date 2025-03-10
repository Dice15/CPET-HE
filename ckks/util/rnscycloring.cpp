
#include "rnscycloring.h"
#include "arithmod.h"
#include "ntt.h"
#include <stdexcept>


namespace cpet
{
	RnsCycloRing::RnsCycloRing() : RnsCycloRing(0, std::shared_ptr<std::vector<Modulus>>()) {}

	RnsCycloRing::RnsCycloRing(PolyModulus poly_modulus, std::shared_ptr<std::vector<Modulus>> moduli, uint64_t value)
	{
		poly_modulus_ = poly_modulus;
		moduli_ = moduli;
		congruences_.assign(moduli_->size(), std::vector<uint64_t>(poly_modulus_.degree()));

		for (uint64_t m = 0; m < moduli_->size(); m++)
		{
			for (auto& coeff : congruences_[m])
			{
				coeff = mod(value, moduli_->at(m).value());
			}
		}

		ntt_form_ = false;
	}

	RnsCycloRing::RnsCycloRing(RnsCycloRing& other)
	{
		poly_modulus_ = other.poly_modulus_;
		moduli_ = other.moduli_;
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

		congruences_[congruence_index][coeff_index] = mod(value, moduli_->at(congruence_index).value());
	}

	void RnsCycloRing::assign(PolyModulus poly_modulus, std::shared_ptr<std::vector<Modulus>> moduli, uint64_t value)
	{
		poly_modulus_ = poly_modulus;
		moduli_ = moduli;
		congruences_.assign(moduli_->size(), std::vector<uint64_t>(poly_modulus_.degree()));

		for (uint64_t m = 0; m < moduli_->size(); m++)
		{
			for (auto& coeff : congruences_[m])
			{
				coeff = mod(value, moduli_->at(m).value());
			}
		}

		ntt_form_ = false;
	}

	const PolyModulus& RnsCycloRing::poly_modulus() const
	{
		return poly_modulus_;
	}

	std::shared_ptr<std::vector<Modulus>> RnsCycloRing::moduli() const
	{
		return moduli_;
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
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || moduli_ != other.moduli_)
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.moduli_ = moduli_;
		destination.congruences_.resize(congruences_.size(), std::vector<uint64_t>(poly_modulus_.degree()));
		destination.ntt_form_ = ntt_form_;

		for (uint64_t m = 0; m < congruences_.size(); m++)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); i++)
			{
				destination.congruences_[m][i] = add_mod(congruences_[m][i], other.congruences_[m][i], moduli_->at(m).value());
			}
		}
	}

	void RnsCycloRing::sub(const RnsCycloRing& other, RnsCycloRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || moduli_ != other.moduli_)
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.moduli_ = moduli_;
		destination.congruences_.resize(congruences_.size(), std::vector<uint64_t>(poly_modulus_.degree()));
		destination.ntt_form_ = ntt_form_;

		for (uint64_t m = 0; m < congruences_.size(); m++)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); i++)
			{
				destination.congruences_[m][i] = sub_mod(congruences_[m][i], other.congruences_[m][i], moduli_->at(m).value());
			}
		}
	}

	void RnsCycloRing::mul(const RnsCycloRing& other, RnsCycloRing& destination) const
	{
		if (poly_modulus_.degree() != other.poly_modulus_.degree() || moduli_ != other.moduli_)
		{
			throw std::invalid_argument("Rings must have matching degree, modulus and reduction polynomial.");
		}

		if (!ntt_form_ || !other.ntt_form_)
		{
			throw std::invalid_argument("Ring of CKKS scheme must be in NTT form.");
		}

		destination.poly_modulus_ = poly_modulus_;
		destination.moduli_ = moduli_;
		destination.congruences_.resize(congruences_.size(), std::vector<uint64_t>(poly_modulus_.degree()));
		destination.ntt_form_ = ntt_form_;

		for (uint64_t m = 0; m < congruences_.size(); m++)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); i++)
			{
				destination.congruences_[m][i] = mul_mod(congruences_[m][i], other.congruences_[m][i], moduli_->at(m).value());
			}
		}
	}

	void RnsCycloRing::negate(RnsCycloRing& destination) const
	{
		destination.poly_modulus_ = poly_modulus_;
		destination.moduli_ = moduli_;
		destination.congruences_.resize(congruences_.size());
		destination.ntt_form_ = ntt_form_;

		for (uint64_t m = 0; m < congruences_.size(); m++)
		{
			for (uint64_t i = 0; i < congruences_[m].size(); i++)
			{
				destination.congruences_[m][i] = negate_mod(congruences_[m][i], moduli_->at(m).value());
			}
		}
	}
}