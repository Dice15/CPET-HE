
#include "ciphertext.h"


namespace cpet
{
	Ciphertext::Ciphertext() : Text() {}

	Ciphertext::Ciphertext(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		uint64_t dimension,
		RnsCycloRing::Form form,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(scale, poly_modulus_degree, basis, dimension, form, ntt_handler)
	{}

	void Ciphertext::tensor_with(const Ciphertext& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (dimension() != 2 || other.dimension() != 2)
		{
			throw std::invalid_argument("Only size of text 2 can be tensor product.");
		}

		rings_.resize(3, rings_[1]);

		RnsCycloRing c0_k1 = rings_[0];
		RnsCycloRing& d0 = rings_[0];
		RnsCycloRing& d1 = rings_[1];
		RnsCycloRing& d2 = rings_[2];
	
		// C0 * K0
		d0.mul_inplace(other.rings_[0]);

		// (C0 * K1) + (C1 * K0)
		c0_k1.mul_inplace(other.rings_[1]);
		d1.mul_inplace(other.rings_[0]);
		d1.add_inplace(c0_k1);

		// C1 * K1
		d2.mul_inplace(other.rings_[1]);
		
		scale_ *= other.scale_;
	}
}