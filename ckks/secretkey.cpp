
#include "secretkey.h"


namespace cpet
{
	SecretKey::SecretKey() : Text() {}

	SecretKey::SecretKey(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		RnsCycloRing::Form form,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(1.0, poly_modulus_degree, basis, 2, form, ntt_handler)
	{}
}