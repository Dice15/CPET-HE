
#include "switchkey.h"


namespace cpet
{
	SwitchKey::SwitchKey() : Text() {}

	SwitchKey::SwitchKey(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		uint64_t dimension,
		RnsCycloRing::Form form,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(1.0, poly_modulus_degree, basis, dimension, form, ntt_handler)
	{}
}