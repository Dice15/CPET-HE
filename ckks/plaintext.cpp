
#include "plaintext.h"


namespace cpet
{
	Plaintext::Plaintext() : Text() {}

	Plaintext::Plaintext(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		RnsCycloRing::Form form,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(scale, poly_modulus_degree, basis, 1, form, ntt_handler)
	{}
}