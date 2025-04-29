#pragma once

#include "text.h"
#include <cstdint>


namespace cpet
{
	class Plaintext : public Text
	{
	public:
		Plaintext();

		Plaintext(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			RnsCycloRing::Form form,
			const std::shared_ptr<const NTT>& ntt_handler
		);
	};
}