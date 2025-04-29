#pragma once

#include "text.h"
#include <cstdint>


namespace cpet
{
	class Ciphertext : public Text
	{
	public:
		Ciphertext();

		Ciphertext(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			uint64_t dimension,
			RnsCycloRing::Form form,
			const std::shared_ptr<const NTT>& ntt_handler
		);

		void tensor_with(const Ciphertext& other);
	};
}