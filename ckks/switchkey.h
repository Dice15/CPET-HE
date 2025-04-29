#pragma once

#include "text.h"
#include <cstdint>


namespace cpet
{
	class SwitchKey : public Text
	{
	public:
		SwitchKey();

		SwitchKey(
			uint64_t poly_modulus_degree,
			const Basis& basis,
			uint64_t dimension,
			RnsCycloRing::Form form,
			const std::shared_ptr<const NTT>& ntt_handler
		);
	};
}