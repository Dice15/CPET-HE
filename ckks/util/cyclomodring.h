#pragma once

#include "polymodulus.h"
#include "modulus.h"
#include <cstdint>
#include <vector>


namespace cpet
{
	class CycloModRing
	{
	public:
		CycloModRing();

		CycloModRing(PolyModulus poly_modulus, Modulus modulus, uint64_t value = 0ULL);

		CycloModRing(CycloModRing &other);

		uint64_t operator()(uint64_t index) const;

		void operator()(uint64_t index, uint64_t value);

		void assign(PolyModulus poly_modulus_degree, Modulus modulus, uint64_t value = 0ULL);

		const PolyModulus& poly_modulus() const;

		const Modulus& modulus() const;

		void set_ntt_form();

		void set_normal_form();

		void add(const CycloModRing& other, CycloModRing &destination) const;

		void sub(const CycloModRing& other, CycloModRing& destination) const;

		void mul(const CycloModRing& other, CycloModRing& destination) const;

		void negate(CycloModRing& destination) const;

	private:
		PolyModulus poly_modulus_;

		Modulus modulus_;

		std::vector<uint64_t> coeffs_;

		bool ntt_form_;
	};
}