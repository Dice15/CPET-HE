#pragma once

#include "polymodulus.h"
#include "modulus.h"
#include <cstdint>
#include <vector>
#include <memory>


namespace cpet
{
	class RnsCycloRing
	{
	public:
		RnsCycloRing();

		RnsCycloRing(PolyModulus poly_modulus, std::shared_ptr<std::vector<Modulus>> moduli, uint64_t value = 0ULL);

		RnsCycloRing(RnsCycloRing&other);

		uint64_t operator()(uint64_t congruence_index, uint64_t coeff_index) const;

		void operator()(uint64_t congruence_index, uint64_t coeff_index, uint64_t value);

		void assign(PolyModulus poly_modulus_degree, std::shared_ptr<std::vector<Modulus>> moduli, uint64_t value = 0ULL);

		const PolyModulus& poly_modulus() const;

		std::shared_ptr<std::vector<Modulus>> moduli() const;

		void set_ntt_form();

		void set_normal_form();

		void add(const RnsCycloRing& other, RnsCycloRing&destination) const;

		void sub(const RnsCycloRing& other, RnsCycloRing& destination) const;

		void mul(const RnsCycloRing& other, RnsCycloRing& destination) const;

		void negate(RnsCycloRing& destination) const;

	private:
		PolyModulus poly_modulus_;

		//ModulusChain moduli_begin_;

		//ModulusChain moduli_end_;

		std::shared_ptr<std::vector<Modulus>> moduli_;

		std::vector<std::vector<uint64_t>> congruences_;

		bool ntt_form_;
	};
}