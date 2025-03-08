#pragma once

#include "polymodulus.h"
#include <cstdint>
#include <complex>
#include <vector>


namespace cpet
{
	class CycloRing
	{
	public:
		CycloRing();

		CycloRing(PolyModulus poly_modulus, const std::complex<double_t>& value = { 0.0, 0.0 });

		CycloRing(const CycloRing& other);

		const std::complex<double_t>& operator()(uint64_t index) const;

		void operator()(uint64_t index, const std::complex<double_t>& value);

		void assign(PolyModulus poly_modulus_degree, const std::complex<double_t>& value = { 0.0, 0.0 });

		const PolyModulus& poly_modulus() const;

		uint64_t slot_count() const;

		void set_ifft_form();

		void set_normal_form();

	private:
		PolyModulus poly_modulus_;

		uint64_t slot_count_;

		std::vector<std::complex<double_t>> coeffs_;

		bool ifft_form_;
	};
}