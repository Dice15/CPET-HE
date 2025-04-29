#pragma once

#include "polymodulus.h"
#include "fft.h"
#include <cstdint>
#include <complex>
#include <vector>
#include <memory>

namespace cpet
{
	class CycloRing
	{
	public:
		enum class Form : bool { slot, coeff };


	public:
		CycloRing();

		CycloRing(
			uint64_t poly_modulus_degree,
			const std::shared_ptr<const FFT>& fft_handler,
			CycloRing::Form form
		);

		void operator()(const CycloRing& other);

		const std::complex<double_t>& get(uint64_t coeff_idx) const;

		void set(uint64_t coeff_idx, const std::complex<double_t>& value);

		uint64_t poly_modulus_degree() const;

		CycloRing::Form form() const;

		void slot_to_coeff(double_t scale);

		void coeff_to_slot(double_t scale);

	private:
		uint64_t poly_modulus_degree_;

		std::vector<std::complex<double_t>> coeffs_;

		std::shared_ptr<const FFT> fft_handler_;

		CycloRing::Form form_;
	};
}