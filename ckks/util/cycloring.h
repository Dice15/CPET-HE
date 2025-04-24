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
		CycloRing() = default;

		CycloRing(
			const PolyModulus& poly_modulus, 
			const std::shared_ptr<const FFT>& fft_handler, 
			const std::complex<double_t>& value = { 0.0, 0.0 });

		const std::complex<double_t>& operator()(uint64_t index) const;

		void operator()(uint64_t index, const std::complex<double_t>& value);

		void assign(
			const PolyModulus& poly_modulus,
			const std::shared_ptr<const FFT>& fft_handler,
			const std::complex<double_t>& value = { 0.0, 0.0 });

		uint64_t poly_modulus_degree() const;

		void set_ifft_form(double_t scale);

		void set_normal_form(double_t scale);

	private:
		uint64_t poly_modulus_degree_;

		std::vector<std::complex<double_t>> coeffs_;

		std::shared_ptr<const FFT> fft_handler_;

		bool ifft_form_;
	};
}