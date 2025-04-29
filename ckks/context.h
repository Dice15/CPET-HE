#pragma once

#include "util/basis.h"
#include "util/fft.h"
#include "util/ntt.h"
#include <cstdint>
#include <vector>
#include <memory>


namespace cpet
{
	class Context
	{
	public:
		Context();

		Context(
			uint64_t scale,
			uint64_t poly_modulus_degree,
			const std::vector<uint64_t>& moduli_p_bit_sizes,
			const std::vector<uint64_t>& moduli_q_bit_sizes
		);

		const std::shared_ptr<Context>& operator()() const;

		double_t scale() const;

		uint64_t poly_modulus_degree() const;

		uint64_t slot_count() const;

		const std::vector<uint64_t>& moduli_p_bit_sizes() const;

		const std::vector<uint64_t>& moduli_q_bit_sizes() const;

		uint64_t max_level() const;

		const Basis& basis(uint64_t level) const;

		const Basis& first_basis() const;

		const Basis& last_basis() const;

		const std::shared_ptr<FFT>& fft_handler() const;

		const std::shared_ptr<NTT>& ntt_handler() const;


	private:
		double_t scale_;

		uint64_t poly_modulus_degree_;

		uint64_t slot_count_;

		std::vector<uint64_t> moduli_p_bit_sizes_;

		std::vector<uint64_t> moduli_q_bit_sizes_;

		uint64_t max_level_;

		std::vector<Basis> bases_;

		std::shared_ptr<FFT> fft_handler_;

		std::shared_ptr<NTT> ntt_handler_;

		std::shared_ptr<Context> this_;
	};
}