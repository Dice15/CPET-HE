#pragma once

#include "polymodulus.h"
#include "basis.h"
#include "ntt.h"
#include <cstdint>
#include <vector>
#include <memory>


namespace cpet
{
	class RnsCycloRing
	{
	public:
		enum class Form : bool { slot, coeff };


	public:
		RnsCycloRing();

		RnsCycloRing(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			RnsCycloRing::Form default_form,
			const std::shared_ptr<const NTT>& ntt_handler
		);

		void operator()(const RnsCycloRing& other);

		uint64_t get(uint64_t coeff_idx) const;

		uint64_t get(uint64_t basis_idx, uint64_t coeff_idx) const;

		void set(uint64_t coeff_idx, uint64_t value);

		void set(uint64_t basis_idx, uint64_t coeff_idx, uint64_t value);

		double_t scale() const;

		uint64_t poly_modulus_degree() const;

		const Basis& basis() const;

		RnsCycloRing::Form form() const;

		void modulus_reduction();

		void rescale();

		void convert_scale_force(double_t scale);

		void convert_basis_force(Basis::Type type);

		void convert_basis_approximate(Basis::Type type);

		void coeff_to_slot();

		void slot_to_coeff();

		void add_inplace(const RnsCycloRing& other);

		void sub_inplace(const RnsCycloRing& other);

		void mul_inplace(const RnsCycloRing& other);

		void mul_inplace(uint64_t value);

		void negate_inplace();


	private:
		void basis_raising_approximate();

		void basis_reduction_approximate();

		double_t scale_;

		uint64_t poly_modulus_degree_;

		Basis basis_;

		std::vector<std::vector<uint64_t>> rns_coeffs_;

		RnsCycloRing::Form form_;

		std::shared_ptr<const NTT> ntt_handler_;
	};
}