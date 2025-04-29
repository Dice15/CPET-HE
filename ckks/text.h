#pragma once

#include "util/rnscycloring.h"
#include <cstdint>
#include <vector>


namespace cpet
{
	class Text
	{
	public:
		Text();

		Text(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			uint64_t dimension,
			RnsCycloRing::Form default_form,
			const std::shared_ptr<const NTT>& ntt_handler
		);

		const RnsCycloRing& get(uint64_t index) const;

		void set(uint64_t index, const RnsCycloRing& ring);
		
		double_t scale() const;

		uint64_t poly_modulus_degree() const;

		const Basis& basis() const;

		uint64_t dimension() const;

		RnsCycloRing::Form form() const;

		void modulus_reduction();

		void rescale();

		void convert_scale_force(double_t scale);

		void convert_basis_force(Basis::Type type);

		void convert_basis_approximate(Basis::Type type);

		void coeff_to_slot();

		void slot_to_coeff();

		void add_inplace(const Text& other);

		void sub_inplace(const Text& other);

		void mul_inplace(const Text& other);


	protected:
		double_t scale_;

		uint64_t poly_modulus_degree_;

		Basis basis_;

		std::vector<RnsCycloRing> rings_;

		RnsCycloRing::Form form_;

		std::shared_ptr<const NTT> ntt_handler_;
	};
}