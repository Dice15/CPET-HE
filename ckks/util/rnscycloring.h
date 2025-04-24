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
		RnsCycloRing();

		RnsCycloRing(
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler
		);

		RnsCycloRing(
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler,
			int64_t value
		);

		void operator()(const RnsCycloRing& other);

		int64_t get_coeff(uint64_t coeff_idx) const;

		void set_coeff(uint64_t coeff_idx, int64_t value);

		uint64_t get_rns_coeff(uint64_t rns_idx, uint64_t coeff_idx) const;

		void set_rns_coeff(uint64_t rns_idx, uint64_t coeff_idx, uint64_t value);

		uint64_t poly_modulus_degree() const;

		const Basis& get_basis() const;

		void set_basis(const Basis* const basis);

		void set_ntt_form();

		void set_normal_form();

		void add_inplace(const RnsCycloRing& other);

		void sub_inplace(const RnsCycloRing& other);

		void mul_inplace(const RnsCycloRing& other);

		void negate_inplace();

		//void modulus_reduction();

		//void rescale();

	private:
		uint64_t poly_modulus_degree_;

		const Basis* basis_;

		std::vector<std::vector<uint64_t>> rns_coeffs_;

		std::shared_ptr<const NTT> ntt_handler_;

		bool ntt_form_;
	};
}



/*void fast_basis_conversion(
	uint64_t basis_a_size,
	uint64_t basis_b_size,
	std::vector<uint64_t>::const_iterator basis_a,
	std::vector<uint64_t>::const_iterator basis_b,
	std::vector<uint64_t>::const_iterator inv_a_hats_a,
	std::vector<std::vector<uint64_t>>::const_iterator a_hats_b,
	uint64_t congruence_size,
	const std::vector<std::vector<uint64_t>>& congruences_a,
	std::vector<std::vector<uint64_t>>& destination) const;

void approximate_modulus_raising(
	uint64_t basis_p_size,
	uint64_t basis_q_size,
	std::vector<uint64_t>::const_iterator basis_p,
	std::vector<uint64_t>::const_iterator basis_q,
	std::vector<uint64_t>::const_iterator inv_q_hats_q,
	std::vector<std::vector<uint64_t>>::const_iterator q_hats_p,
	uint64_t congruence_size,
	const std::vector<std::vector<uint64_t>>& congruences_q,
	std::vector<std::vector<uint64_t>>& destination) const;

void approximate_modulus_reduction(
	uint64_t basis_p_size,
	uint64_t basis_q_size,
	std::vector<uint64_t>::const_iterator basis_p,
	std::vector<uint64_t>::const_iterator basis_q,
	std::vector<uint64_t>::const_iterator inv_P_q,
	std::vector<uint64_t>::const_iterator inv_p_hats_p,
	std::vector<std::vector<uint64_t>>::const_iterator p_hats_q,
	uint64_t congruence_size,
	const std::vector<std::vector<uint64_t>>& congruences_d,
	std::vector<std::vector<uint64_t>>& destination) const;*/

	//void convert_basis(Basis::basis_type basis_type);