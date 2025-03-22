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
		RnsCycloRing() = delete;

		RnsCycloRing(
			const PolyModulus& poly_modulus, 
			const Basis& basis, 
			double_t scale,
			const std::shared_ptr<const NTT>& ntt_handler, 
			uint64_t value = 0ULL);

		int64_t operator()(uint64_t index) const;

		void operator()(uint64_t index, int64_t value);

		void assign(
			const PolyModulus& poly_modulus,
			const Basis& basis,
			double_t scale,
			const std::shared_ptr<const NTT>& ntt_handler,
			uint64_t value = 0ULL);

		uint64_t poly_modulus_degree() const;

		double_t get_scale() const;

		void set_scale(double_t scale);

		const Basis& basis() const; // 삭제 필요

		void set_ntt_form();

		void set_normal_form();

		void add_inplace(const RnsCycloRing& other);

		void add(const RnsCycloRing& other, RnsCycloRing& destination) const;

		void sub_inplace(const RnsCycloRing& other);

		void sub(const RnsCycloRing& other, RnsCycloRing& destination) const;

		void mul_inplace(const RnsCycloRing& other);

		void mul(const RnsCycloRing& other, RnsCycloRing& destination) const;

		void negate_inplace();

		void modulus_reduction();

		void rescale();


	private:
		uint64_t poly_modulus_degree_;

		Basis basis_;

		double_t scale_;

		std::vector<std::vector<uint64_t>> congruences_;

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