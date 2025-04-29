#pragma once

#include <cstdint>
#include <vector>
#include <memory>


namespace cpet
{
	class Basis
	{
	public:
		enum class Type : bool { PQ, Q };


		struct Constant
		{
			uint64_t p_begin_;

			uint64_t p_end_;

			uint64_t q_begin_;

			uint64_t q_end_;

			std::vector<uint64_t> moduli_;

			std::vector<uint64_t> P_mod_p_and_q_;

			std::vector<uint64_t> inv_P_mod_q_;

			std::vector<uint64_t> inv_p_hats_mod_p;

			std::vector<std::vector<uint64_t>> p_hats_mod_q_;

			std::vector<std::vector<uint64_t>> inv_q_hats_mod_q_by_lev_;

			std::vector<std::vector<std::vector<uint64_t>>> q_hats_mod_p_by_lev;	
		};

	
	public:
		Basis();

		Basis(
			uint64_t poly_modulus_degree,
			Basis::Type default_type,
			const std::vector<uint64_t>& moduli_p_bit_sizes,
			const std::vector<uint64_t>& moduli_q_bit_sizes
		);

		bool operator==(const Basis& other) const;

		bool operator!=(const Basis& other) const;

		bool operator<(const Basis& other) const;

		bool operator>(const Basis& other) const;

		const std::vector<uint64_t>& get_moduli() const;

		uint64_t begin() const;

		uint64_t end() const;

		uint64_t size() const;

		uint64_t level() const;

		void drop_basis();

		void convert_basis(Basis::Type type);

		Basis::Type get_basis_type() const;

		uint64_t capacity() const;

		uint64_t at(uint64_t basis_idx) const;

		uint64_t P_mod_p_and_q(uint64_t basis_idx) const;

		uint64_t inv_P_mod_q(uint64_t basis_q_idx) const;

		uint64_t inv_p_hats_mod_p(uint64_t basis_p_idx) const;

		uint64_t p_hats_mod_q(uint64_t basis_p_idx, uint64_t basis_q_idx) const;

		uint64_t inv_q_hats_mod_q(uint64_t basis_q_idx) const;

		uint64_t q_hats_mod_p(uint64_t basis_q_idx, uint64_t basis_p_idx) const;


	private:
		std::shared_ptr<const Constant> constant_;

		Basis::Type basis_type_;

		uint64_t drop_q_count_;
	};
}