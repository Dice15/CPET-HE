#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace cpet
{
	class Basis
	{
	public:
		enum class BasisType : uint16_t
		{
			basis_q = 0,
			basis_pq = 1,
		};


		struct Constant
		{
			uint64_t p_begin_;

			uint64_t p_end_;

			uint64_t q_begin_;

			uint64_t q_end_;

			std::vector<uint64_t> moduli_;

			std::vector<uint64_t> P_mod_ql_;

			std::vector<uint64_t> inv_P_mod_ql_;

			std::vector<uint64_t> inv_p_hats_mod_p;

			std::vector<std::vector<uint64_t>> p_hats_mod_ql_;

			std::vector<std::vector<uint64_t>> inv_ql_hats_mod_ql_;

			std::vector<std::vector<std::vector<uint64_t>>> ql_hats_mod_p;	
		};

	
	public:
		Basis();

		Basis(
			uint64_t poly_modulus_degree,
			const std::vector<uint64_t>& moduli_p_bit_sizes,
			const std::vector<uint64_t>& moduli_q_bit_sizes
		);

		bool operator==(const Basis& other) const;

		bool operator!=(const Basis& other) const;

		BasisType get_type() const;

		void set_type(BasisType basis_type);

		const std::vector<uint64_t>& get_moduli() const;

		uint64_t begin() const;

		uint64_t end() const;

		uint64_t size() const;

		uint64_t capacity() const;

		uint64_t at(uint64_t index) const;

		void drop_q();

		const std::vector<uint64_t>& P_mod_ql() const;

		const std::vector<uint64_t>& inv_P_mod_ql() const;

		const std::vector<uint64_t>& inv_p_hats_mod_p() const;

		const std::vector<std::vector<uint64_t>>& p_hats_mod_ql() const;

		const std::vector<uint64_t>& inv_ql_hats_mod_ql() const;

		const std::vector<std::vector<uint64_t>>& ql_hats_mod_p() const;


	private:
		std::shared_ptr<const Constant> constant_;

		BasisType basis_type_;

		uint64_t drop_q_count_;
	};
}