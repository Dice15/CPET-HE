#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace cpet
{
	class Basis
	{
	public:
		enum class basis_type : uint16_t;
		class constant;
		class iterator;
		

		enum class basis_type : uint16_t
		{
			basis_d = 0,
			basis_q = 1
		};


		class constant
		{
		public:
			constant(
				uint64_t poly_modulus_degree,
				const std::vector<uint64_t>& basis_p_bit_sizes,
				const std::vector<uint64_t>& basis_q_bit_sizes);

		private:
			uint64_t basis_p_size_;

			uint64_t basis_q_size_;

			uint64_t basis_d_size_;

			std::vector<uint64_t> basis_p_;

			std::vector<uint64_t> basis_q_;

			std::vector<uint64_t> basis_d_;

			std::vector<uint64_t> inv_P_q_l_;

			std::vector<uint64_t> inv_p_hats_p_;

			std::vector<std::vector<uint64_t>> p_hats_q_l_;

			std::vector<std::vector<uint64_t>> inv_q_l_hats_q_l_;

			std::vector<std::vector<std::vector<uint64_t>>> q_l_hats_p_;	

			friend class Basis;

			friend class Basis::iterator;
		};

	
	public:
		Basis() = default;

		Basis(std::shared_ptr<const constant> const_data);

		bool operator==(const Basis& other) const;

		bool operator!=(const Basis& other) const;

		uint64_t operator[](uint64_t index) const;

		uint64_t begin() const;

		uint64_t end() const;

		uint64_t size() const;

		void convert_basis(basis_type basis);

		void pop_basis_q();

		uint64_t basis_p_size() const;

		uint64_t basis_q_size() const;

		uint64_t basis_d_size() const;

		const std::vector<uint64_t>& basis_p() const;

		const std::vector<uint64_t>& basis_q() const;

		const std::vector<uint64_t>& basis_d() const;

		const std::vector<uint64_t>& inv_P_q() const;

		const std::vector<uint64_t>& inv_p_hats_p() const;

		const std::vector<std::vector<uint64_t>>& p_hats_q() const;

		const std::vector<uint64_t>& inv_q_hats_q() const;

		const std::vector<std::vector<uint64_t>>& q_hats_p() const;


	private:
		std::shared_ptr<const constant> const_data_;

		uint64_t pop_count_of_basis_q_;

		basis_type curr_basis_;
	};
}