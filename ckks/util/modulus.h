#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace cpet
{
	class Basis
	{
	public:
		class ConstData
		{
		public:
			ConstData();

			ConstData(
				uint64_t poly_modulus_degree,
				const std::vector<uint64_t>& basis_p_bit_sizes,
				const std::vector<uint64_t>& basis_q_bit_sizes);

			uint64_t basis_p_size() const;

			uint64_t basis_q_size() const;

			std::vector<uint64_t>::const_iterator basis_p() const;

			std::vector<uint64_t>::const_iterator basis_q() const;

			std::vector<uint64_t>::const_iterator inv_P_q() const;

			std::vector<uint64_t>::const_iterator inv_p_hats_p() const;

			std::vector<std::vector<uint64_t>>::const_iterator p_hats_q() const;

			std::vector<std::vector<uint64_t>>::const_iterator inv_q_l_hats_q() const;

			std::vector<std::vector<std::vector<uint64_t>>>::const_iterator q_l_hats_p() const;

		private:
			std::vector<uint64_t> basis_p_;

			std::vector<uint64_t> basis_q_;

			std::vector<uint64_t> inv_P_q_;

			std::vector<uint64_t> inv_p_hats_p_;

			std::vector<std::vector<uint64_t>> p_hats_q_;

			std::vector<std::vector<uint64_t>> inv_q_l_hats_q_;

			std::vector<std::vector<std::vector<uint64_t>>> q_l_hats_p_;	
		};

		enum class basis_type : uint16_t
		{
			basis_d = 0,
			basis_q = 1
		};

		enum class iterator_type : uint16_t
		{
			begin = 0,
			end = 1
		};

		class basis_iterator
		{
		public:
			basis_iterator(
				std::shared_ptr<const ConstData> const_data,
				basis_type basis,
				uint64_t pop_modulus_of_basis_q_count,
				iterator_type iter);
			
			uint64_t operator-(basis_iterator other) const;

			basis_iterator operator+(uint64_t n) const;

			basis_iterator& operator++();

			basis_iterator operator++(int);

			basis_iterator& operator--();

			basis_iterator operator--(int);

			bool operator==(const basis_iterator& other) const;

			bool operator!=(const basis_iterator& other) const;

			const uint64_t& operator*() const;

		private:
			std::shared_ptr<const Basis::ConstData> const_data_;

			uint64_t begin_index_;
			
			uint64_t end_index_;

			uint64_t curr_index_;
		};
	
	public:
		Basis();

		Basis(std::shared_ptr<const ConstData> const_data);

		Basis(const Basis& basis);

		bool operator==(const Basis& other) const;

		bool operator!=(const Basis& other) const;

		uint64_t size() const;

		basis_iterator begin() const;

		basis_iterator end() const;

		void convert_basis(basis_type basis);

		basis_type current_basis() const;

		void pop_modulus_of_basis_q();

		uint64_t basis_p_size() const;

		uint64_t basis_q_size() const;

		std::vector<uint64_t>::const_iterator basis_p() const;

		std::vector<uint64_t>::const_iterator basis_q() const;

		std::vector<uint64_t>::const_iterator inv_P_q() const;

		std::vector<uint64_t>::const_iterator inv_p_hats_p() const;

		std::vector<std::vector<uint64_t>>::const_iterator p_hats_q() const;

		std::vector<uint64_t>::const_iterator inv_q_hats_q() const;

		std::vector<std::vector<uint64_t>>::const_iterator q_hats_p() const;


	private:
		std::shared_ptr<const ConstData> const_data_;

		uint64_t pop_modulus_of_basis_q_count_;

		basis_type curr_basis_;
	};
}