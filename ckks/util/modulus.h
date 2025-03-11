#pragma once

#include <cstdint>
#include <vector>
#include <memory>

namespace cpet
{
	class Modulus
	{
	public:
		Modulus();

		Modulus(uint64_t modulus);

		uint64_t value() const;

	private:
		uint64_t modulus_;
	};

	class ModulusChain 
	{
	public:
		class iterator
		{
		public:
			iterator(std::shared_ptr<const std::vector<uint64_t>> moduli, uint64_t index, uint64_t rbegin_index, uint64_t rend_index);

			iterator& operator++();

			iterator operator++(int);

			iterator& operator--();

			iterator operator--(int);

			bool operator==(const iterator& other) const;

			bool operator!=(const iterator& other) const;

			const uint64_t& operator*() const;

			const uint64_t* operator->() const;

		private:
			std::shared_ptr<const std::vector<uint64_t>> moduli_;

			uint64_t index_;

			uint64_t rbegin_index_;
			
			uint64_t rend_index_;
		};

	public:
		ModulusChain();

		ModulusChain(uint64_t poly_modulus_degree, const std::vector<uint64_t>& modulus_bit_sizes, const std::vector<uint64_t>& modulus_counts);

		ModulusChain(const ModulusChain& modulus_chain);

		uint64_t size() const;

		void set_range(uint64_t begin, uint64_t end);

		uint64_t rsize() const;
		
		iterator rbegin() const;

		iterator rend() const;

	private:
		std::shared_ptr<const std::vector<uint64_t>> moduli_;

		uint64_t rbegin_index_;

		uint64_t rend_index_;
	};
}