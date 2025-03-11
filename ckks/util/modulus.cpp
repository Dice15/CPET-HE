
#include "modulus.h"
#include "numeric.h"
#include <stdexcept>


namespace cpet
{
	// Modulus
    Modulus::Modulus() :Modulus(2) {}

    Modulus::Modulus(uint64_t modulus)
    {
        if (is_over_60_bit(modulus) || !is_prime(modulus) || modulus < 2)
        {
            throw std::invalid_argument("The modulus must be a prime number greater than or equal to 2 and less than or equal to 60 bits.");
        }

        modulus_ = modulus;
    }

	uint64_t Modulus::value() const
    {
        return modulus_;
    }


	// Modulus chain iterator
	ModulusChain::iterator::iterator(
		std::shared_ptr<const std::vector<uint64_t>> moduli, uint64_t index, uint64_t rbegin_index, uint64_t rend_index)
		: moduli_(moduli), index_(index), rbegin_index_(rbegin_index), rend_index_(rend_index) {}

	ModulusChain::iterator& ModulusChain::iterator::operator++()
	{
		if (index_ >= rend_index_)
		{
			throw std::out_of_range("Iterator cannot be incremented beyond the valid range");
		}

		index_++;
		return *this;
	}

	ModulusChain::iterator ModulusChain::iterator::operator++(int) 
	{ 
		iterator temp = *this;
		++(*this);
		return temp;
	}

	ModulusChain::iterator& ModulusChain::iterator::operator--()
	{ 
		if (index_ <= rbegin_index_)
		{
			throw std::out_of_range("Iterator cannot be decremented below the valid range");
		}

		index_--;
		return *this; 
	}

	ModulusChain::iterator ModulusChain::iterator::operator--(int)
	{
		iterator temp = *this;
		--(*this);
		return temp;
	}

	bool ModulusChain::iterator::operator==(const iterator& other) const
	{ 
		return index_ == other.index_;
	}

	bool ModulusChain::iterator::operator!=(const iterator& other) const
	{
		return index_ != other.index_;
	}

	const uint64_t& ModulusChain::iterator::operator*() const
	{
		if (index_ < rbegin_index_ || index_ >= rend_index_)
		{
			throw std::out_of_range("Dereferencing iterator out of valid range");
		}

		return moduli_->at(index_);
	}

	const uint64_t* ModulusChain::iterator::operator->() const
	{ 
		if (index_ < rbegin_index_ || index_ >= rend_index_)
		{
			throw std::out_of_range("Accessing iterator out of valid range");
		}

		return &moduli_->at(index_); 
	}


	// Modulus chain
	ModulusChain::ModulusChain() :
		moduli_(std::make_shared<std::vector<uint64_t>>()), rbegin_index_(0), rend_index_(0) {}

	ModulusChain::ModulusChain(
		uint64_t poly_modulus_degree, const std::vector<uint64_t>& modulus_bit_sizes, const std::vector<uint64_t>& modulus_counts)
	{
		if (modulus_bit_sizes.size() != modulus_counts.size())
		{
			throw std::invalid_argument("Modulus's bit sizes and counts must have the same size.");
		}
	
		std::vector<uint64_t> moduli;
		std::vector<uint64_t> new_moduli;
		uint64_t chain_size = 0;

		for (const auto count : modulus_counts)
		{
			chain_size += count;
		}

		moduli.reserve(chain_size);

		for (uint64_t i = 0; i < modulus_counts.size(); i++)
		{
			try_primes(poly_modulus_degree * 2, modulus_bit_sizes[i], modulus_counts[i], new_moduli);

			moduli.insert(moduli.end(), new_moduli.begin(), new_moduli.end());
		}

		moduli_ = std::make_shared<std::vector<uint64_t>>(moduli);
		rbegin_index_ = 0;
		rend_index_ = moduli_->size();
	}

	ModulusChain::ModulusChain(const ModulusChain& modulus_chain) :
		moduli_(modulus_chain.moduli_), rbegin_index_(modulus_chain.rbegin_index_), rend_index_(modulus_chain.rend_index_) {}

	uint64_t ModulusChain::size() const
	{
		return moduli_->size();
	}

	void ModulusChain::set_range(uint64_t begin_index, uint64_t end_index)
	{
		if (begin_index < 0 || end_index > moduli_->size() || begin_index > end_index)
		{
			throw std::invalid_argument("Out of range");
		}

		rbegin_index_ = begin_index;
		rend_index_ = end_index;
	}

	uint64_t ModulusChain::rsize() const
	{
		return rend_index_ - rbegin_index_;
	}

	ModulusChain::iterator ModulusChain::rbegin() const
	{
		return iterator(moduli_, rbegin_index_, rbegin_index_, rend_index_);
	}

	ModulusChain::iterator ModulusChain::rend() const
	{
		return iterator(moduli_, rend_index_, rbegin_index_, rend_index_);
	}
}