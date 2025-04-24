
#include "text.h"
#include <stdexcept>


namespace cpet
{
	Text::Text() :
		scale_(1.0),
		basis_(nullptr),
		rings_({})
	{}

	Text::Text(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		uint64_t size
	) :
		scale_(scale),
		basis_(&basis),
		rings_(std::vector<RnsCycloRing>(size, RnsCycloRing(poly_modulus_degree, basis, ntt_handler)))
	{}

	const RnsCycloRing& Text::get_ring(uint64_t slot_idx) const
	{
		if (rings_.size() <= slot_idx)
		{
			throw std::out_of_range("Out of range.");
		}

		return rings_[slot_idx];
	}

	void Text::set_ring(uint64_t slot_idx, const RnsCycloRing& ring)
	{
		if (rings_.size() <= slot_idx)
		{
			throw std::out_of_range("Out of range.");
		}

		if (*basis_ != ring.get_basis())
		{
			throw std::out_of_range("Basis is mismatched.");
		}

		rings_[slot_idx](ring);
	}

	double_t Text::get_scale() const
	{
		return scale_;
	}

	void Text::set_scale(double_t scale)
	{
		scale_ = scale;
	}

	const Basis* Text::get_basis() const
	{
		return basis_;
	}

	void Text::set_basis(Basis basis)
	{
		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].set_basis(basis_);
		}
	}

	uint64_t Text::size() const
	{
		return rings_.size();
	}

	void Text::set_ntt_form()
	{
		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].set_ntt_form();
		}
	}

	void Text::set_normal_form()
	{
		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].set_normal_form();
		}
	}



	/// <summary>
	/// ///////////////////
	/// </summary>



	Plaintext::Plaintext() : Text() {}

	Plaintext::Plaintext(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(scale, poly_modulus_degree, basis, ntt_handler, 1)
	{}

	void Plaintext::add_with(const Plaintext& other)
	{
		if (scale_ != other.scale_)
		{
			throw std::invalid_argument("Scale of texts must be matched.");
		}

		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].add_inplace(other.rings_[slot_idx]);
		}
	}

	void Plaintext::sub_with(const Plaintext& other)
	{
		if (scale_ != other.scale_)
		{
			throw std::invalid_argument("Scale of texts must be matched.");
		}

		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].sub_inplace(other.rings_[slot_idx]);
		}
	}

	void Plaintext::tensor_with(const Plaintext& other)
	{
		rings_[0].mul_inplace(other.rings_[0]);
		scale_ *= other.scale_;
	}



	/// <summary>
	/// ///////////////////
	/// </summary>



	Ciphertext::Ciphertext() : Text() {}

	Ciphertext::Ciphertext(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(scale, poly_modulus_degree, basis, ntt_handler, 2)
	{}

	void Ciphertext::add_with(const Ciphertext& other)
	{
		if (scale_ != other.scale_)
		{
			throw std::invalid_argument("Scale of texts must be matched.");
		}

		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].add_inplace(other.rings_[slot_idx]);
		}
	}

	void Ciphertext::sub_with(const Ciphertext& other)
	{
		if (scale_ != other.scale_)
		{
			throw std::invalid_argument("Scale of texts must be matched.");
		}

		for (uint64_t slot_idx = 0; slot_idx < rings_.size(); ++slot_idx)
		{
			rings_[slot_idx].sub_inplace(other.rings_[slot_idx]);
		}
	}

	void Ciphertext::tensor_with(const Ciphertext& other)
	{
		if (rings_.size() != 2 || other.rings_.size() != 2)
		{
			throw std::invalid_argument("Only size of text 2 can be tensor product.");
		}

		RnsCycloRing ring_0 = rings_[0];

		rings_.resize(3);
		rings_[0].mul_inplace(other.rings_[0]);
		rings_[1].mul_inplace(other.rings_[0]);
		ring_0.mul_inplace(other.rings_[1]);
		rings_[1].add_inplace(ring_0);
		rings_[2].mul_inplace(other.rings_[0]);

		scale_ *= other.scale_;
	}



	/// <summary>
	/// ///////////////////
	/// </summary>



	SecretKey::SecretKey() : Text() {}

	SecretKey::SecretKey(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(1.0, poly_modulus_degree, basis, ntt_handler, 2)
	{}



	/// <summary>
	/// ///////////////////
	/// </summary>



	PublicKey::PublicKey() : Text() {}

	PublicKey::PublicKey(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(1.0, poly_modulus_degree, basis, ntt_handler, 2)
	{}



	/// <summary>
	/// ///////////////////
	/// </summary>



	EvaluateKey::EvaluateKey() : Text() {}

	EvaluateKey::EvaluateKey(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		Text(1.0, poly_modulus_degree, basis, ntt_handler, 2)
	{}
}