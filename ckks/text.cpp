
#include "text.h"
#include <stdexcept>


namespace cpet
{
	Text::Text() :
		scale_(0.0),
		poly_modulus_degree_(0),
		basis_(Basis()),
		rings_({}),
		form_(RnsCycloRing::Form::coeff),
		ntt_handler_(nullptr)
	{}

	Text::Text(
		double_t scale,
		uint64_t poly_modulus_degree,
		const Basis& basis,
		uint64_t dimension,
		RnsCycloRing::Form default_form,
		const std::shared_ptr<const NTT>& ntt_handler
	) :
		scale_(scale),
		poly_modulus_degree_(poly_modulus_degree),
		basis_(basis),
		rings_(std::vector<RnsCycloRing>(dimension, RnsCycloRing(scale, poly_modulus_degree, basis, default_form, ntt_handler))),
		form_(default_form),
		ntt_handler_(ntt_handler)
	{}

	const RnsCycloRing& Text::get(uint64_t index) const
	{
		if (dimension() <= index)
		{
			throw std::out_of_range("Index out of range.");
		}

		return rings_[index];
	}

	void Text::set(uint64_t index, const RnsCycloRing& ring)
	{
		if (dimension() <= index)
		{
			throw std::out_of_range("Index out of range.");
		}

		if (scale_ != ring.scale() || poly_modulus_degree_ != ring.poly_modulus_degree() || basis_ != ring.basis() || form_ != ring.form())
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		// TODO
	}

	double_t Text::scale() const
	{
		return scale_;
	}

	uint64_t Text::poly_modulus_degree() const
	{
		return poly_modulus_degree_;
	}

	const Basis& Text::basis() const
	{
		return basis_;
	}

	uint64_t Text::dimension() const
	{
		return rings_.size();
	}

	RnsCycloRing::Form Text::form() const
	{
		return form_;
	}

	void Text::coeff_to_slot()
	{
		for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
		{
			rings_[dim_idx].coeff_to_slot();
		}
	}

	void Text::modulus_reduction()
	{
		for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
		{
			rings_[dim_idx].modulus_reduction();
		}

		basis_.drop_basis();
	}

	void Text::rescale()
	{
		for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
		{
			rings_[dim_idx].rescale();
		}

		scale_ /= basis_.at(basis_.end() - 1);
		basis_.drop_basis();
	}

	void Text::convert_scale_force(double_t scale)
	{
		scale_ = scale;
	}

	void Text::convert_basis_force(Basis::Type type)
	{
		for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
		{
			rings_[dim_idx].convert_basis_force(type);
		}

		basis_.convert_basis(type);
	}

	void Text::convert_basis_approximate(Basis::Type type)
	{
		for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
		{
			rings_[dim_idx].convert_basis_approximate(type);
		}

		basis_.convert_basis(type);
	}

	void Text::slot_to_coeff()
	{
		for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
		{
			rings_[dim_idx].slot_to_coeff();
		}
	}

	void Text::add_inplace(const Text& other)
	{
		if (scale_ != other.scale_ || poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (dimension() == other.dimension())
		{
			for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
			{
				rings_[dim_idx].add_inplace(other.rings_[dim_idx]);
			}
		}
		else if (dimension() < other.dimension())
		{
			rings_.resize(other.dimension(), RnsCycloRing(scale_, poly_modulus_degree_, basis_, form_, ntt_handler_));

			for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
			{
				rings_[dim_idx].add_inplace(other.rings_[dim_idx]);
			}
		}
		else
		{
			for (uint64_t dim_idx = 0; dim_idx < other.dimension(); ++dim_idx)
			{
				rings_[dim_idx].add_inplace(other.rings_[dim_idx]);
			}
		}
	}

	void Text::sub_inplace(const Text& other)
	{
		if (scale_ != other.scale_ || poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (dimension() == other.dimension())
		{
			for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
			{
				rings_[dim_idx].sub_inplace(other.rings_[dim_idx]);
			}
		}
		else if (dimension() < other.dimension())
		{
			rings_.resize(other.dimension(), RnsCycloRing(scale_, poly_modulus_degree_, basis_, form_, ntt_handler_));

			for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
			{
				rings_[dim_idx].sub_inplace(other.rings_[dim_idx]);
			}
		}
		else
		{
			for (uint64_t dim_idx = 0; dim_idx < other.dimension(); ++dim_idx)
			{
				rings_[dim_idx].sub_inplace(other.rings_[dim_idx]);
			}
		}
	}

	void Text::mul_inplace(const Text& other)
	{
		if (poly_modulus_degree_ != other.poly_modulus_degree_ || basis_ != other.basis_ || form_ != other.form_)
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (dimension() == other.dimension())
		{
			for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
			{
				rings_[dim_idx].mul_inplace(other.rings_[dim_idx]);
			}
		}
		else if (dimension() < other.dimension())
		{
			rings_.resize(other.dimension(), RnsCycloRing(scale_, poly_modulus_degree_, basis_, form_, ntt_handler_));

			for (uint64_t dim_idx = 0; dim_idx < dimension(); ++dim_idx)
			{
				rings_[dim_idx].mul_inplace(other.rings_[dim_idx]);
			}
		}
		else
		{
			for (uint64_t dim_idx = 0; dim_idx < other.dimension(); ++dim_idx)
			{
				rings_[dim_idx].mul_inplace(other.rings_[dim_idx]);
			}
		}

		scale_ *= other.scale_;
	}
}