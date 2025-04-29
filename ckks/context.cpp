
#include "context.h"


namespace cpet
{
	Context::Context() :
		scale_(1.0),
		poly_modulus_degree_(0),
		slot_count_(0),
		moduli_p_bit_sizes_({}),
		moduli_q_bit_sizes_({}),
		max_level_(0),
		bases_({}),
		fft_handler_(nullptr),
		ntt_handler_(nullptr)
	{}

	Context::Context(
		uint64_t scale,
		uint64_t poly_modulus_degree,
		const std::vector<uint64_t>& moduli_p_bit_sizes,
		const std::vector<uint64_t>& moduli_q_bit_sizes
	) :
		scale_(scale),
		poly_modulus_degree_(poly_modulus_degree),
		slot_count_(poly_modulus_degree >> 1ULL),
		moduli_p_bit_sizes_(moduli_p_bit_sizes),
		moduli_q_bit_sizes_(moduli_q_bit_sizes),
		max_level_(moduli_q_bit_sizes.size())
	{
		Basis basis(poly_modulus_degree_ * 2, Basis::Type::Q, moduli_p_bit_sizes_, moduli_q_bit_sizes_);

		bases_.resize(max_level_ + 1);

		for (uint64_t lev = max_level_; lev > 0; --lev)
		{
			bases_[lev] = basis;

			if (lev > 1)
			{
				basis.drop_basis();
			}
		}

		fft_handler_ = std::make_shared<FFT>(FFT(poly_modulus_degree));

		// TODO: NTT에는 그냥 모듈리만 보내는게 나을 듯?
		bases_.back().convert_basis(Basis::Type::PQ);
		ntt_handler_ = std::make_shared<NTT>(NTT(poly_modulus_degree, bases_.back()));
		bases_.back().convert_basis(Basis::Type::Q);

		this_ = std::make_shared<Context>(this);
	}

	const std::shared_ptr<Context>& Context::operator()() const
	{
		return this_;
	}

	double_t Context::scale() const
	{
		return scale_;
	}

	uint64_t Context::poly_modulus_degree() const
	{
		return poly_modulus_degree_;
	}

	uint64_t Context::slot_count() const
	{
		return slot_count_;
	}

	const std::vector<uint64_t>& Context::moduli_p_bit_sizes() const
	{
		return moduli_p_bit_sizes_;
	}

	const std::vector<uint64_t>& Context::moduli_q_bit_sizes() const
	{
		return moduli_q_bit_sizes_;
	}

	uint64_t Context::max_level() const
	{
		return max_level_;
	}

	const Basis& Context::basis(uint64_t level) const
	{
		return bases_[level];
	}

	const Basis& Context::first_basis() const
	{
		return bases_.front();
	}

	const Basis& Context::last_basis() const
	{
		return bases_.back();
	}

	const std::shared_ptr<FFT>& Context::fft_handler() const
	{
		return fft_handler_;
	}

	const std::shared_ptr<NTT>& Context::ntt_handler() const
	{
		return ntt_handler_;
	}
}