
#include "keygenerator.h"
#include "util/distribution.h"


namespace cpet
{
	KeyGenerator::KeyGenerator() :
		sk_(SecretKey())
	{}

	KeyGenerator::KeyGenerator(
		const Context& context
	) :
		context_(context()),
		sk_(SecretKey(context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::coeff, context_->ntt_handler()))
	{
		RnsCycloRing o(1.0, context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::coeff, context_->ntt_handler());
		RnsCycloRing s(1.0, context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::coeff, context_->ntt_handler());

		o.set(0, 1);
		sample_from_hamming_dist(s);

		sk_.set(0, o);
		sk_.set(1, s);
		sk_.coeff_to_slot();
	}

	void KeyGenerator::secret_key(SecretKey& destination) const
	{
		destination = sk_;
	}

	void KeyGenerator::create_public_key(PublicKey& destination) const
	{
		const RnsCycloRing& s = sk_.get(1);
		RnsCycloRing a(1.0, context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::coeff, context_->ntt_handler());
		RnsCycloRing b(1.0, context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::coeff, context_->ntt_handler());
		RnsCycloRing e(1.0, context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::coeff, context_->ntt_handler());

		sample_from_uniform_dist(a);
		sample_from_gaussian_dist(e);

		a.coeff_to_slot();
		e.coeff_to_slot();
		b = a;
		b.negate_inplace();
		b.mul_inplace(s);
		b.add_inplace(e);

		PublicKey pk(context_->poly_modulus_degree(), context_->first_basis(), RnsCycloRing::Form::slot, context_->ntt_handler());
		pk.set(0, b);
		pk.set(1, a);
	}

	void KeyGenerator::create_evaluate_key(SwitchKey& destination) const
	{
		const RnsCycloRing& s = sk_.get(1);
		RnsCycloRing s_sq = s;

		s_sq.mul_inplace(s);

		internal_create_switch_key(s_sq, s, destination);
	}

	void KeyGenerator::internal_create_switch_key(const cpet::RnsCycloRing& s1, const cpet::RnsCycloRing& s2, SwitchKey& destination) const
	{
		if (s1.scale() != s2.scale() || s1.poly_modulus_degree() != s2.poly_modulus_degree() || s1.basis() != s2.basis() || s1.form() != s2.form())
		{
			throw std::out_of_range("Parameter is mismatched.");
		}

		if (s1.basis() != context_->first_basis() || s2.basis() != context_->first_basis())
		{
			throw std::out_of_range("Basis is mismatched with context.");
		}

		Basis basis = context_->first_basis();
		basis.convert_basis(Basis::Type::PQ);

		RnsCycloRing a(1.0, context_->poly_modulus_degree(), basis, RnsCycloRing::Form::coeff, context_->ntt_handler());
		RnsCycloRing b(1.0, context_->poly_modulus_degree(), basis, RnsCycloRing::Form::coeff, context_->ntt_handler());
		RnsCycloRing e(1.0, context_->poly_modulus_degree(), basis, RnsCycloRing::Form::coeff, context_->ntt_handler());
		RnsCycloRing P(1.0, context_->poly_modulus_degree(), basis, RnsCycloRing::Form::coeff, context_->ntt_handler());

		sample_from_uniform_dist(a);
		sample_from_gaussian_dist(e);

		for (uint64_t basis_idx = basis.begin(); basis_idx < basis.end(); ++basis_idx)
		{
			P.set(basis_idx, 0, basis.P_mod_p_and_q(basis_idx));
		}

		a.coeff_to_slot();
		e.coeff_to_slot();
		P.coeff_to_slot();
		b = a;
		b.negate_inplace();
		b.mul_inplace(s2);
		b.add_inplace(e);
		P.mul_inplace(s1);
		b.add_inplace(P);

		SwitchKey swk(context_->poly_modulus_degree(), context_->first_basis(), 2, RnsCycloRing::Form::slot, context_->ntt_handler());
		swk.set(0, b);
		swk.set(1, a);
	}
}