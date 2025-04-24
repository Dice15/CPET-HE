#pragma once

#include "util/rnscycloring.h"
#include <cstdint>
#include <vector>


namespace cpet
{
	class Text
	{
	public:
		Text();

		Text(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler,
			uint64_t size
		);

		const RnsCycloRing& get_ring(uint64_t slot_idx) const;

		void set_ring(uint64_t slot_idx, const RnsCycloRing& ring);
		
		double_t get_scale() const;

		void set_scale(double_t scale);

		const Basis* get_basis() const;

		void set_basis(Basis basis);

		uint64_t size() const;

		void set_ntt_form();

		void set_normal_form();


	protected:
		double_t scale_;

		const Basis* basis_;

		std::vector<RnsCycloRing> rings_;
	};


	class Plaintext : public Text
	{
	public:
		Plaintext();

		Plaintext(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler
		);

		void add_with(const Plaintext& other);

		void sub_with(const Plaintext& other);

		void tensor_with(const Plaintext& other);
	};


	class Ciphertext : public Text
	{
	public:
		Ciphertext();

		Ciphertext(
			double_t scale,
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler
		);

		void add_with(const Ciphertext& other);

		void sub_with(const Ciphertext& other);

		void tensor_with(const Ciphertext& other);
	};


	class SecretKey : public Text
	{
	public:
		SecretKey();

		SecretKey(
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler
		);
	};


	class PublicKey : public Text
	{
	public:
		PublicKey();

		PublicKey(
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler
		);
	};


	class EvaluateKey : public Text
	{
	public:
		EvaluateKey();

		EvaluateKey(
			uint64_t poly_modulus_degree,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler
		);
	};
}