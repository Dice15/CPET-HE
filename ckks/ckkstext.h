#pragma once

#include "util/rnscycloring.h"
#include <cstdint>
#include <vector>


namespace cpet
{
	class CKKStext
	{
	public:
		CKKStext();

		CKKStext(
			double_t scale,
			uint64_t row_size,
			uint64_t col_size,
			const PolyModulus& poly_modulus,
			const Basis& basis,
			const std::shared_ptr<const NTT>& ntt_handler);

		int64_t operator()(uint64_t row, uint64_t col, uint64_t index) const;

		void operator()(uint64_t row, uint64_t col, uint64_t index, int64_t value);
		
		double_t get_scale() const;

		void set_scale(double_t scale);

		uint64_t row_size() const;

		uint64_t col_size() const;

		void set_ntt_form();

		void set_normal_form();

		void add(const CKKStext& other);

		void sub(const CKKStext& other);

		void mul(const CKKStext& other);

		void modulus_reduction();

		void rescale();

	private:
		uint64_t row_size_;

		uint64_t col_size_;

		std::vector<std::vector<RnsCycloRing>> matrix_;
	};
}
