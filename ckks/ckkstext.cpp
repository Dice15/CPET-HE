
#include "ckkstext.h"
#include <stdexcept>


namespace cpet
{
	CKKStext::CKKStext() :
		row_size_(0), col_size_(0), matrix_(std::vector<std::vector<RnsCycloRing>>()) {}

	CKKStext::CKKStext(
		double_t scale,
		uint64_t row_size,
		uint64_t col_size,
		const PolyModulus& poly_modulus,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler)
	{
		if (basis.curr_basis() != Basis::basis_type::basis_q)
		{
			throw std::out_of_range("CKKStext's default basis must be basis q.");
		}

		row_size_ = row_size;
		col_size_ = col_size;
		matrix_.assign(row_size, std::vector<RnsCycloRing>(col_size, RnsCycloRing(poly_modulus, basis, scale, ntt_handler)));
	}

	int64_t CKKStext::operator()(uint64_t row, uint64_t col, uint64_t index) const
	{
		if (row >= row_size_ || col >= col_size_)
		{
			throw std::out_of_range("Out of range.");
		}

		return matrix_[row][col](index);
	}

	void CKKStext::operator()(uint64_t row, uint64_t col, uint64_t index, int64_t value)
	{
		if (row >= row_size_ || col >= col_size_)
		{
			throw std::out_of_range("Out of range.");
		}

		matrix_[row][col](index, value);
	}

	double_t CKKStext::get_scale() const
	{
		if (row_size_ == 0 || col_size_ == 0)
		{
			return 1.0;
		}

		return matrix_[0][0].get_scale();
	}

	void CKKStext::set_scale(double_t scale)
	{
		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].set_scale(scale);
			}
		}
	}

	uint64_t CKKStext::row_size() const
	{
		return row_size_;
	}

	uint64_t CKKStext::col_size() const
	{
		return col_size_;
	}

	void CKKStext::set_ntt_form()
	{
		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].set_ntt_form();
			}
		}
	}

	void CKKStext::set_normal_form()
	{
		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].set_normal_form();
			}
		}
	}

	void CKKStext::add(const CKKStext& other)
	{
		if (row_size_ != other.row_size_ || col_size_ != other.col_size_)
		{
			throw std::invalid_argument("Matrixs must have matching dimension.");
		}

		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].add_inplace(other.matrix_[r][c]);
			}
		}
	}

	void CKKStext::sub(const CKKStext& other)
	{
		if (row_size_ != other.row_size_ || col_size_ != other.col_size_)
		{
			throw std::invalid_argument("Matrixs must have matching dimension.");
		}

		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].sub_inplace(other.matrix_[r][c]);
			}
		}
	}

	void CKKStext::mul(const CKKStext& other)
	{
		if (col_size_ != other.row_size_)
		{
			throw std::invalid_argument("Matrixs must have matching dimension.");
		}

		std::vector<std::vector<RnsCycloRing>> result(row_size_);
		RnsCycloRing temp = matrix_[0][0];

		for (uint64_t r = 0; r < row_size_; r++)
		{
			result[r].reserve(other.col_size_);

			for (uint64_t k = 0; k < other.col_size_; k++)
			{
				matrix_[r][0].mul(other.matrix_[0][k], temp);
				result[r].emplace_back(std::move(temp));

				for (uint64_t c = 1; c < col_size_; c++)
				{
					matrix_[r][c].mul(other.matrix_[c][k], temp);
					result[r][k].add_inplace(temp);
				}
			}
		}

		matrix_ = std::move(result);
	}

	void CKKStext::modulus_reduction()
	{
		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].modulus_reduction();
			}
		}
	}

	void CKKStext::rescale()
	{
		for (uint64_t r = 0; r < row_size_; r++)
		{
			for (uint64_t c = 0; c < col_size_; c++)
			{
				matrix_[r][c].rescale();
			}
		}
	}
}