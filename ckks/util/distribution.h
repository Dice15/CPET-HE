#pragma once

#include "rnscycloring.h"
#include <vector>
#include <cstdint>
#include <random>
#include <numeric>


namespace cpet
{
	void sample_ring_from_uniform_dist(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		RnsCycloRing& destination
	) {
		std::random_device rd;
		std::mt19937 rand(rd());
		destination = std::move(RnsCycloRing(poly_modulus_degree, basis, ntt_handler));

		for (uint64_t rns_idx = basis.begin(); rns_idx < basis.end(); rns_idx++)
		{
			// [0, q-1] ���� �յ� ���� ����.
			uint64_t min = 0ULL;
			uint64_t max = basis.at(rns_idx) - 1ULL;
			std::uniform_int_distribution<uint64_t> dist(min, max);

			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++)
			{
				destination.set_rns_coeff(rns_idx, coeff_idx, dist(rand));
			}
		}
	}

	// TODO: �̻� ����þ� ���� ���̺귯�� �Ǵ� ���� ������ �ؾ���. ����� �ӽ� �ڵ���.
	void sample_ring_from_gaussian_dist(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		RnsCycloRing& destination
	) {
		std::random_device rd;
		std::mt19937 rand(rd());
		destination = std::move(RnsCycloRing(poly_modulus_degree, basis, ntt_handler));

		// ǥ�������� 3.2�� ����þ� ������ ����
		double_t e = 0.0;
		double_t sd = 3.2;
		std::normal_distribution<double_t> dist(e, sd);

		for (uint64_t rns_idx = basis.begin(); rns_idx < basis.end(); rns_idx++)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++)
			{
				double_t r = dist(rand);

				if (std::signbit(r))
				{
					r = negate_mod(static_cast<uint64_t>(std::abs(std::llround(r))), basis.at(rns_idx));
				}

				destination.set_rns_coeff(rns_idx, coeff_idx, r);
			}
		}
	}

	void sample_ring_from_hamming_dist(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		RnsCycloRing& destination
	) {
		std::random_device rd;
		std::mt19937 rand(rd());
		destination = std::move(RnsCycloRing(poly_modulus_degree, basis, ntt_handler));

		// h = 64
		uint64_t h = std::min(64ULL, poly_modulus_degree);

		// -1, +1 �߿��� �����ϵ��� �յ� ���� ����. (0�� -1 ���).
		std::uniform_int_distribution<uint64_t> dist(0, 1);

		for (uint64_t rns_idx = basis.begin(); rns_idx < basis.end(); rns_idx++)
		{
			// ���׽� ��� �ε����� �����ϰ� ����.
			std::vector<uint64_t> coeff_idx(poly_modulus_degree);
			std::iota(coeff_idx.begin(), coeff_idx.end(), 0);   // [0, 1, ..., n - 1]
			std::shuffle(coeff_idx.begin(), coeff_idx.end(), rand);

			for (uint64_t i = 0; i < h; i++)
			{
				uint64_t r = dist(rand);

				if (r == 0ULL)
				{
					r = basis.at(rns_idx) - 1ULL;
				}

				destination.set_rns_coeff(rns_idx, coeff_idx[i], r);
			}
		}
	}

	void sample_ring_from_zero_one_dist(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		RnsCycloRing& destination
	) {
		std::random_device rd;
		std::mt19937 rand(rd());
		destination = std::move(RnsCycloRing(poly_modulus_degree, basis, ntt_handler));

		// [0.0, 1.0] �Ǽ� �յ� ���� ����.
		std::uniform_real_distribution<double_t> dist(0.0, 1.0);

		for (uint64_t rns_idx = basis.begin(); rns_idx < basis.end(); rns_idx++)
		{
			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++)
			{
				uint64_t r = dist(rand);

				if (r < 0.25)
				{
					r = basis.at(rns_idx) - 1ULL;
				}
				else if (r < 0.5)
				{
					r = 1ULL;
				}
				else
				{
					r= 0ULL;
				}

				destination.set_rns_coeff(rns_idx, coeff_idx, r);
			}
		}
	}
}