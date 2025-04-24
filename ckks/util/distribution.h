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
			// [0, q-1] 정수 균등 분포 설정.
			uint64_t min = 0ULL;
			uint64_t max = basis.at(rns_idx) - 1ULL;
			std::uniform_int_distribution<uint64_t> dist(min, max);

			for (uint64_t coeff_idx = 0; coeff_idx < poly_modulus_degree; coeff_idx++)
			{
				destination.set_rns_coeff(rns_idx, coeff_idx, dist(rand));
			}
		}
	}

	// TODO: 이산 가우시안 분포 라이브러리 또는 직접 구현을 해야함. 현재는 임시 코드임.
	void sample_ring_from_gaussian_dist(
		uint64_t poly_modulus_degree,
		const Basis& basis,
		const std::shared_ptr<const NTT>& ntt_handler,
		RnsCycloRing& destination
	) {
		std::random_device rd;
		std::mt19937 rand(rd());
		destination = std::move(RnsCycloRing(poly_modulus_degree, basis, ntt_handler));

		// 표준편차가 3.2인 가우시안 분포로 설정
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

		// -1, +1 중에서 선택하도록 균등 분포 설정. (0을 -1 취급).
		std::uniform_int_distribution<uint64_t> dist(0, 1);

		for (uint64_t rns_idx = basis.begin(); rns_idx < basis.end(); rns_idx++)
		{
			// 다항식 계수 인덱스를 랜덤하게 섞음.
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

		// [0.0, 1.0] 실수 균등 분포 설정.
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