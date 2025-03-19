#pragma once

#include "polymodulus.h"
#include <cstdint>
#include <vector>
#include <complex>


namespace cpet
{
    class FFT
    {
    public:
        FFT() = default;

        FFT(const PolyModulus& poly_modulus);

        void variant_canonical_embedding(std::vector<std::complex<double_t>>& ring) const;

        void inverse_variant_canonical_embedding(std::vector<std::complex<double_t>>& vector) const;

    private:
        void compute_root_of_unity(uint64_t n, uint64_t exponent, std::complex<double_t>& destination) const;

        void compute_inverse_root_of_unity(uint64_t n, uint64_t exponent, std::complex<double_t>& destination) const;

        void fft(std::vector<std::complex<double_t>>& ring, bool inverse = false) const;

        void inverse_fft(std::vector<std::complex<double_t>>& vector) const;

        void fft_negacyclic(std::vector<std::complex<double_t>>& ring) const;

        void inverse_fft_negacyclic(std::vector<std::complex<double_t>>& vector) const;

        double_t PI = 3.14159265358979323846;

        uint64_t poly_modulus_degree_;

        std::vector<uint64_t> bit_reversal_table_;

        std::vector<std::complex<double_t>> zeta_powers_;

        std::vector<std::complex<double_t>> omega_powers_;

        std::vector<std::complex<double_t>> inv_zeta_powers_;

        std::vector<std::complex<double_t>> inv_omega_powers_;
    };
}