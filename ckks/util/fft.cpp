
#include "fft.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>  

namespace cpet
{
    FFT::FFT(const PolyModulus& poly_modulus) :poly_modulus_degree_(poly_modulus.degree())
    {
        uint64_t poly_modulus_degree_n = poly_modulus_degree_;
        uint64_t poly_modulus_degree_2n = poly_modulus_degree_n << 1;


        // Negacyclic FFT and IFFT (These transformations use 2n-th primitive root(ζ))
        // 
        // X[k] = ∑x[j]*ζ^(j(2k+1)) = ∑x[j]*ζ^(j)*ζ^(2jk), for all k,j ∈ { 0, ... , n-1 }, where n is poly modulus degree
        // 
        // So, compute powers of ζ, ω, ζ⁻¹, ω⁻¹( where ω=ζ²and ω⁻¹=ζ⁻²are root for standard FFT and IFFT).
        zeta_powers_.resize(poly_modulus_degree_n);
        inv_zeta_powers_.resize(poly_modulus_degree_n);
        omega_powers_.resize(poly_modulus_degree_n);
        inv_omega_powers_.resize(poly_modulus_degree_n);

        for (uint64_t j = 0; j < poly_modulus_degree_n; ++j)
        {
            compute_root_of_unity(poly_modulus_degree_2n, j, zeta_powers_[j]);
            compute_inverse_root_of_unity(poly_modulus_degree_2n, j, inv_zeta_powers_[j]);
            compute_root_of_unity(poly_modulus_degree_n, j, omega_powers_[j]);
            compute_inverse_root_of_unity(poly_modulus_degree_n, j, inv_omega_powers_[j]);
        }


        // Create bit-reversal table
        bit_reversal_table_.resize(poly_modulus_degree_n);

        for (uint64_t i = 1, j = 0; i < poly_modulus_degree_n; ++i)
        {
            uint64_t bit = poly_modulus_degree_n >> 1;

            while (j & bit)
            {
                j -= bit;
                bit >>= 1;
            }

            j += bit;

            bit_reversal_table_[i] = j;
        }
    }

    void FFT::variant_canonical_embedding(std::vector<std::complex<double_t>>& ring, double_t scale) const
    {
        if (ring.size() != poly_modulus_degree_)
        {
            throw std::invalid_argument("Ring's poly modulus degree is mismatched with fft handler.");
        }


        fft_negacyclic(ring);


        // [a₀, conj(aₙ₋₁), ... , aₙ₋₁, conj(a₀)] -> [a₀, a₁, ... aₙ₋₁, 0, ... , 0, 0]
        for (uint64_t j = 0; j < poly_modulus_degree_; j += 2)
        {
            ring[j >> 1] = ring[j] / scale;
        }

        for (uint64_t j = (poly_modulus_degree_ >> 1); j < poly_modulus_degree_; ++j)
        {
            ring[j] = 0.0;
        }
    }

    void FFT::inverse_variant_canonical_embedding(std::vector<std::complex<double_t>>& vector, double_t scale) const
    {
        if (vector.size() != poly_modulus_degree_)
        {
            throw std::invalid_argument("Vector's size is mismatched with fft handler.");
        }


        // [a₀, a₁, ... aₙ₋₁, 0, ... , 0, 0] -> [a₀, conj(aₙ₋₁), ... , aₙ₋₁, conj(a₀)]
        for (uint64_t j = (poly_modulus_degree_ >> 1) - 1; j > 0; --j)
        {
            vector[j * 2] = vector[j];
        }

        for (uint64_t j = 0; j < poly_modulus_degree_; j += 2)
        {
            vector[poly_modulus_degree_ - j - 1] = std::conj(vector[j]);
        }

        inverse_fft_negacyclic(vector);

        for (uint64_t i = 0; i < poly_modulus_degree_; i++)
        {
            vector[i] = std::complex<double_t>(vector[i].real() * scale, 0.0);
        }
    }

    void FFT::compute_root_of_unity(uint64_t n, uint64_t exponent, std::complex<double_t>& destination) const
    {
        destination = std::exp(
            std::complex<double_t>(0.0, (2.0 * PI * static_cast<double_t>(exponent)) / static_cast<double_t>(n))
        );
    }

    void FFT::compute_inverse_root_of_unity(uint64_t n, uint64_t exponent, std::complex<double_t>& destination) const
    {
        destination = std::exp(
            std::complex<double_t>(0.0, (-2.0 * PI * static_cast<double_t>(exponent)) / static_cast<double_t>(n))
        );
    }

    void FFT::fft(std::vector<std::complex<double_t>>& ring, bool inverse) const
    {
        const std::vector<std::complex<double_t>>& omega_powers = inverse ? inv_omega_powers_ : omega_powers_;


        // Rearrange ring's coeff using fft table.
        for (uint64_t i = 1; i < poly_modulus_degree_; ++i)
        {
            uint64_t j = bit_reversal_table_[i];

            if (i < j)
            {
                std::swap(ring[i], ring[j]);
            }
        }

        // X[k] = ∑x[j]*ω^(jk), for all k,j ∈ { 0, ... , n-1 }, where n is poly modulus degree.
        // We will use Cooley-Tukey algorithm for FFT.
        for (uint64_t len = 2; len <= poly_modulus_degree_; len <<= 1)
        {
            uint64_t exponent_len = poly_modulus_degree_ / len;

            for (uint64_t i = 0; i < poly_modulus_degree_; i += len)
            {
                uint64_t exponent = 0;

                for (uint64_t j = 0; j < (len >> 1); j++)
                {
                    std::complex<double_t> u = ring[i + j];
                    std::complex<double_t> v = ring[i + j + (len >> 1ULL)] * omega_powers[exponent];

                    ring[i + j] = u + v;
                    ring[i + j + (len >> 1ULL)] = u - v;

                    exponent += exponent_len;
                }
            }
        }
    }

    void FFT::inverse_fft(std::vector<std::complex<double_t>>& vector) const
    {
        // X[k] = 1/n * {∑x[j]*ω^(-jk)}, for all k,j ∈ { 0, ... , n-1 }, where n is poly modulus degree.
        fft(vector, true);


        // Scaling coeff to 1/n * coeff, where n is poly modulus degree.
        double_t n_inv = 1.0 / static_cast<double_t>(poly_modulus_degree_);

        for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
        {
            vector[i] *= n_inv;
        }
    }

    void FFT::fft_negacyclic(std::vector<std::complex<double_t>>& ring) const
    {
        // Negacyclic FFT (This transformation uses 2n-th primitive root(ζ))
        // 
        // X[k] = ∑x[j]*ζ^(j(2k+1)) = ∑x[j]*ζ^(j)*ζ^(2jk), for all k,j ∈ { 0, ... , n-1 }, where n is poly modulus degree
        // 
        // X[0] = x[0]*ζ^(0*1) + x[1]*ζ^(1*1) + x[2]*ζ^(2*1) ... + x[n-1]*ζ^((n-1)*1)
        // X[1] = x[0]*ζ^(0*5) + x[1]*ζ^(1*3) + x[2]*ζ^(2*3) ... + x[n-1]*ζ^((n-1)*3)
        // X[2] = x[0]*ζ^(0*9) + x[1]*ζ^(1*5) + x[2]*ζ^(2*5) ... + x[n-1]*ζ^((n-1)*5)
        // ...
        // X[n-1] = x[0]*ζ^(0*(2n-1)) + x[1]*ζ^(1*(2n-1)) + x[2]*ζ^(2*(2n-1)) ... + x[n-1]*ζ^((n-1)*(2n-1))
        // 
        // So, multiply ζʲ to ring's all coeff and use ζ²as Standard FFT's ω.
        for (uint64_t j = 0; j < poly_modulus_degree_; ++j)
        {
            ring[j] *= zeta_powers_[j];
        }

        fft(ring);
    }

    void FFT::inverse_fft_negacyclic(std::vector<std::complex<double_t>>& vector) const
    {
        // Negacyclic IFFT (This transformation uses 2n-th primitive root(ζ))
        // 
        // x[k] = 1/n * {∑X[j]*ζ^(-j*(2k+1))} = 1/n * {∑X[j]*ζ^(-j)*ζ^(-2jk)}, for all k,j ∈ {0, ... , n-1}, where n is poly modulus degree
        // 
        // x[0] = 1/n * {X[0]*ζ^-(0*1) + X[1]*ζ^-(1*1) + X[2]*ζ^-(2*1) ... + X[n-1]*ζ^-((n-1)*1)}
        // x[1] = 1/n * {X[0]*ζ^-(0*3) + X[1]*ζ^-(1*3) + X[2]*ζ^-(2*3) ... + X[n-1]*ζ^-((n-1)*3)}
        // x[2] = 1/n * {X[0]*ζ^-(0*5) + X[1]*ζ^-(1*5) + X[2]*ζ^-(2*5) ... + X[n-1]*ζ^-((n-1)*5)}
        // ...
        // x[n-1] = 1/n * {X[0]*ζ^-(0*(2n-1)) + X[1]*ζ^-(1*(2n-1)) + X[2]*ζ^-(2*(2n-1)) ... + X[n-1]*ζ^-((n-1)*(2n-1))}
        // 
        // So, use ζ⁻²as standard IFFT's ω and multiply ζ⁻ʲ to ring's all coeff.
        inverse_fft(vector);

        for (uint64_t j = 0; j < poly_modulus_degree_; ++j)
        {
            vector[j] *= inv_zeta_powers_[j];
        }
    }
}