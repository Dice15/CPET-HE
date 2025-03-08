
#include "fft.h"
#include "arithmod.h"
#include <stdexcept>
#include <cmath>
#include <complex>
#include <iostream>

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace cpet
{
	void compute_minimal_primitive_root_powers(const PolyModulus& poly_modulus, std::vector<std::complex<double_t>>& destination)
	{
        // root^(n/2) = -1, where n is poly modulus degree.
        destination.reserve(poly_modulus.degree());

        for (uint64_t i = 0; i < poly_modulus.degree(); i++)
        {
            destination.push_back(std::exp(
                std::complex<double_t>(0.0, (2.0 * PI * i) / static_cast<double_t>(poly_modulus.degree()))
            ));
        }
	}

    void compute_inverse_minimal_primitive_root_powers(const PolyModulus& poly_modulus, std::vector<std::complex<double_t>>& destination)
    {
        // root^(n/2) = -1, where n is poly modulus degree.
        destination.reserve(poly_modulus.degree());

        for (uint64_t i = 0; i < poly_modulus.degree(); i++)
        {
            destination.push_back(std::exp(
                std::complex<double_t>(0.0, (-2.0 * PI * i) / static_cast<double_t>(poly_modulus.degree()))
            ));
        }
    }

    void create_fft_table(const PolyModulus& poly_modulus, std::vector<uint64_t>& destination)
    {
        // Assign fft table
        destination.assign(poly_modulus.degree(), 0ULL);


        // Bit-reversal
        for (uint64_t i = 1, j = 0; i < poly_modulus.degree(); i++)
        {
            uint64_t bit = poly_modulus.degree() >> 1;

            while (j & bit)
            {
                j -= bit;
                bit >>= 1;
            }

            j += bit;

            destination[i] = j;
        }
    }

    void fft(CycloRing& ring, const std::vector<std::complex<double_t>>& omega_powers)
    {
        // Get ring data.
        const PolyModulus& poly_modulus = ring.poly_modulus();


        // Rearrange ring's coeff using fft table.
        std::vector<uint64_t> fft_table;
        create_fft_table(poly_modulus, fft_table);

        for (uint64_t i = 1; i < poly_modulus.degree(); i++)
        {
            uint64_t j = fft_table[i];

            if (i < j)
            {
                std::complex<double_t> coeff_i = ring(i);
                ring(i, ring(j));
                ring(j, coeff_i);
            }
        }
 
        // X[k] = ¢²x[j]*¥ø^(jk), for all k,j ¡ô {0, ... , n-1}, where n is poly modulus degree.
        // We will use Cooley-Tukey algorithm for FFT.
        for (uint64_t len = 2; len <= poly_modulus.degree(); len <<= 1)
        {
            uint64_t exponent_len = poly_modulus.degree() / len;

            for (uint64_t i = 0; i < poly_modulus.degree(); i += len)
            {
                uint64_t exponent = 0;

                for (uint64_t j = 0; j < (len >> 1); j++)
                {
                    std::complex<double_t> u = ring(i + j);
                    std::complex<double_t> v = ring(i + j + (len >> 1ULL)) * omega_powers[exponent];

                    ring(i + j, u + v);
                    ring(i + j + (len >> 1ULL), u - v);

                    exponent += exponent_len;
                }
            }
        }
    }

    void inverse_fft(CycloRing& ring, const std::vector<std::complex<double_t>>& inverse_omega_powers)
    {
        // Get ring data.
        const PolyModulus& poly_modulus = ring.poly_modulus();


        // X[k] = 1/n * {¢²x[j]*¥ø^(-jk)}, for all k,j ¡ô {0, ... , n-1}, where n is poly modulus degree.
        // So, use ¥ø^-1 as Standard NTT's ¥ø for IFFT.
        fft(ring, inverse_omega_powers);


        // Scaling coeff to 1/n * coeff, where n is poly modulus degree.
        double inv_n = 1.0 / static_cast<double_t>(ring.poly_modulus().degree());

        for (uint64_t i = 0; i < ring.poly_modulus().degree(); i++)
        {
            ring(i, ring(i) * inv_n);
        }
    }

    void fft_negacyclic(CycloRing& ring)
    {
        // Get ring data.
        const PolyModulus& poly_modulus_n = ring.poly_modulus();


        // Find 2n-th primitive root(¥æ). 
        std::vector<std::complex<double_t>> zeta_powers;
        const PolyModulus poly_modulus_2n(2 * poly_modulus_n.degree());

        compute_minimal_primitive_root_powers(poly_modulus_2n, zeta_powers);


        // Negacyclic FFT
        // X[k] = ¢²x[j]*¥æ^(j(2k+1)) = ¢²x[j]*¥æ^j*(¥æ^2)^jk, for all k,j ¡ô {0, ... , n-1}
        // 
        // X[0] = x[0]*¥æ^(0*1) + x[1]*¥æ^(1*1) + x[2]*¥æ^(2*1) ... + x[n-1]*¥æ^((n-1)*1)
        // X[1] = x[0]*¥æ^(0*3) + x[1]*¥æ^(1*3) + x[2]*¥æ^(2*3) ... + x[n-1]*¥æ^((n-1)*3)
        // X[2] = x[0]*¥æ^(0*5) + x[1]*¥æ^(1*5) + x[2]*¥æ^(2*5) ... + x[n-1]*¥æ^((n-1)*5)
        // ...
        // X[n-1] = x[0]*¥æ^(0*(2n-1)) + x[1]*¥æ^(1*(2n-1)) + x[2]*¥æ^(2*(2n-1)) ... + x[n-1]*¥æ^((n-1)*(2n-1))
        // 
        // So, multiply ¥æ^j to ring's all coeff and use ¥æ©÷as Standard FFT's ¥ø.
        uint64_t zeta_exponent = 0;

        for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
        {
            ring(j, ring(j) * zeta_powers[zeta_exponent]);
            zeta_exponent = zeta_exponent + 1;
        }

        std::vector<std::complex<double_t>> omega_powers(poly_modulus_n.degree());

        for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
        {
            omega_powers[j] = zeta_powers[j * 2];
        }

        fft(ring, omega_powers);
    }

    void inverse_fft_negacyclic(CycloRing& ring)
    {
        // Get ring data.
        const PolyModulus& poly_modulus_n = ring.poly_modulus();


        // Find 2n-th primitive root(¥æ).
        std::vector<std::complex<double_t>> inverse_zeta_powers;
        const PolyModulus poly_modulus_2n(2 * poly_modulus_n.degree());
        
        compute_inverse_minimal_primitive_root_powers(poly_modulus_2n, inverse_zeta_powers);

        //std::cout << "2n-th primitive root for fft negacyclic: " << inverse_zeta_powers[1] << "\n";

        // Negacyclic IFFT
        // x[k] = 1/n * {¢²X[j]*¥æ^(-j(2k+1))} = 1/n * {¢²X[j]*¥æ^-j*(¥æ^2)^jk}, for all k,j ¡ô {0, ... , n-1}
        // 
        // x[0] = 1/n * {X[0]*¥æ^-(0*1) + X[1]*¥æ^-(1*1) + X[2]*¥æ^-(2*1) ... + X[n-1]*¥æ^-((n-1)*1)}
        // x[1] = 1/n * {X[0]*¥æ^-(0*3) + X[1]*¥æ^-(1*3) + X[2]*¥æ^-(2*3) ... + X[n-1]*¥æ^-((n-1)*3)}
        // x[2] = 1/n * {X[0]*¥æ^-(0*5) + X[1]*¥æ^-(1*5) + X[2]*¥æ^-(2*5) ... + X[n-1]*¥æ^-((n-1)*5)}
        // ...
        // x[n-1] = 1/n * {X[0]*¥æ^-(0*(2n-1)) + X[1]*¥æ^-(1*(2n-1)) + X[2]*¥æ^-(2*(2n-1)) ... + X[n-1]*¥æ^-((n-1)*(2n-1))}
        // 
        // So, use ¥æ©÷as Standard IFFT's ¥ø and multiply ¥æ^-j to ring's all coeff.
        std::vector<std::complex<double_t>> inverse_omega_powers(poly_modulus_n.degree());

        for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
        {
            inverse_omega_powers[j] = inverse_zeta_powers[j * 2];
        }

        inverse_fft(ring, inverse_omega_powers);

        uint64_t inv_zeta_exponent = 0;

        for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
        {
            ring(j, ring(j) * inverse_zeta_powers[inv_zeta_exponent]);
            inv_zeta_exponent = inv_zeta_exponent + 1;
        }
    }

    void variant_canonical_embedding(CycloRing& ring)
    {
        // Get ring data.
        const PolyModulus& poly_modulus = ring.poly_modulus();

        fft_negacyclic(ring);
    }

    void inverse_variant_canonical_embedding(CycloRing& ring)
    {
        // Get ring data.
        const PolyModulus& poly_modulus = ring.poly_modulus();

        for (uint64_t j = 0; j < poly_modulus.degree() / 2; j++)
        {
            ring(poly_modulus.degree() - 1 - j, std::conj(ring(j)));
        }

        inverse_fft_negacyclic(ring);

        for (uint64_t i = 0; i < poly_modulus.degree(); i++)
        {
            ring(i, std::complex<double_t>(ring(i).real(), 0.0));
        }

    }
}