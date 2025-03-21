
#include "ntt.h"
#include "arithmod.h"
#include <stdexcept>
#include <random>


namespace cpet
{
    NTT::NTT(const PolyModulus& poly_modulus, const Basis& basis) :poly_modulus_degree_(poly_modulus.degree()), basis_(basis.basis_d())
    {
        uint64_t basis_size = basis_.size();
        uint64_t poly_modulus_degree_n = poly_modulus_degree_;
        uint64_t poly_modulus_degree_2n = poly_modulus_degree_n << 1;


        // Negacyclic NTT and INTT (These transformations use 2n-th primitive root(ζ))
        // 
        // X[k] = ∑x[j]*ζ^(j(2k+1)) = ∑x[j]*ζ^(j)*ζ^(2jk), for all k,j ∈ { 0, ... , n-1 }, where n is poly modulus degree
        // 
        // So, compute powers of ζ, ω, ζ⁻¹, ω⁻¹( where ω=ζ²and ω⁻¹=ζ⁻²are root for standard NTT and INTT).
        zeta_powers_by_basis_.resize(basis_size, std::vector<uint64_t>(poly_modulus_degree_n));
        inv_zeta_powers_by_basis_.resize(basis_size, std::vector<uint64_t>(poly_modulus_degree_n));
        omega_powers_by_basis_.resize(basis_size, std::vector<uint64_t>(poly_modulus_degree_n));
        inv_omega_powers_by_basis_.resize(basis_size, std::vector<uint64_t>(poly_modulus_degree_n));

        for (uint64_t b = 0; b < basis_size; ++b)
        {
            uint64_t zeta;

            if (!try_minimal_primitive_root(poly_modulus_degree_2n, basis_[b], zeta))
            {
                throw std::logic_error("Failed to find primitive root for NTT negacyclic.");
            }

            uint64_t inv_zeta = inverse_mod(zeta, basis_[b]);
            uint64_t omaga = mul_mod(zeta, zeta, basis_[b]);
            uint64_t inv_omaga = mul_mod(inv_zeta, inv_zeta, basis_[b]);

            zeta_powers_by_basis_[b][0] = 1;
            inv_zeta_powers_by_basis_[b][0] = 1;
            omega_powers_by_basis_[b][0] = 1;
            inv_omega_powers_by_basis_[b][0] = 1;

            for (uint64_t j = 1; j < poly_modulus_degree_n; ++j)
            {
                zeta_powers_by_basis_[b][j] = mul_mod(zeta_powers_by_basis_[b][j - 1], zeta, basis_[b]);
                inv_zeta_powers_by_basis_[b][j] = mul_mod(inv_zeta_powers_by_basis_[b][j - 1], inv_zeta, basis_[b]);
                omega_powers_by_basis_[b][j] = mul_mod(omega_powers_by_basis_[b][j - 1], omaga, basis_[b]);
                inv_omega_powers_by_basis_[b][j] = mul_mod(inv_omega_powers_by_basis_[b][j - 1], inv_omaga, basis_[b]);
            }
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

    bool NTT::is_primitive_root(uint64_t n, uint64_t root, uint64_t prime) const
    {
        if (root >= prime)
        {
            throw std::out_of_range("operand");
        }

        if (root == 0)
        {
            return false;
        }

        // root^(n/2) ≡ -1 ≡ n - 1 (mod p).
        return pow_mod(root, n >> 1ULL, prime) == (prime - 1ULL);
    }

	bool NTT::try_primitive_root(uint64_t n, uint64_t prime, uint64_t& destination) const
	{
        uint64_t group_size = prime - 1ULL;
        uint64_t quotient_group_size = group_size / n;

        // Group size must be divisible by poly modulus degree.
        if (group_size - quotient_group_size * n != 0ULL)
        {
            return false;
        }

        std::random_device rand;
        uint32_t attempt_counter = 0;
        uint32_t attempt_counter_max = 100;

        do
        {
            attempt_counter++;

            // Create random number g on 64-bit integer.
            // if h = g^quotient_group_size ≡ 1 (mod prime), then h is prime n-th root.
            destination = mod((static_cast<uint64_t>(rand()) << 32ULL) | static_cast<uint64_t>(rand()), prime);      
            destination = pow_mod(destination, quotient_group_size, prime);

        } while (!is_primitive_root(n, destination, prime) && (attempt_counter < attempt_counter_max));

        return is_primitive_root(n, destination, prime);
	}

    bool NTT::try_minimal_primitive_root(uint64_t n, uint64_t prime, uint64_t& destination) const
    {
        uint64_t root;

        if (!try_primitive_root(n, prime, root))
        {
            return false;
        }

        uint64_t root_sq = mul_mod(root, root, prime);
        uint64_t candidate_minimal_root = root;

        // For all odd power of roots can be primitive root.
        for (size_t i = 0; i < n; i += 2)
        {
            if (candidate_minimal_root < root)
            {
                root = candidate_minimal_root;
            }

            // Candidate minimal root is always odd power of root.
            candidate_minimal_root = mul_mod(candidate_minimal_root, root_sq, prime);
        }

        destination = root;

        return true;
    }

    void NTT::ntt(std::vector<std::vector<uint64_t>>& rings, uint64_t basis_begin, uint64_t basis_end, bool inverse) const
    {
        const std::vector<std::vector<uint64_t>>& omega_powers = inverse ? inv_omega_powers_by_basis_ : omega_powers_by_basis_;


        // Rearrange ring's coeff using ntt table.
        for (uint64_t b = basis_begin; b < basis_end; ++b)
        {
            for (uint64_t i = 1; i < poly_modulus_degree_; i++)
            {
                uint64_t j = bit_reversal_table_[i];

                if (i < j)
                {
                    std::swap(rings[b][i], rings[b][j]);
                }
            }
        }


        // X[k] = ∑x[j]*ω^(jk), for all k,j ∈ {0, ... , n-1}
        // We will use Cooley-Tukey algorithm for NTT.
        for (uint64_t b = basis_begin; b < basis_end; ++b)
        {
            for (uint64_t len = 2; len <= poly_modulus_degree_; len <<= 1)
            {
                uint64_t exponent_len = poly_modulus_degree_ / len;

                for (uint64_t i = 0; i < poly_modulus_degree_; i += len)
                {
                    uint64_t exponent = 0;

                    for (size_t j = 0; j < (len >> 1); j++)
                    {
                        uint64_t u = rings[b][i + j];
                        uint64_t v = mul_mod(rings[b][i + j + (len >> 1)], omega_powers[b][exponent], basis_[b]);

                        rings[b][i + j] = add_mod(u, v, basis_[b]);
                        rings[b][i + j + (len >> 1)] = sub_mod(u, v, basis_[b]);

                        exponent += exponent_len;
                    }
                }
            }
        }
    }

    void NTT::inverse_ntt(std::vector<std::vector<uint64_t>>& vectors, uint64_t basis_begin, uint64_t basis_end) const
    {
        // X[k] = 1/n * {∑x[j]*ω^(-jk)}, for all k,j ∈ {0, ... , n-1}.
        ntt(vectors, basis_begin, basis_end, true);


        // Scaling coeff to 1/n * coeff, where n is poly modulus degree.
        for (uint64_t b = basis_begin; b < basis_end; ++b)
        {
            uint64_t n_inv = inverse_mod(poly_modulus_degree_, basis_[b]);

            for (size_t i = 0; i < poly_modulus_degree_; i++)
            {
                vectors[b][i] = mul_mod(vectors[b][i], n_inv, basis_[b]);
            }
        }
    }

    void NTT::ntt_negacyclic(std::vector<std::vector<uint64_t>>& rings, uint64_t basis_begin, uint64_t basis_end) const
    {
        // Negacyclic NTT (This transformation uses 2n-th primitive root(ζ))
        // 
        // X[k] = ∑x[j]*ζ^(j(2k+1)) = ∑x[j]*ζ^(j)*ζ^(2jk), for all k,j ∈ { 0, ... , n-1 }, where n is poly modulus degree
        // 
        // X[0] = x[0]*ζ^(0*1) + x[1]*ζ^(1*1) + x[2]*ζ^(2*1) ... + x[n-1]*ζ^((n-1)*1)
        // X[1] = x[0]*ζ^(0*5) + x[1]*ζ^(1*3) + x[2]*ζ^(2*3) ... + x[n-1]*ζ^((n-1)*3)
        // X[2] = x[0]*ζ^(0*9) + x[1]*ζ^(1*5) + x[2]*ζ^(2*5) ... + x[n-1]*ζ^((n-1)*5)
        // ...
        // X[n-1] = x[0]*ζ^(0*(2n-1)) + x[1]*ζ^(1*(2n-1)) + x[2]*ζ^(2*(2n-1)) ... + x[n-1]*ζ^((n-1)*(2n-1))
        // 
        // So, multiply ζʲ to ring's all coeff and use ζ²as Standard NTT's ω.
        for (uint64_t b = basis_begin; b < basis_end; ++b)
        {
            for (uint64_t i = 0; i < poly_modulus_degree_; ++i)
            {
                rings[b][i] = mul_mod(rings[b][i], zeta_powers_by_basis_[b][i], basis_[b]);
            }
        }

        ntt(rings, basis_begin, basis_end);
    }

    void NTT::inverse_ntt_negacyclic(std::vector<std::vector<uint64_t>>& vectors, uint64_t basis_begin, uint64_t basis_end) const
    {
        // Negacyclic INTT (This transformation uses 2n-th primitive root(ζ))
        // 
        // x[k] = 1/n * {∑X[j]*ζ^(-j*(2k+1))} = 1/n * {∑X[j]*ζ^(-j)*ζ^(-2jk)}, for all k,j ∈ {0, ... , n-1}, where n is poly modulus degree
        // 
        // x[0] = 1/n * {X[0]*ζ^-(0*1) + X[1]*ζ^-(1*1) + X[2]*ζ^-(2*1) ... + X[n-1]*ζ^-((n-1)*1)}
        // x[1] = 1/n * {X[0]*ζ^-(0*3) + X[1]*ζ^-(1*3) + X[2]*ζ^-(2*3) ... + X[n-1]*ζ^-((n-1)*3)}
        // x[2] = 1/n * {X[0]*ζ^-(0*5) + X[1]*ζ^-(1*5) + X[2]*ζ^-(2*5) ... + X[n-1]*ζ^-((n-1)*5)}
        // ...
        // x[n-1] = 1/n * {X[0]*ζ^-(0*(2n-1)) + X[1]*ζ^-(1*(2n-1)) + X[2]*ζ^-(2*(2n-1)) ... + X[n-1]*ζ^-((n-1)*(2n-1))}
        // 
        // So, use ζ⁻²as standard INTT's ω and multiply ζ⁻ʲ to ring's all coeff.
        inverse_ntt(vectors, basis_begin, basis_end);

        for (uint64_t b = basis_begin; b < basis_end; ++b)
        {
            for (uint64_t j = 0; j < poly_modulus_degree_; ++j)
            {
                vectors[b][j] = mul_mod(vectors[b][j], inv_zeta_powers_by_basis_[b][j], basis_[b]);
            }
        }
    }
}