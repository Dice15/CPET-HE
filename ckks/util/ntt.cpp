
#include "ntt.h"
#include "arithmod.h"
#include <stdexcept>
#include <random>
#include <iostream>


namespace cpet
{
	bool is_primitive_root(uint64_t root, const PolyModulus& poly_modulus, const Modulus& modulus)
	{
        if (root >= modulus.value())
        {
            throw std::out_of_range("operand");
        }

        if (root == 0)
        {
            return false;
        }

        // root^(n/2) ≡ -1 ≡ n - 1 (mod p), where n is poly modulus degree and p is modulus.
        return pow_mod(root, poly_modulus.degree() >> 1ULL, modulus.value()) == (modulus.value() - 1);
	}

	bool try_primitive_root(const PolyModulus& poly_modulus, const Modulus& modulus, uint64_t& destination)
	{
        uint64_t group_size = modulus.value() - 1;
        uint64_t quotient_group_size = group_size / poly_modulus.degree();

        // Group size must be divisible by poly modulus degree.
        if (group_size - quotient_group_size * poly_modulus.degree() != 0)
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
            // if h = g^quotient_group_size ≡ 1 (mod modulus), then h is modulus degree-th root.
            destination = mod((static_cast<uint64_t>(rand()) << 32) | static_cast<uint64_t>(rand()), modulus.value());      
            destination = pow_mod(destination, quotient_group_size, modulus.value());

        } while (!is_primitive_root(destination, poly_modulus.degree(), modulus) && (attempt_counter < attempt_counter_max));

        return is_primitive_root(destination, poly_modulus, modulus);
	}

    bool try_minimal_primitive_root(const PolyModulus& poly_modulus, const Modulus& modulus, uint64_t& destination)
    {
        uint64_t root;

        if (!try_primitive_root(poly_modulus, modulus, root))
        {
            return false;
        }

        uint64_t root_sq = mul_mod(root, root, modulus.value());
        uint64_t candidate_minimal_root = root;

        // For all odd power of roots can be primitive root.
        for (size_t i = 0; i < poly_modulus.degree(); i += 2)
        {
            if (candidate_minimal_root < root)
            {
                root = candidate_minimal_root;
            }

            // Candidate minimal root is always odd power of root.
            candidate_minimal_root = mul_mod(candidate_minimal_root, root_sq, modulus.value());
        }

        destination = root;

        return true;
    }

    bool try_inverse_minimal_primitive_root(const PolyModulus& poly_modulus, const Modulus& modulus, uint64_t& destination)
    {
        uint64_t root;

        try_minimal_primitive_root(poly_modulus, modulus, root);

        destination = inverse_mod(root, modulus.value());

        return true;
    }

    void create_ntt_table(const PolyModulus& poly_modulus, std::vector<uint64_t>& destination)
    {
        // Assign ntt table
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

    void ntt(RnsCycloRing& ring, const std::vector<uint64_t>& omegas)
    {
        // Get ring data.
        const auto& poly_modulus = ring.poly_modulus();
        const auto& moduli = *ring.moduli();


        // Rearrange ring's coeff using ntt table.
        std::vector<uint64_t> ntt_table;
        create_ntt_table(poly_modulus, ntt_table);

        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            for (uint64_t i = 1; i < poly_modulus.degree(); i++)
            {
                uint64_t j = ntt_table[i];

                if (i < j)
                {
                    uint64_t coeff = ring(m, i);
                    ring(m, i, ring(m, j));
                    ring(m, j, coeff);

                }
            }
        }


        // X[k] = ∑x[j]*ω^(jk), for all k,j ∈ {0, ... , n-1}
        // We will use Cooley-Tukey algorithm for NTT.
        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            for (uint64_t len = 2; len <= poly_modulus.degree(); len <<= 1)
            {
                uint64_t wlen = pow_mod(omegas[m], poly_modulus.degree() / len, moduli[m].value());

                for (size_t i = 0; i < poly_modulus.degree(); i += len)
                {
                    uint64_t w = 1;

                    for (size_t j = 0; j < (len >> 1); j++)
                    {
                        uint64_t u = ring(m, i + j);
                        uint64_t v = mul_mod(ring(m, i + j + (len >> 1)), w, moduli[m].value());

                        ring(m, i + j, add_mod(u, v, moduli[m].value()));
                        ring(m, i + j + (len >> 1), sub_mod(u, v, moduli[m].value()));

                        w = mul_mod(w, wlen, moduli[m].value());
                    }
                }
            }
        }    
    }

    void inverse_ntt(RnsCycloRing& ring, const std::vector<uint64_t>& inv_omegas)
    {
        // Get ring data.
        const auto& poly_modulus = ring.poly_modulus();
        const auto& moduli = *ring.moduli();


        // X[k] = 1/n * {∑x[j]*ω^(-jk)}, for all k,j ∈ {0, ... , n-1}
        // So, use ω^-1 as Standard NTT's ω for INTT.
        ntt(ring, inv_omegas);


        // Scaling coeff to 1/n * coeff
        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            uint64_t n_inv = inverse_mod(poly_modulus.degree(), moduli[m].value());

            for (size_t i = 0; i < poly_modulus.degree(); i++)
            {
                ring(m, i, mul_mod(ring(m, i), n_inv, moduli[m].value()));
            }
        }
    }

    void ntt_negacyclic(RnsCycloRing& ring)
    {
        // Get ring data.
        const auto& poly_modulus_n = ring.poly_modulus();
        const auto& moduli = *ring.moduli();


        // Find 2n-th primitive root(ζ). 
        std::vector<uint64_t> zetas(moduli.size());
        const PolyModulus poly_modulus_2n(2 * poly_modulus_n.degree());

        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            if (!try_minimal_primitive_root(poly_modulus_2n, moduli[m], zetas[m]))
            {
                throw std::logic_error("Failed to find minimal primitive root for NTT negacyclic.");
            }
        }


        //std::cout << "\n2n-th primitive root for ntt negacyclic: " << zeta << "\n";

        // Negacyclic NTT
        // X[k] = ∑x[j]*ζ^(j(2k+1)) = ∑x[j]*ζ^j*(ζ^2)^jk, for all k,j ∈ {0, ... , n-1}
        // 
        // X[0] = x[0]*ζ^(0*1) + x[1]*ζ^(1*1) + x[2]*ζ^(2*1) ... + x[n-1]*ζ^((n-1)*1)
        // X[1] = x[0]*ζ^(0*3) + x[1]*ζ^(1*3) + x[2]*ζ^(2*3) ... + x[n-1]*ζ^((n-1)*3)
        // X[2] = x[0]*ζ^(0*5) + x[1]*ζ^(1*5) + x[2]*ζ^(2*5) ... + x[n-1]*ζ^((n-1)*5)
        // ...
        // X[n-1] = x[0]*ζ^(0*(2n-1)) + x[1]*ζ^(1*(2n-1)) + x[2]*ζ^(2*(2n-1)) ... + x[n-1]*ζ^((n-1)*(2n-1))
        // 
        // So, multiply ζ^j to ring's all coeff and use ζ²as Standard NTT's ω.
        std::vector<uint64_t> omegas(moduli.size());

        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            uint64_t zeta_pow = 1;

            for (uint64_t i = 0; i < poly_modulus_n.degree(); i++)
            {
                ring(m, i, mul_mod(ring(m, i), zeta_pow, moduli[m].value()));
                zeta_pow = mul_mod(zeta_pow, zetas[m], moduli[m].value());
            }

            omegas[m] = pow_mod(zetas[m], 2, moduli[m].value());
        }

        ntt(ring, omegas);
    }

    void inverse_ntt_negacyclic(RnsCycloRing& ring)
    {
        // Get ring data.
        const auto& poly_modulus_n = ring.poly_modulus();
        const auto& moduli = *ring.moduli();


        // Find 2n-th primitive root(ζ).
        std::vector<uint64_t> inv_zetas(moduli.size());
        const PolyModulus poly_modulus_2n(2 * poly_modulus_n.degree());

        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            if (!try_inverse_minimal_primitive_root(poly_modulus_2n, moduli[m], inv_zetas[m]))
            {
                throw std::logic_error("Failed to find minimal primitive root for NTT negacyclic.");
            }
        }


        // Negacyclic INTT
        // X[k] = 1/n * {∑x[j]*ζ^(-j(2k+1))} = 1/n * {∑x[j]*ζ^-j*(ζ^2)^jk}, for all k,j ∈ {0, ... , n-1}
        // 
        // X[0] = 1/n * {x[0]*ζ^-(0*1) + x[1]*ζ^-(1*1) + x[2]*ζ^-(2*1) ... + x[n-1]*ζ^-((n-1)*1)}
        // X[1] = 1/n * {x[0]*ζ^-(0*3) + x[1]*ζ^-(1*3) + x[2]*ζ^-(2*3) ... + x[n-1]*ζ^-((n-1)*3)}
        // X[2] = 1/n * {x[0]*ζ^-(0*5) + x[1]*ζ^-(1*5) + x[2]*ζ^-(2*5) ... + x[n-1]*ζ^-((n-1)*5)}
        // ...
        // x[n-1] = 1/n * {X[0]*ζ^-(0*(2n-1)) + X[1]*ζ^-(1*(2n-1)) + X[2]*ζ^-(2*(2n-1)) ... + X[n-1]*ζ^-((n-1)*(2n-1))}
        // 
        // So, use ζ²as Standard INTT's ω and multiply ζ^-j to ring's all coeff.
        std::vector<uint64_t> inv_omegas(moduli.size());

        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            inv_omegas[m] = pow_mod(inv_zetas[m], 2, moduli[m].value());
        }

        inverse_ntt(ring, inv_omegas);

        for (uint64_t m = 0; m < moduli.size(); m++)
        {
            uint64_t inv_zeta_pow = 1;

            for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
            {
                ring(m, j, mul_mod(ring(m, j), inv_zeta_pow, moduli[m].value()));
                inv_zeta_pow = mul_mod(inv_zeta_pow, inv_zetas[m], moduli[m].value());
            }
        }
    }
}