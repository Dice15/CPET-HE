
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

    void ntt(CycloModRing& ring, uint64_t omega)
    {
        // Get ring data.
        const PolyModulus& poly_modulus = ring.poly_modulus();
        const Modulus& modulus = ring.modulus();


        // Rearrange ring's coeff using ntt table.
        std::vector<uint64_t> ntt_table;
        create_ntt_table(poly_modulus, ntt_table);

        for (uint64_t i = 1; i < poly_modulus.degree(); i++)
        {
            uint64_t j = ntt_table[i];

            if (i < j)
            {
                uint64_t coeff_i = ring(i);
                ring(i, ring(j));
                ring(j, coeff_i);
            }
        }

       
        // X[k] = ∑x[j]*ω^(jk), for all k,j ∈ {0, ... , n-1}
        // We will use Cooley-Tukey algorithm for NTT.
        for (uint64_t len = 2; len <= poly_modulus.degree(); len <<= 1) 
        {
            uint64_t wlen = pow_mod(omega, poly_modulus.degree() / len, modulus.value());

            for (size_t i = 0; i < poly_modulus.degree(); i += len)
            {
                uint64_t w = 1;

                for (size_t j = 0; j < (len >> 1); j++) 
                {
                    uint64_t u = ring(i + j);
                    uint64_t v = mul_mod(ring(i + j + (len >> 1)), w, modulus.value());

                    ring(i + j, add_mod(u, v, modulus.value()));
                    ring(i + j + (len >> 1), sub_mod(u, v, modulus.value()));

                    w = mul_mod(w, wlen, modulus.value());
                }
            }
        }
    }

    void inverse_ntt(CycloModRing& ring, uint64_t omega)
    {
        // Get ring data.
        const PolyModulus& poly_modulus = ring.poly_modulus();
        const Modulus& modulus = ring.modulus();


        // X[k] = 1/n * {∑x[j]*ω^(-jk)}, for all k,j ∈ {0, ... , n-1}
        // So, use ω^-1 as Standard NTT's ω for INTT.
        uint64_t inv_omega = inverse_mod(omega, modulus.value());
        ntt(ring, inv_omega);


        // Scaling coeff to 1/n * coeff
        uint64_t n_inv = inverse_mod(poly_modulus.degree(), modulus.value());

        for (size_t i = 0; i < poly_modulus.degree(); i++)
        {
            ring(i, mul_mod(ring(i), n_inv, modulus.value()));
        }
    }

    void ntt_negacyclic(CycloModRing& ring)
    {
        // Get ring data.
        const PolyModulus& poly_modulus_n = ring.poly_modulus();
        const Modulus& modulus = ring.modulus();


        // Find 2n-th primitive root(ζ). 
        uint64_t zeta;
        const PolyModulus poly_modulus_2n(2 * poly_modulus_n.degree());

        if (!try_minimal_primitive_root(poly_modulus_2n, ring.modulus(), zeta))
        {
            throw std::logic_error("Failed to find minimal primitive root for NTT negacyclic.");
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
        uint64_t zeta_pow = 1;

        for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
        {
            ring(j, mul_mod(ring(j), zeta_pow, modulus.value()));
            zeta_pow = mul_mod(zeta_pow, zeta, modulus.value());
        }

        uint64_t omega = pow_mod(zeta, 2, modulus.value());
        ntt(ring, omega);
    }

    void inverse_ntt_negacyclic(CycloModRing& ring)
    {
        // Get ring data.
        const PolyModulus& poly_modulus_n = ring.poly_modulus();
        const Modulus& modulus = ring.modulus();


        // Find 2n-th primitive root(ζ).
        uint64_t zeta;
        const PolyModulus poly_modulus_2n(2 * poly_modulus_n.degree());

        if (!try_minimal_primitive_root(poly_modulus_2n, modulus, zeta))
        {
            throw std::logic_error("Failed to find minimal primitive root for inverse NTT negacyclic.");
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
        uint64_t omega = pow_mod(zeta, 2, modulus.value());
        inverse_ntt(ring, omega);

        uint64_t inv_zeta = inverse_mod(zeta, modulus.value());
        uint64_t inv_zeta_pow = 1;

        for (uint64_t j = 0; j < poly_modulus_n.degree(); j++)
        {
            ring(j, mul_mod(ring(j), inv_zeta_pow, modulus.value()));
            inv_zeta_pow = mul_mod(inv_zeta_pow, inv_zeta, modulus.value());
        }
    }
}