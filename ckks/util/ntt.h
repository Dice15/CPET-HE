#pragma once

#include "polymodulus.h"
#include "basis.h"
#include <cstdint>
#include <vector>


namespace cpet
{
    class NTT
    {
    public:
        NTT() = default;

        NTT(const PolyModulus& poly_modulus, const Basis& basis);

        void ntt_negacyclic(std::vector<std::vector<uint64_t>>& rings, uint64_t basis_begin, uint64_t basis_end) const;

        void inverse_ntt_negacyclic(std::vector<std::vector<uint64_t>>& vectors, uint64_t basis_begin, uint64_t basis_end) const;

    private:
        bool is_primitive_root(uint64_t n, uint64_t root, uint64_t prime) const;

        bool try_primitive_root(uint64_t n, uint64_t prime, uint64_t& destination) const;

        bool try_minimal_primitive_root(uint64_t n, uint64_t prime, uint64_t& destination) const;

        void ntt(std::vector<std::vector<uint64_t>>& rings, uint64_t basis_begin, uint64_t basis_end, bool inverse = false) const;

        void inverse_ntt(std::vector<std::vector<uint64_t>>& vectors, uint64_t basis_begin, uint64_t basis_end) const;

        uint64_t poly_modulus_degree_;

        std::vector<uint64_t> basis_;

        std::vector<uint64_t> bit_reversal_table_;

        std::vector<std::vector<uint64_t>> zeta_powers_by_basis_;

        std::vector<std::vector<uint64_t>> omega_powers_by_basis_;

        std::vector<std::vector<uint64_t>> inv_zeta_powers_by_basis_;

        std::vector<std::vector<uint64_t>> inv_omega_powers_by_basis_;
    };
}