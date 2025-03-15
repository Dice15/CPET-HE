#pragma once

#include <cstdint>
#include <vector>
#include "rnscycloring.h"
#include "polymodulus.h"
#include "modulus.h"


namespace cpet
{
    bool is_primitive_root(uint64_t root, const PolyModulus& poly_modulus, uint64_t modulus);

    bool try_primitive_root(const PolyModulus& poly_modulus, uint64_t modulus, uint64_t& destination);

    bool try_minimal_primitive_root(const PolyModulus& poly_modulus, uint64_t modulus, uint64_t& destination);

    bool try_inverse_minimal_primitive_root(const PolyModulus& poly_modulus, uint64_t modulus, uint64_t& destination);

    void create_ntt_table(const PolyModulus& poly_modulus, std::vector<uint64_t>& destination);

    void ntt(RnsCycloRing& ring, const std::vector<uint64_t>& omegas);

    void inverse_ntt(RnsCycloRing& ring, const std::vector<uint64_t>& inv_omegas);

    void ntt_negacyclic(RnsCycloRing& ring);

    void inverse_ntt_negacyclic(RnsCycloRing& ring);
}