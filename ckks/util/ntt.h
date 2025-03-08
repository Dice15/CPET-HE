#pragma once

#include <cstdint>
#include <vector>
#include "cyclomodring.h"
#include "polymodulus.h"
#include "modulus.h"


namespace cpet
{
    bool is_primitive_root(uint64_t root, const PolyModulus& poly_modulus, const Modulus& modulus);

    bool try_primitive_root(const PolyModulus& poly_modulus, const Modulus& modulus, uint64_t& destination);

    bool try_minimal_primitive_root(const PolyModulus& poly_modulus, const Modulus& modulus, uint64_t& destination);

    void create_ntt_table(const PolyModulus& poly_modulus, std::vector<uint64_t>& destination);

    void ntt_negacyclic(CycloModRing& ring);

    void inverse_ntt_negacyclic(CycloModRing& ring);
}