#pragma once

#include <cstdint>
#include <vector>
#include "cycloring.h"
#include "polymodulus.h"


namespace cpet
{
    void compute_minimal_primitive_root_powers(const PolyModulus& poly_modulus, std::vector<std::complex<double_t>>& destination);

    void create_fft_table(const PolyModulus& poly_modulus, std::vector<uint64_t>& destination);

    void inverse_variant_canonical_embedding(CycloRing& ring);

    void variant_canonical_embedding(CycloRing& ring);
}