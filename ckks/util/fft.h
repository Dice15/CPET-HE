#pragma once

#include <cstdint>
#include <vector>
#include "cycloring.h"
#include "polymodulus.h"


namespace cpet
{
    void compute_minimal_primitive_root_powers(const PolyModulus& poly_modulus, std::vector<std::complex<double_t>>& destination);

    void create_fft_table(const PolyModulus& poly_modulus, std::vector<uint64_t>& destination);

    void fft(CycloRing& ring, const std::vector<std::complex<double_t>>& omega_powers);

    void inverse_fft(CycloRing& ring, const std::vector<std::complex<double_t>>& inv_omega_powers);

    void fft_negacyclic(CycloRing& ring);

    void inverse_fft_negacyclic(CycloRing& ring);

    void inverse_variant_canonical_embedding(CycloRing& ring);

    void variant_canonical_embedding(CycloRing& ring);
}