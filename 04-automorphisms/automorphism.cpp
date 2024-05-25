#include <cstdint>
#include <cstring>
#include <vector>

#include "automorphisms.h"
#include "spqlios/arithmetic/vec_znx_arithmetic.h"

void apply_automorphism_on_plaintext(const MODULE* module, int64_t p, uint64_t k,  //
                                     int64_t* res_mu,                              //
                                     const int64_t* in_mu) {
  uint8_t* tmp_space = get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module, message_limbs, message_limbs));
  // TODO: 1. Apply ZnX automorphism
  //       2. Base 2k normalize
}

void create_keyswitch(const MODULE* module, int64_t p, uint64_t k, VMP_PMAT* autom_ks_a, VMP_PMAT* autom_ks_b,
                      const int64_t* skey, const SVP_PPOL* skey_dft) {
  uint8_t* tmp_space = get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, autom_nrows, autom_ncols));
  std::vector<int64_t> autom_s(N);
  // TODO: Apply ZnX automorphism on skey and store in autom_s
  vec_znx_automorphism(module, p,             //
                       autom_s.data(), 1, N,  //
                       skey, 1, N);

  DECLARE_3D_ZERO_MATRIX(q_a, int64_t, autom_nrows, autom_ncols, N);
  DECLARE_3D_ZERO_MATRIX(q_b, int64_t, autom_nrows, autom_ncols, N);
  DECLARE_3D_ZERO_MATRIX(mu_b, int64_t, autom_nrows, autom_ncols, N);
  for (uint64_t j = 0; j < autom_nrows; ++j) {
    memcpy(mu_b[j][j], autom_s.data(), N * sizeof(int64_t));
  }
  for (uint64_t j = 0; j < autom_nrows; ++j) {
    // TODO: Apply rlwe_encrypt
  }

  // TODO: Apply s to q_a and q_b to create prepared matrices by using vmp_prepare_contiguous
}

void apply_automorphism(const MODULE* module, int64_t p, uint64_t k,  //
                        int64_t* res_a, int64_t* res_b,               //
                        const int64_t* in_a, const int64_t* in_b, const VMP_PMAT* autom_ks_a,
                        const VMP_PMAT* autom_ks_b) {
  std::vector<int64_t> autom_a(ell * N);
  std::vector<int64_t> autom_b(ell * N);
  // TODO: Apply ZnX automorphism to in_a and in_b

  VEC_ZNX_DFT* autom_a_dft = vec_znx_dft_alloc(module, autom_ncols);  // a * automorhism(s)
  VEC_ZNX_DFT* temp_dft = vec_znx_dft_alloc(module, autom_ncols);
  VEC_ZNX_BIG* temp_big = (VEC_ZNX_BIG*)temp_dft;  // alias
  uint8_t* tmp_space =
      get_tmp_space(std::max(vmp_apply_dft_to_dft_tmp_bytes(module, autom_ncols, ell, autom_nrows, autom_ncols),
                             vec_znx_big_normalize_base2k_tmp_bytes(module, ell, ell)));
  // TODO: Apply ZnX dft to autom_a, then apply dft_to_dft and store in temp_dft

  // TODO: Apply idft to temp_dft and store in temp_big


  // Next Step: create the output for b
  // TODO: 1. sub tmp_b to the b part of res (because of this, the final result is not normalized)
  //       2. Then normalized it in base2k and store in res_b

  // Next Step: create the output for a
  // TODO: 1. apply dft_to_dft to autom_a_dft
  //       2. Apply idft
  //       3. Negate the result by using vec_znx_big_sub with the first input NULL and size 0
  //       4. Finally, normalized it in base2k and store in res_a

  free(temp_dft);
}
