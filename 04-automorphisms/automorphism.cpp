#include <cstdint>
#include <cstring>
#include <vector>

#include "automorphisms.h"
#include "spqlios/arithmetic/vec_znx_arithmetic.h"

void apply_automorphism_on_plaintext(const MODULE* module, int64_t p, uint64_t k,  //
                                     int64_t* res_mu,                              //
                                     const int64_t* in_mu) {
  uint8_t* tmp_space = get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module));
  vec_znx_automorphism(module, p,                 //
                       res_mu, message_limbs, N,  //
                       in_mu, message_limbs, N);
  vec_znx_normalize_base2k(module, K,
                           res_mu, message_limbs, N,
                           res_mu, message_limbs, N,
                           tmp_space);
}

void create_keyswitch(const MODULE* module, int64_t p, uint64_t k, VMP_PMAT* autom_ks_a, VMP_PMAT* autom_ks_b,
                      const int64_t* skey, const SVP_PPOL* skey_dft) {
  uint8_t* tmp_space = get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, autom_nrows, autom_ncols));
  std::vector<int64_t> autom_s(N);
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
    rlwe_encrypt(module, k, *q_a[j], *q_b[j], autom_ncols, *mu_b[j], j+1, skey_dft);
  }

  vmp_prepare_contiguous(module, autom_ks_a, **q_a, autom_nrows, autom_ncols, tmp_space);
  vmp_prepare_contiguous(module, autom_ks_b, **q_b, autom_nrows, autom_ncols, tmp_space);
  DECLARE_3D_ZERO_MATRIX(decr_b, int64_t, autom_nrows, autom_ncols, N);
  //for (uint64_t j = 0; j < autom_nrows; ++j) {
  //  rlwe_decrypt(module, k, *decr_b[j], *q_a[j], *q_b[j], autom_ncols, skey_dft);
  //}

}

void apply_automorphism(const MODULE* module, int64_t p, uint64_t k,  //
                        int64_t* res_a, int64_t* res_b,               //
                        const int64_t* in_a, const int64_t* in_b, const VMP_PMAT* autom_ks_a,
                        const VMP_PMAT* autom_ks_b) {
  // autom_a = apply public automorphism the a part of in_rlwe
  // autom_b = apply public automorphism the b part of in_rlwe
  std::vector<int64_t> autom_a(ell * N);
  std::vector<int64_t> autom_b(ell * N);
  vec_znx_automorphism(module, p,               //
                       autom_a.data(), ell, N,  //
                       in_a, ell, N);
  vec_znx_automorphism(module, p,               //
                       autom_b.data(), ell, N,  //
                       in_b, ell, N);

  // autom_a_s = vmp(tmp_a, autom_ks)  -- apply idft
  VEC_ZNX_DFT* autom_a_dft = vec_znx_dft_alloc(module, autom_ncols);  // a * automorhism(s)
  VEC_ZNX_DFT* temp_dft = vec_znx_dft_alloc(module, autom_ncols);
  VEC_ZNX_BIG* temp_big = (VEC_ZNX_BIG*)temp_dft;  // alias
  uint8_t* tmp_space = get_tmp_space(std::max(
                                         vmp_apply_dft_to_dft_tmp_bytes(module, autom_ncols, ell, autom_nrows, autom_ncols),
                                         vec_znx_big_normalize_base2k_tmp_bytes(module)));
  vec_znx_dft(module, //
              autom_a_dft, ell, //
              autom_a.data(), ell, N);
  vmp_apply_dft_to_dft(module, //
                       temp_dft, autom_ncols, //
                       autom_a_dft, ell, //
                       autom_ks_b, autom_nrows, autom_ncols,
                       tmp_space);

  vec_znx_idft(module,                 //
               temp_big, autom_ncols,  //
               temp_dft, autom_ncols,  //
               tmp_space);
  // res = big_normalize(autom_a_s)

  // sub tmp_b to the b part of res (because of this, the final result is not normalized)
  vec_znx_big_sub_small_a(module,                 // N
                          temp_big, autom_ncols,  // res
                          autom_b.data(), ell, N, //
                          temp_big, autom_ncols);

  vec_znx_big_normalize_base2k(module, k, //
                               res_b, ell, N,  //
                               temp_big, autom_ncols, //
                               tmp_space);

  vmp_apply_dft_to_dft(module, //
                       temp_dft, autom_ncols, //
                       autom_a_dft, ell, //
                       autom_ks_a, autom_nrows, autom_ncols, //
                       tmp_space);
  vec_znx_idft(module,                 //
               temp_big, autom_ncols,  //
               temp_dft, autom_ncols,  //
               tmp_space);

  vec_znx_big_sub(module,                 // N
                  temp_big, autom_ncols,  // res
                  NULL, 0, //
                  temp_big, autom_ncols);
  vec_znx_big_normalize_base2k(module, k, //
                               res_a, ell, N, //
                               temp_big, autom_ncols, //
                               tmp_space);

  free(temp_dft);
  free(autom_a_dft);
}
