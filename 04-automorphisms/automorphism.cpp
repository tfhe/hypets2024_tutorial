#include <cstdint>

#include "spqlios/arithmetic/vec_znx_arithmetic.h"

constexpr uint64_t N = 64;  // fixed ring dimension
constexpr uint64_t K = 16;  // fixed limb size (2^16 bits)
constexpr uint64_t message_limbs = 4; // number of message limbs
constexpr uint64_t ell = 5; // size of input and output ringlwe
constexpr uint64_t autom_nrows = 5; // nrows of the autom matrix
constexpr uint64_t autom_ncols = 6; // nrows of the autom matrix

void rlwe_encrypt(const MODULE* module, uint64_t k, //
                  int64_t* a, int64_t* b, uint64_t nlimbs, //
                  const int64_t* mu);

double rlwe_decrypt(const MODULE* module, uint64_t k, //
                    int64_t* mu,
                    int64_t* a, int64_t* b, uint64_t nlimbs);

void apply_automorphism_on_plaintext(const MODULE* module, int64_t p, uint64_t k,  //
                        int64_t* res_mu,          //
                        const int64_t* in_mu);

void create_keyswitch(const MODULE* module, int64_t p, uint64_t k,  //
                      VMP_PMAT* autom_ks_a, VMP_PMAT* autom_ks_b);


void apply_automorphism(const MODULE* module, int64_t p, uint64_t k,  //
                        int64_t* res_a, int64_t* res_b,          //
                        const int64_t* in_a, const int64_t* in_b,
                        const VMP_PMAT* autom_ks_a, const VMP_PMAT* autom_ks_b) {
  const uint64_t nn = module_get_n(module);
  REQUIRE_DRAMATICALLY(res_size == std::max(in_size, autom_key_ncols), "bug: general res size not supported!");
  // largest size of in_a that can be taken into account
  const uint64_t a_size = std::min((in_size + 1) / 2, autom_key_nrows);
  // largest size of in_b that can be taken into account
  const uint64_t b_size = in_size / 2;
  const uint64_t a_res_size = (res_size + 1) / 2;
  const uint64_t b_res_size = res_size / 2;

  // autom_a = apply public automorphism the a part of in_rlwe
  // autom_b = apply public automorphism the b part of in_rlwe
  std::vector<int64_t> autom_a(a_size * nn);
  std::vector<int64_t> autom_b(b_size * nn);
  vec_znx_automorphism(module, p,                   //
                       autom_a.data(), a_size, nn,  //
                       in_rlwe, a_size, 2 * nn);
  vec_znx_automorphism(module, p,                   //
                       autom_b.data(), b_size, nn,  //
                       in_rlwe + nn, b_size, 2 * nn);

  // autom_a_s = vmp(tmp_a, autom_ks)  -- apply idft
  VEC_ZNX_DFT* autom_a_s_dft = vec_znx_dft_alloc(module, autom_key_ncols);  // a * automorhism(s)
  VEC_ZNX_BIG* autom_a_s_big = (VEC_ZNX_BIG*)autom_a_s_dft;                 // alias
  uint8_t* tmp_space =
      get_tmp_space(vmp_apply_dft_tmp_bytes(module, autom_key_ncols, a_size, autom_key_nrows, autom_key_ncols));
  vmp_apply_dft(module,                                      //
                autom_a_s_dft, autom_key_ncols,              //
                autom_a.data(), a_size, nn,                  //
                autom_ks, autom_key_nrows, autom_key_ncols,  //
                tmp_space);
  vec_znx_idft(module,                          //
               autom_a_s_big, autom_key_ncols,  //
               autom_a_s_dft, autom_key_ncols,  //
               tmp_space);
  // res = big_normalize(autom_a_s)
  vec_znx_big_range_normalize_base2k(module, k,                             //
                                     res_rlwe, a_res_size, 2 * nn,          // a part of res (even pos)
                                     autom_a_s_big, 0, autom_key_ncols, 2,  // a part (even range)
                                     tmp_space);
  vec_znx_big_range_normalize_base2k(module, k,                             //
                                     res_rlwe + nn, b_res_size, 2 * nn,     // b part of res (odd pos)
                                     autom_a_s_big, 1, autom_key_ncols, 2,  // b part (odd range)
                                     tmp_space);
  // add tmp_b to the b part of res (because of this, the final result is not normalized)
  vec_znx_add(module, res_rlwe + nn, b_res_size, 2 * nn, res_rlwe + nn, b_res_size, 2 * nn, autom_b.data(), b_size, nn);

  free(autom_a_s_dft);
}
