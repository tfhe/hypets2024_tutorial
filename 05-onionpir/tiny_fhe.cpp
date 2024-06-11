#include "tiny_fhe.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>

#include "spqlios/arithmetic/vec_znx_arithmetic.h"

rng& randgen() {
  static thread_local rng gen;
  return gen;
}

void random_centered_reduced(uint64_t n, uint64_t log2_base2k, int64_t* res) {
  const int64_t bound = 1l << (log2_base2k - 1);
  std::uniform_int_distribution<int64_t> dist_base2k(-bound, bound - 1);
  rng& r = randgen();
  for (uint64_t i = 0; i < n; ++i) {
    res[i] = dist_base2k(r);
  }
}

void random_log2bound_symmetric(uint64_t n, uint64_t log2bound, int64_t* res) {
  const int64_t bound = 1l << log2bound;
  std::uniform_int_distribution<int64_t> dist(-bound, bound);
  rng& r = randgen();
  for (uint64_t i = 0; i < n; ++i) {
    res[i] = dist(r);
  }
}

void random_normal(uint64_t n, double stddev, int64_t* res) {
  std::normal_distribution<double> dist_normal(0, stddev);
  rng& r = randgen();
  for (uint64_t i = 0; i < n; ++i) {
    res[i] = (int64_t)round(dist_normal(r));
  }
}

void random_binary(uint64_t n, int64_t* res) {
  static std::uniform_int_distribution<int64_t> dist_binary(0, 1);
  rng& r = randgen();
  for (uint64_t i = 0; i < n; ++i) {
    res[i] = dist_binary(r);
  }
}

void round_polynomial_noise(uint64_t nn,              // ring dimension
                            int64_t k,                // message is k-normalized
                            int64_t message_basebit,  //
                            double* rem, int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a,
                            uint64_t a_size, uint64_t a_sl) {
  REQUIRE_DRAMATICALLY(k <= 62, "k is too large for this function");
  REQUIRE_DRAMATICALLY(message_basebit < k * int64_t(a_size), "message is too large: no rem");
  const int64_t _2pk = 1l << k;
  const int64_t _2pkm = 1l << (k - 1);
  int64_t jstart = message_basebit / k;
  uint64_t lrem = (k + ((k - message_basebit) % k)) % k;
  for (uint64_t i = 0; i < nn; ++i) {
    double rem_j = 0;
    int64_t carry_j = 0;
    int64_t j;
    // from a_size downto jstart excluded: the incoming carry_j on block j
    // is the amount we need to add to block j so that the contribution
    // of the remainder: sum(a_t 2^{-K(t+1)}) for t>j
    // is between [0 and 2^-Kt[
    // all a_j shall be between [0 and 2^K[
    for (j = a_size - 1; j > jstart; --j) {
      int64_t aj = a[j * a_sl + i] + carry_j;
      if (aj < 0) {
        aj += _2pk;  // makes it in [0,_2pk[
        carry_j = -1;
      } else {
        // aj = aj; is already in [0,_2pk[
        carry_j = 0;
      }
      rem_j += aj * pow(2., -k * (j + 1));
      if (j < int64_t(res_size)) res[i + j * res_sl] = 0;
    }
    // from jstart downto 0 included: the incoming carry_j on block j
    // is the amount we need to add to block j so that the contribution
    // of the remainder is between [2^-(Kt-1) and 2^-(Kt-1)[
    // all a_j shall be between [-2^(K-1) and 2^(K-1)[
    {
      REQUIRE_DRAMATICALLY(j == jstart, "bug");
      int64_t aj = a[j * a_sl + i] + carry_j;
      if (lrem == 0) {
        if (aj < -_2pkm) {
          aj += _2pk;
          carry_j = -1;
        } else {
          // aj = aj;
          carry_j = 0;
        }
        rem_j += aj * pow(2., -k * (j + 1));
        if (j < int64_t(res_size)) res[i + j * res_sl] = 0;
      } else {
        // eliminate/round the lrem lsbs
        int64_t aj_lo = (aj << (64 - lrem)) >> (64 - lrem);
        rem_j += aj_lo * pow(2., -k * (j + 1));
        aj -= aj_lo;
        // ensure it is still centered
        if (aj < -_2pkm) {
          aj += _2pk;
          carry_j = -1;
        } else if (aj >= _2pkm) {
          aj -= _2pk;
          carry_j = 1;
        } else {
          // aj = aj;
          carry_j = 0;
        }
        if (j < int64_t(res_size)) res[j * res_sl + i] = aj;
      }
      j--;
    }
    for (; j >= 0; j--) {
      int64_t aj = a[j * a_sl + i] + carry_j;
      if (aj < -_2pkm) {
        aj += _2pk;
        carry_j = -1;
      } else if (aj >= _2pkm) {
        aj -= _2pk;
        carry_j = 1;
      } else {
        // aj = aj;
        carry_j = 0;
      }
      if (j < int64_t(res_size)) res[j * res_sl + i] = aj;
    }
    if (rem != nullptr) rem[i] = rem_j;
  }
  // fill some zeros if the result is larger
  if (res_size > a_size) {
    for (uint64_t j = a_size; j < res_size; ++j) {
      memset(res + j * res_sl, 0, nn * sizeof(int64_t));
    }
  }
}

void torus_coeffs_values(uint64_t nn, uint64_t k, double* res, const int64_t* a, uint64_t a_size, uint64_t a_sl) {
  round_polynomial_noise(nn, k, 0, res, nullptr, 0, 0, a, a_size, a_sl);
}

void noise_values(uint64_t nn,               // ring dimension
                  uint64_t k,                // message is k-normalized
                  uint64_t message_basebit,  // message basis
                  double* res, int64_t* a, uint64_t a_size, uint64_t a_sl) {
  round_polynomial_noise(nn, k, message_basebit, res, nullptr, 0, 0, a, a_size, a_sl);
}

void round_polynomial(uint64_t nn,               // ring dimension
                      uint64_t k,                // message is k-normalized
                      uint64_t message_basebit,  //
                      int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a, uint64_t a_size,
                      uint64_t a_sl) {
  round_polynomial_noise(nn, k, message_basebit, nullptr, res, res_size, res_sl, a, a_size, a_sl);
}

void apply_automorphism(const MODULE* module, int64_t p, uint64_t k,  //
                        int64_t* res_rlwe, uint64_t res_size,         //
                        const int64_t* in_rlwe, uint64_t in_size,     //
                        const VMP_PMAT* autom_ks, uint64_t autom_key_nrows, uint64_t autom_key_ncols) {
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

void rlwe_encrypt_base2k(const MODULE* module,                                    //
                         uint64_t log2_base2k,                                    // output base 2^K
                         int64_t* new_b, uint64_t new_b_size, uint64_t new_b_sl,  // b part of rlwe
                         const SVP_PPOL* s,                                       // secret key
                         const int64_t* a, uint64_t a_size, uint64_t a_sl,        // a part of rlwe
                         const int64_t* phi, uint64_t phi_size, uint64_t phi_sl   // message + noise
) {
  assert(phi_size <= a_size);
  assert(phi_size <= new_b_size);

  const uint64_t res_size = a_size < new_b_size ? a_size : new_b_size;

  // dft(a.s)
  VEC_ZNX_DFT* res_dft = vec_znx_dft_alloc(module, res_size);
  svp_apply_dft(module, res_dft, res_size, s, a, res_size, a_sl);

  // allocate temporary space for idft and big normalize
  const uint64_t idft_bytes = vec_znx_idft_tmp_bytes(module);
  const uint64_t norm_bytes = vec_znx_big_normalize_base2k_tmp_bytes(module);
  uint8_t* tmp = get_tmp_space(idft_bytes > norm_bytes ? idft_bytes : norm_bytes);

  // a.s
  VEC_ZNX_BIG* res_big = (VEC_ZNX_BIG*)res_dft;  // inplace
  vec_znx_idft(module, res_big, res_size, res_dft, res_size, tmp);

  // a.s + phi
  vec_znx_big_add_small(module, res_big, res_size, res_big, res_size, phi, phi_size, phi_sl);

  // normalized a.s + phi
  vec_znx_big_normalize_base2k(module, log2_base2k, new_b, new_b_size, new_b_sl, res_big, res_size, tmp);

  free(res_big);
}

void rlwe_phase_base2k(const MODULE* module,                              //
                       uint64_t log2_base2k,                              // output base 2^K
                       int64_t* phi, uint64_t phi_size, uint64_t phi_sl,  // decrypted phase
                       const SVP_PPOL* s,                                 // secret key
                       const int64_t* a, uint64_t a_size, uint64_t a_sl,  // a part of rlwe
                       const int64_t* b, uint64_t b_size, uint64_t b_sl   // message + noise
) {
  const uint64_t work_size = a_size < b_size ? b_size : a_size;

  // dft(a.s)
  VEC_ZNX_DFT* res_dft = vec_znx_dft_alloc(module, work_size);
  svp_apply_dft(module,           //
                res_dft, a_size,  //
                s,                // s
                a, a_size, a_sl);

  // allocate temporary space for idft and big normalize
  const uint64_t idft_bytes = vec_znx_idft_tmp_bytes(module);
  const uint64_t norm_bytes = vec_znx_big_normalize_base2k_tmp_bytes(module);
  uint8_t* tmp = get_tmp_space(idft_bytes > norm_bytes ? idft_bytes : norm_bytes);

  // a.s
  VEC_ZNX_BIG* res_big = (VEC_ZNX_BIG*)res_dft;  // inplace
  vec_znx_idft(module,                           //
               res_big, a_size,                  //
               res_dft, a_size, tmp);

  // b - a.s
  vec_znx_big_sub_small_a(module, res_big, work_size,  //
                          b, b_size, b_sl,             //
                          res_big, a_size);

  // normalized b - a.s
  vec_znx_big_normalize_base2k(module, log2_base2k,    //
                               phi, phi_size, phi_sl,  //
                               res_big, work_size, tmp);

  free(res_big);
}

void rlwe_trace_expand_1step(const MODULE* module,                      //
                             int64_t autom_p, int64_t rotate_p,         //
                             uint64_t log2_base2k,                      //
                             int64_t* res0_rlwe, uint64_t res0_size,    //
                             int64_t* res1_rlwe, uint64_t res1_size,    //
                             const int64_t* in_rlwe, uint64_t in_size,  //
                             const VMP_PMAT* autom_ks, uint64_t autom_key_nrows, uint64_t autom_key_ncols) {
  const uint64_t nn = module_get_n(module);

  const uint64_t tmp_size = std::max(in_size, autom_key_ncols);  // limitation of apply_automorphism
  const uint64_t a_tmp_size = (tmp_size + 1) / 2;
  const uint64_t b_tmp_size = tmp_size / 2;

  std::vector<int64_t> tmp_raw(2 * nn * tmp_size);
  int64_t* tmp = tmp_raw.data();
  int64_t* tmp1 = tmp_raw.data() + nn * tmp_size;

  // res0 = in(X^p)
  apply_automorphism(module, autom_p, log2_base2k, tmp, tmp_size, in_rlwe, in_size, autom_ks, autom_key_nrows,
                     autom_key_ncols);

  // res1 = rotate(in - in(X^p))
  vec_znx_sub(module, tmp1, tmp_size, nn, in_rlwe, in_size, nn, tmp, tmp_size, nn);
  vec_znx_rotate(module, rotate_p, tmp1, tmp_size, nn, tmp1, tmp_size, nn);

  // normalize res1
  const uint64_t a_res1_size = (res1_size + 1) / 2;
  const uint64_t b_res1_size = res1_size / 2;

  vec_znx_normalize_base2k(module, log2_base2k, res1_rlwe, a_res1_size, 2 * nn, tmp1, a_tmp_size, 2 * nn,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module)));

  vec_znx_normalize_base2k(module, log2_base2k, res1_rlwe + nn, b_res1_size, 2 * nn, tmp1 + nn, b_tmp_size, 2 * nn,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module)));

  // res0 = in + in(X^p)
  vec_znx_add(module, tmp1, tmp_size, nn, in_rlwe, in_size, nn, tmp, tmp_size, nn);

  // normalize res0
  const uint64_t a_res0_size = (res0_size + 1) / 2;
  const uint64_t b_res0_size = res0_size / 2;

  vec_znx_normalize_base2k(module, log2_base2k, res0_rlwe, a_res0_size, 2 * nn, tmp1, a_tmp_size, 2 * nn,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module)));

  vec_znx_normalize_base2k(module, log2_base2k, res0_rlwe + nn, b_res0_size, 2 * nn, tmp1 + nn, b_tmp_size, 2 * nn,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module)));
}

void rlwe_trace_expand_1step_no_rotation(const MODULE* module,                      //
                                         int64_t autom_p,                           //
                                         uint64_t log2_base2k,                      //
                                         int64_t* res0_rlwe, uint64_t res0_size,    //
                                         const int64_t* in_rlwe, uint64_t in_size,  //
                                         const VMP_PMAT* autom_ks, uint64_t autom_key_nrows, uint64_t autom_key_ncols) {
  const uint64_t nn = module_get_n(module);

  const uint64_t tmp_size = std::max(in_size, autom_key_ncols);  // limitation of apply_automorphism
  std::vector<int64_t> tmp_raw(nn * tmp_size);
  int64_t* tmp = tmp_raw.data();

  // res0 = in(X^p)
  apply_automorphism(module, autom_p, log2_base2k, tmp, tmp_size, in_rlwe, in_size, autom_ks, autom_key_nrows,
                     autom_key_ncols);

  // res0 = in + in(X^p)
  vec_znx_add(module, tmp, tmp_size, nn, in_rlwe, in_size, nn, tmp, tmp_size, nn);

  // normalize res0
  const uint64_t a_tmp_size = (tmp_size + 1) / 2;
  const uint64_t b_tmp_size = tmp_size / 2;
  const uint64_t a_res0_size = (res0_size + 1) / 2;
  const uint64_t b_res0_size = res0_size / 2;

  vec_znx_normalize_base2k(module, log2_base2k, res0_rlwe, a_res0_size, 2 * nn, tmp, a_tmp_size, 2 * nn,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module)));

  vec_znx_normalize_base2k(module, log2_base2k, res0_rlwe + nn, b_res0_size, 2 * nn, tmp + nn, b_tmp_size, 2 * nn,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module)));
}
