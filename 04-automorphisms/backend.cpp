#include <cstring>
#include <vector>

#include "automorphisms.h"

uint8_t* get_tmp_space(uint64_t bytes);

void rlwe_encrypt(const MODULE* module, uint64_t k,         //
                  int64_t* a, int64_t* b, uint64_t nlimbs,  //
                  const int64_t* mu, const SVP_PPOL* skey) {
  std::vector<Int64VecN> tmp(nlimbs);
  VEC_ZNX_DFT* tmp_dft = vec_znx_dft_alloc(module, nlimbs);
  VEC_ZNX_BIG* tmp_big = (VEC_ZNX_BIG*)tmp_dft;
  uint64_t tmp_bytes = vec_znx_big_normalize_base2k_tmp_bytes(module, nlimbs, nlimbs);

  // randomize a completely
  for (uint64_t i = 0; i < nlimbs; ++i) {
    random_centered_reduced(N, K, a + i * N);
  }
  // generate a small noise in the last limb of mu
  memset(tmp.data(), 0, N * nlimbs * sizeof(int64_t));
  random_log2bound_symmetric(N, 0, tmp[nlimbs - 1]);
  vec_znx_add(module,                  //
              *tmp.data(), nlimbs, N,  //
              *tmp.data(), nlimbs, N,  //
              mu, message_limbs, N);
  // compute as + mu + e
  svp_apply_dft(module,           //
                tmp_dft, nlimbs,  //
                skey,             //
                a, nlimbs, N);
  vec_znx_idft_tmp_a(module,           //
                     tmp_big, nlimbs,  //
                     tmp_dft, nlimbs);
  vec_znx_big_add_small(module,           //
                        tmp_big, nlimbs,  //
                        tmp_big, nlimbs,  //
                        *tmp.data(), nlimbs, N);
  vec_znx_big_normalize_base2k(module, k,        //
                               b, nlimbs, N,     //
                               tmp_big, nlimbs,  //
                               get_tmp_space(tmp_bytes));
  free(tmp_dft);
}

double rlwe_decrypt(const MODULE* module, uint64_t k,                     //
                    int64_t* mu,                                          //
                    const int64_t* a, const int64_t* b, uint64_t nlimbs,  //
                    const SVP_PPOL* skey) {
  std::vector<Int64VecN> tmp(nlimbs);
  std::vector<double> rem(N);
  VEC_ZNX_DFT* tmp_dft = vec_znx_dft_alloc(module, nlimbs);
  VEC_ZNX_BIG* tmp_big = (VEC_ZNX_BIG*)tmp_dft;
  uint64_t tmp_bytes = vec_znx_big_normalize_base2k_tmp_bytes(module, nlimbs, nlimbs);

  // compute b - a.s
  svp_apply_dft(module,           //
                tmp_dft, nlimbs,  //
                skey,             //
                a, nlimbs, N);
  vec_znx_idft_tmp_a(module,           //
                     tmp_big, nlimbs,  //
                     tmp_dft, nlimbs);
  vec_znx_big_sub_small_a(module,           //
                          tmp_big, nlimbs,  //
                          b, nlimbs, N,     //
                          tmp_big, nlimbs);
  // we could normalize directly to the final output, but here,
  // let's dig deeper and print also some noise amplitude
  vec_znx_big_normalize_base2k(module, k,               //
                               *tmp.data(), nlimbs, N,  //
                               tmp_big, nlimbs,         //
                               get_tmp_space(tmp_bytes));
  round_polynomial_noise(N, k, k * message_limbs,  //
                         rem.data(),               //
                         mu, message_limbs, N,     //
                         *tmp.data(), nlimbs, N);
  double noise_ampl = 0;
  for (uint64_t i = 0; i < N; ++i) {
    double d = fabs(rem[i]);
    if (d > noise_ampl) noise_ampl = d;
  }
  std::cerr << "noise log2: " << log2(noise_ampl) << std::endl;
  free(tmp_dft);
  return log2(noise_ampl);
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

/** create and get some reusable tmp space */
uint8_t* get_tmp_space(const uint64_t bytes) {
  static __thread uint64_t size = 0;  // current tmp size
  static __thread uint8_t* space = nullptr;
  if (bytes > size) {
    free(space);
    if (size == 0) size = 1L << 20;  // minimal size
    // double the tmp size until it fits
    while (size < bytes) size <<= 1;
    // realloc the new size (free + alloc is enough for tmp space)
    space = (uint8_t*)aligned_alloc(64, size);
    REQUIRE_DRAMATICALLY(space != nullptr, "Out of memory");
  }
  return space;
}

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
