#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "cmath"
#include "gtest/gtest.h"
#include "onionpir.h"
#include "tiny_fhe.h"
#include "spqlios/arithmetic/vec_znx_arithmetic.h"

/** @file this test file mimics what a user would actually write: only public api */

TEST(hello, hello) { ASSERT_EQ(0, 0); }

TEST(onionpir, encode_decode) {
  MODULE_TYPE mtype = rand() % 2 == 0 ? FFT64 : NTT120;
  MODULE* mod = new_module_info(ONIONPIR_N, mtype);

  // generate random data
  uint8_t* inp = new uint8_t[7 * ONIONPIR_N];
  for (uint64_t i = 0; i < 7 * ONIONPIR_N; ++i) {
    inp[i] = rand() & 255;
  }

  int64_t* encoded = new int64_t[4 * ONIONPIR_N];
  onionpir_encode_plaintext(mod, encoded, inp);

  // check encoded values do not overflow
  for (uint64_t i = 0; i < 4 * ONIONPIR_N; ++i) {
    ASSERT_LT(encoded[i], 1l << 15);
    ASSERT_GE(encoded[i], -(1l << 15));
  }
  // additionally last limb must contain only multiples of 2^8
  for (uint64_t i = 3 * ONIONPIR_N; i < 4 * ONIONPIR_N; ++i) {
    ASSERT_EQ(encoded[i] & 0xFF, 0);
  }

  uint8_t* decoded = new uint8_t[7 * ONIONPIR_N];
  onionpir_decode_plaintext(mod, decoded, encoded);

  // check decode(encode(inp)) == inp
  for (uint64_t i = 0; i < 7 * ONIONPIR_N; ++i) {
    ASSERT_EQ(inp[i], decoded[i]) << i;
  }

  delete[] decoded;
  delete[] encoded;
  delete[] inp;
  delete_module_info(mod);
}

TEST(tiny_fhe, random_centered_reduced) {
  for (const uint64_t n : {1ul << 10, 1ul << 16}) {
    for (const uint64_t log2_base2k : {5, 19}) {
      int64_t* vec = new int64_t[n];
      random_centered_reduced(n, log2_base2k, vec);

      const int64_t bound = 1l << (log2_base2k - 1);
      double mean = 0;
      for (uint64_t i = 0; i < n; ++i) {
        ASSERT_GE(vec[i], -bound);
        ASSERT_LT(vec[i], bound);
        mean += vec[i];
      }
      mean /= n;

      // empirically chosen so that 1000 test executions pass
      ASSERT_GT(mean, -0.1 * bound);
      ASSERT_LT(mean, 0.1 * bound);

      delete[] vec;
    }
  }
}

TEST(tiny_fhe, random_log2bound_symmetric) {
  for (const uint64_t n : {1ul << 10, 1ul << 16}) {
    for (const uint64_t log2bound : {5, 19}) {
      int64_t* vec = new int64_t[n];
      random_log2bound_symmetric(n, log2bound, vec);

      const int64_t bound = 1l << log2bound;
      double mean = 0;
      for (uint64_t i = 0; i < n; ++i) {
        ASSERT_GE(vec[i], -bound);
        ASSERT_LE(vec[i], bound);
        mean += vec[i];
      }
      mean /= n;

      // empirically chosen so that 1000 test executions pass
      ASSERT_GT(mean, -0.1 * bound);
      ASSERT_LT(mean, 0.1 * bound);

      delete[] vec;
    }
  }
}

TEST(tiny_fhe, random_normal) {
  for (const uint64_t n : {1ul << 10, 1ul << 16}) {
    for (const uint64_t log2bound : {5, 19}) {
      const double scale = (double)(1l << log2bound) / 7.14355203435219;  // Prob(-bound <= X <= bound) < 2^-40

      int64_t* vec = new int64_t[n];
      random_normal(n, scale, vec);

      const int64_t bound = 1l << log2bound;
      double mean = 0;
      for (uint64_t i = 0; i < n; ++i) {
        ASSERT_GE(vec[i], -bound);
        ASSERT_LE(vec[i], bound);
        mean += vec[i];
      }
      mean /= n;

      // empirically chosen so that 1000 test executions pass
      ASSERT_GT(mean, -0.1 * bound) << bound;
      ASSERT_LT(mean, 0.1 * bound) << bound;

      delete[] vec;
    }
  }
}

TEST(tiny_fhe, random_binary) {
  for (const uint64_t n : {1ul << 10, 1ul << 16}) {
    int64_t* vec = new int64_t[n];
    random_binary(n, vec);

    double mean = 0;
    for (uint64_t i = 0; i < n; ++i) {
      ASSERT_GE(vec[i], 0);
      ASSERT_LE(vec[i], 1);
      mean += vec[i];
    }
    mean /= n;

    // empirically chosen so that 1000 test executions pass
    ASSERT_GE(mean, 0.44);
    ASSERT_LE(mean, 0.56);

    delete[] vec;
  }
}

TEST(tiny_fhe, rlwe_encrypt_decrypt) {
  const uint64_t nn = ONIONPIR_N;
  MODULE* mod = new_module_info(nn, FFT64);

  const uint64_t log2_base2k = ONIONPIR_K;
  const int64_t base2k = 1l << log2_base2k;

  for (const uint64_t a_size : {1, 2, 5}) {
    for (const uint64_t b_size : {a_size}) {
      // deactivated other dimensions, this test is not accurate
      for (uint64_t phi_size = 1; phi_size <= std::min(a_size, b_size); ++phi_size) {
        // generate random key
        int64_t* sk = new int64_t[nn];
        random_log2bound_symmetric(nn, 0, sk);

        // key to dft domain
        SVP_PPOL* s = new_svp_ppol(mod);
        svp_prepare(mod, s, sk);

        delete[] sk;

        // generate random message
        int64_t* phi = new int64_t[phi_size * nn];
        for (uint64_t i = 0; i < phi_size * nn; i += nn) {
          random_centered_reduced(nn, ONIONPIR_K, phi + i);
        }

        int64_t* a = new int64_t[a_size * nn];
        for (uint64_t i = 0; i < a_size * nn; i += nn) {
          random_centered_reduced(nn, ONIONPIR_K, a + i);
        }

        int64_t* b = new int64_t[b_size * nn];

        // encrypt
        rlwe_encrypt_base2k(mod, log2_base2k, b, b_size, nn, s, a, a_size, nn, phi, phi_size, nn);

        // check b is base2k reduced
        for (uint64_t i = 0; i < b_size * nn; ++i) {
          ASSERT_GE(b[i], -base2k / 2);
          ASSERT_LT(b[i], base2k / 2);
        }

        // decrypt
        int64_t* phi_res = new int64_t[phi_size * nn];
        rlwe_phase_base2k(mod, log2_base2k, phi_res, phi_size, nn, s, a, a_size, nn, b, b_size, nn);

        delete[] b;
        delete[] a;

        // check equality
        for (uint64_t i = 0; i < phi_size * nn; ++i) {
          ASSERT_EQ(phi[i], phi_res[i]);
        }

        delete[] phi_res;
        delete[] phi;

        delete_svp_ppol(s);
      }
    }
  }

  delete_module_info(mod);
}

TEST(tiny_fhe, rlwe_encrypt_decrypt_1) {
  const uint64_t nn = ONIONPIR_N;
  MODULE* mod = new_module_info(nn, FFT64);

  const uint64_t log2_base2k = ONIONPIR_K;
  const int64_t base2k = 1l << log2_base2k;

  for (const uint64_t ab_size : {1, 2, 5}) {
    for (uint64_t phi_size = 1; phi_size <= ab_size; ++phi_size) {
      // generate random key
      int64_t* sk = new int64_t[nn];
      random_log2bound_symmetric(nn, 0, sk);

      // key to dft domain
      SVP_PPOL* s = new_svp_ppol(mod);
      svp_prepare(mod, s, sk);

      delete[] sk;

      // generate random message
      int64_t* phi = new int64_t[phi_size * nn];
      for (uint64_t i = 0; i < phi_size * nn; i += nn) {
        random_centered_reduced(nn, ONIONPIR_K, phi + i);
      }

      // allocate a and b elements on same memory region: a0, b0, a1, b1, ..
      int64_t* ab = new int64_t[2 * ab_size * nn];

      // fill a with random
      for (uint64_t i = 0; i < 2 * ab_size * nn; i += 2 * nn) {
        random_centered_reduced(nn, ONIONPIR_K, ab + i);
      }

      // encrypt
      rlwe_encrypt_base2k(mod, log2_base2k, ab + nn, ab_size, 2 * nn, s, ab, ab_size, 2 * nn, phi, phi_size, nn);

      // check b is base2k reduced
      for (uint64_t i = nn; i < 2 * ab_size * nn; i += 2 * nn) {
        for (uint64_t j = 0; j < nn; ++j) {
          ASSERT_GE(ab[i + j], -base2k / 2);
          ASSERT_LT(ab[i + j], base2k / 2);
        }
      }

      // decrypt
      int64_t* phi_res = new int64_t[phi_size * nn];
      rlwe_phase_base2k(mod, log2_base2k, phi_res, phi_size, nn, s, ab, ab_size, 2 * nn, ab + nn, ab_size, 2 * nn);

      delete[] ab;

      // check equality
      for (uint64_t i = 0; i < phi_size * nn; ++i) {
        ASSERT_EQ(phi[i], phi_res[i]);
      }

      delete[] phi_res;
      delete[] phi;

      delete_svp_ppol(s);
    }
  }

  delete_module_info(mod);
}
