#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <cstdio>

#include "cmath"
#include "onionpir.h"
#include "test/testlib/fft64_layouts.h"
#include "test/testlib/polynomial_vector.h"
#include "tiny_fhe.h"

/** this test file can use private functionalities, private api, and the testlib */

/** evaluates one coefficient, centermod 2^-log2pow_mod */
static double simple_eval(int64_t log2pow_mod, int64_t k, const int64_t* a, uint64_t a_size, uint64_t a_sl) {
  double value = 0;
  for (int64_t j = a_size - 1; j >= 0; j--) {
    if (k * (j + 1) <= log2pow_mod) break;
    value += a[j * a_sl] * pow(2., -k * (j + 1));
  }
  double quotient = pow(2., log2pow_mod);
  return value - floor(value * quotient + 0.5) / quotient;
}

TEST(noise, round_polynomial_noise) {
  uint64_t nn = 512;
  MODULE* module = new_module_info(nn, FFT64);
  for (int64_t k : {8, 19}) {
    for (int64_t L : {8, 13, 19, 21}) {
      uint64_t sa = 6;
      uint64_t a_sl = nn + 3;
      // generate a base2k random reduced a
      std::vector<uint8_t> tmp_space(vec_znx_normalize_base2k_tmp_bytes(module));
      znx_vec_i64_layout a(nn, sa, a_sl);
      a.fill_random(50);
      vec_znx_normalize_base2k(module, k, a.data(), sa, a_sl, a.data(), sa, a_sl, tmp_space.data());
      a.data()[0] = -1l << (k - 1);  // one nasty corner case
      a.data()[a_sl] = -1l << (k - 1);
      const thash hash_input = a.content_hash();

      // compare simple-eval
      std::vector<double> evals(nn);
      torus_coeffs_values(nn, k, evals.data(), a.data(), sa, a_sl);
      ASSERT_EQ(a.content_hash(), hash_input);
      for (uint64_t i = 0; i < nn; ++i) {
        double expect = simple_eval(0, k, a.data() + i, sa, a_sl);
        double actual = evals[i];
        double diff = expect - evals[i];
        ASSERT_GE(actual, -0.5);
        ASSERT_LT(actual, 0.5);
        ASSERT_LE(fabs(diff), 1e-10);
      }
      uint64_t sr = sa - 1;
      uint64_t r_sl = 2 * nn;
      znx_vec_i64_layout res(nn, sr, r_sl);
      round_polynomial_noise(nn, k, L,
                             evals.data(),          // rem
                             res.data(), sr, r_sl,  // round
                             a.data(), sa, a_sl);
      ASSERT_EQ(a.content_hash(), hash_input);
      for (uint64_t i = 0; i < nn; ++i) {
        double expect = simple_eval(L, k, a.data() + i, sa, a_sl);
        double actual = evals[i];
        double diff = expect - actual;
        double quotient = pow(2., -L);
        ASSERT_GE(actual, -quotient / 2);
        ASSERT_LT(actual, quotient / 2);
        ASSERT_LE(fabs(diff), quotient * 1e-6);
        // check that res is reduced
        for (uint64_t j = 0; j < sr; ++j) {
          ASSERT_GE(res.data()[i + j * r_sl], -1l << (k - 1));
          ASSERT_LT(res.data()[i + j * r_sl], 1l << (k - 1));
        }
        double expect_full = simple_eval(0, k, a.data() + i, sa, a_sl);
        double expect_res = simple_eval(0, k, res.data() + i, sr, r_sl);
        diff = expect_full - expect_res;
        diff -= floor(diff + 0.5);
        ASSERT_LE(fabs(diff), quotient / 2.);
      }
    }
  }
  delete_module_info(module);
}

struct message_and_noise {
  std::vector<int64_t> decrypted;
  std::vector<double> noise;
  double noise_amplitude;
};

message_and_noise decrypt(const MODULE* module, uint64_t message_basebit, const onionpir_secret_key& skey,
                          const int64_t* rlwe, uint64_t rlwe_size) {
  uint64_t b_size = rlwe_size / 2;
  uint64_t a_size = (rlwe_size + 1) / 2;
  std::vector<int64_t> decrypted(b_size * ONIONPIR_N);
  std::vector<double> noise(ONIONPIR_N);
  rlwe_phase_base2k((MODULE*)module, ONIONPIR_K,               //
                    decrypted.data(), b_size, ONIONPIR_N,      //
                    skey.ppol_s,                               //
                    rlwe, a_size, 2 * ONIONPIR_N,              //
                    rlwe + ONIONPIR_N, b_size, 2 * ONIONPIR_N  //
  );
  uint64_t final_size = (message_basebit + ONIONPIR_K - 1) / ONIONPIR_K;
  round_polynomial_noise(ONIONPIR_N, ONIONPIR_K, message_basebit,  //
                         noise.data(),                             //
                         decrypted.data(), final_size, ONIONPIR_N, decrypted.data(), b_size, ONIONPIR_N);
  decrypted.resize(final_size * ONIONPIR_N);
  double max_noise = 0;
  for (uint64_t i = 0; i < ONIONPIR_N; ++i) {
    double n = fabs(noise[i]);
    if (n > max_noise) max_noise = n;
  }
  message_and_noise res;
  res.decrypted = std::move(decrypted);
  res.noise = std::move(noise);
  res.noise_amplitude = max_noise;
  return res;
}

void check_1D_zero(int64_t* v, uint64_t size) {
  for (uint64_t i = 0; i < size; ++i) {
    ASSERT_EQ(v[i], 0) << "non-zero-pos: " << i;
  }
}

void check_2D_zero(int64_t* v, uint64_t rows, uint64_t cols) {
  for (uint64_t i = 0; i < rows; ++i) {
    for (uint64_t j = 0; j < cols; ++j) {
      ASSERT_EQ(v[i * cols + j], 0) << "non-zero-pos: " << i << "," << j;
    }
  }
}

void check_normalized(int64_t* v, uint64_t size, uint64_t log2_base2k) {
  const int64_t lb = -(1l << (log2_base2k - 1));
  const int64_t ub = (1l << (log2_base2k - 1));
  for (uint64_t i = 0; i < size; ++i) {
    ASSERT_GE(v[i], lb);
    ASSERT_LT(v[i], ub);
  }
}

TEST(onionpir, generate_query_phase1) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  uint64_t row = 11;
  onionpir_input_query query;
  onionpir_secret_key skey;
  // TODO keygen -> make a separate function?
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_query_phase1(module.get(), row, query, skey);
  // verify that the phase of the query is correct
  ASSERT_EQ(query.phase1.size(), ONIONPIR_query1_ncols * ONIONPIR_N);
  message_and_noise decr = decrypt(module.get(), 80, skey, query.phase1.data(), ONIONPIR_query1_ncols);
  double noise_bound = pow(2., -115);  // TODO
  ASSERT_GT(decr.noise_amplitude, 0.);
  ASSERT_LE(decr.noise_amplitude, noise_bound);
  int64_t EXPECTED_VALUE = (1L << ONIONPIR_K) / ONIONPIR_N;
  for (uint64_t i = 0; i < 4; ++i) {
    int64_t& v = decr.decrypted[(i + 1) * ONIONPIR_N + (4 * row + i)];
    ASSERT_EQ(v, EXPECTED_VALUE);
    v = 0;
  }
  // the rest should be zero
  check_2D_zero(decr.decrypted.data(), 5, ONIONPIR_N);
}

TEST(onionpir, generate_queryexp_phase1) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  uint64_t row = 7;
  onionpir_expanded_query query;
  onionpir_secret_key skey;
  // TODO keygen -> make a separate function?
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_queryexp_phase1(module.get(), row, query, skey);
  // verify that the phase of the query is correct
  ASSERT_NE(query.query_exp_phase1, nullptr);
#ifdef TEST_MODE
  double noise_bound = pow(2., -87);  // TODO
  ASSERT_EQ(query.query_exp_phase1_raw.size(), ONIONPIR_query1exp_nrows * ONIONPIR_query1exp_ncols * ONIONPIR_N);
  int64_t* q = query.query_exp_phase1_raw.data();
  for (uint64_t r = 0; r < ONIONPIR_query1exp_nrows; ++r) {
    message_and_noise decr =
        decrypt(module.get(), 64, skey, q + r * ONIONPIR_query1exp_ncols * ONIONPIR_N, ONIONPIR_query1exp_ncols);
    ASSERT_GT(decr.noise_amplitude, 0.);
    ASSERT_LE(decr.noise_amplitude, noise_bound);
    for (uint64_t i = 0; i < 4; ++i) {
      if (r == 4 * row + i) {
        int64_t& v = decr.decrypted[i * ONIONPIR_N + 0];
        ASSERT_EQ(v, 1);
        v = 0;
      }
    }
    // the rest should be zero
    check_2D_zero(decr.decrypted.data(), 4, ONIONPIR_N);
  }
#endif
}

static void check_prgsw_of(const MODULE* module, const std::vector<int64_t> prgsw, const std::vector<int64_t> mu,
                           const onionpir_secret_key& skey, const uint64_t nrows, const uint64_t ncols,
                           const uint64_t message_basebit, const double noise_bound) __attribute((unused));

static void check_prgsw_of(const MODULE* module, const std::vector<int64_t> prgsw, const std::vector<int64_t> mu,
                           const onionpir_secret_key& skey, const uint64_t nrows, const uint64_t ncols,
                           const uint64_t message_basebit, const double noise_bound) {
  ASSERT_EQ(prgsw.size(), nrows * ncols * ONIONPIR_N);
  const int64_t* q = prgsw.data();
  const uint64_t ell = (message_basebit + ONIONPIR_K - 1) / ONIONPIR_K;
  for (uint64_t r = 0; r < ell; ++r) {
    message_and_noise decr = decrypt(module, message_basebit, skey, q + r * ncols * ONIONPIR_N, ncols);
    ASSERT_GT(decr.noise_amplitude, 0.);  // noiseless?
    ASSERT_LE(decr.noise_amplitude, noise_bound);
    // row r should be s^2
    for (uint64_t j = 0; j < ONIONPIR_N; ++j) {
      int64_t& v = decr.decrypted[r * ONIONPIR_N + j];
      ASSERT_EQ(v, mu[j]);
      v = 0;
    }
    // the rest should be zero
    check_2D_zero(decr.decrypted.data(), ell, ONIONPIR_N);
  }
}

TEST(onionpir, generate_cloud_key_phase2) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_cloud_key ckey;
  onionpir_secret_key skey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase2(module.get(), ckey, skey);
  // verify that the phase of the query is correct
  ASSERT_NE(ckey.rk_s, nullptr);
#ifdef TEST_MODE
  std::vector<int64_t> s2(ONIONPIR_N);
  znx_small_single_product(module.get(), s2.data(), skey.s.data(), skey.s.data(),
                           get_tmp_space(znx_small_single_product_tmp_bytes(module.get())));
  double noise_bound = pow(2., -116);  // TODO
  check_prgsw_of(module.get(), ckey.rk_s_raw, s2, skey, ONIONPIR_rk_nrows, ONIONPIR_rk_ncols, 96, noise_bound);
#endif
}

TEST(onionpir, generate_cloud_key_phase1) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_cloud_key ckey;
  onionpir_secret_key skey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);
  std::vector<int64_t> minus_autom_s(ONIONPIR_N);
  // verify that the phase of the query is correct
  for (uint64_t i = 0; i < 12; ++i) {
    int64_t p = (4096 >> i) + 1;
    vec_znx_automorphism(module.get(), p,                      //
                         minus_autom_s.data(), 1, ONIONPIR_N,  //
                         skey.s.data(), 1, ONIONPIR_N);
    vec_znx_sub(module.get(),                         // negate
                minus_autom_s.data(), 1, ONIONPIR_N,  //
                minus_autom_s.data(), 0, ONIONPIR_N,  //
                minus_autom_s.data(), 1, ONIONPIR_N);
    ASSERT_NE(ckey.autom.at(p), nullptr);
#ifdef TEST_MODE
    double noise_bound = pow(2., -116);  // TODO
    check_prgsw_of(module.get(), ckey.autom_raw.at(p), minus_autom_s, skey, ONIONPIR_automkey_nrows,
                   ONIONPIR_automkey_ncols, 96, noise_bound);
#endif
  }
}

TEST(tinyfhe, apply_automorphism) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_cloud_key ckey;
  onionpir_secret_key skey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);
  std::vector<int64_t> minus_autom_s(ONIONPIR_N);
  int64_t p = 3;  // ID todo
  // generate some random plaintext
  uint64_t mu_size = 6;
  uint64_t a_size = 8;
  uint64_t b_size = 7;
  uint64_t rlwe_size = a_size + b_size;
  uint64_t res_rlwe_size = 16;
  std::vector<int64_t> message(mu_size * ONIONPIR_N, 0);
  random_centered_reduced(message.size(), ONIONPIR_K, message.data());
  std::vector<int64_t> rlwe(rlwe_size * ONIONPIR_N, 0);
  for (uint64_t i = 0; i < mu_size; ++i) {
    memcpy(rlwe.data() + (2 * i + 1) * ONIONPIR_N, message.data() + i * ONIONPIR_N, ONIONPIR_N * sizeof(int64_t));
  }
  random_log2bound_symmetric(ONIONPIR_N, 2, rlwe.data() + (2 * b_size - 1) * ONIONPIR_N);
  onionpir_rlwe_trivial_encrypt_inplace(module.get(), rlwe.data(), rlwe_size, skey);
  // apply automorphism on message
  std::vector<int64_t> autom_mes(mu_size * ONIONPIR_N);
  vec_znx_automorphism(module.get(), p,                        //
                       autom_mes.data(), mu_size, ONIONPIR_N,  //
                       message.data(), mu_size, ONIONPIR_N);
  // apply automorphism
  std::vector<int64_t> res_rlwe(res_rlwe_size * ONIONPIR_N, 0);
  apply_automorphism(module.get(), p, ONIONPIR_K, res_rlwe.data(), res_rlwe_size, rlwe.data(), rlwe_size,
                     ckey.autom.at(p).get(), ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols);
  // decrypt message
  message_and_noise decr = decrypt(module.get(), mu_size * ONIONPIR_K, skey, res_rlwe.data(), res_rlwe_size);
  ASSERT_GT(decr.noise_amplitude, 0.);  // noiseless?
  ASSERT_LE(decr.noise_amplitude, pow(2., -96));
  std::cerr << "new-noise: " << log2(decr.noise_amplitude) << std::endl;
  for (uint64_t i = 0; i < mu_size; ++i) {
    for (uint64_t j = 0; j < ONIONPIR_N; ++j) {
      int64_t actual = decr.decrypted[i * ONIONPIR_N + j];
      int64_t expect = autom_mes[i * ONIONPIR_N + j];
      ASSERT_EQ(actual, expect);
    }
  }
}

EXPORT void onionpir_online_phase1(const MODULE* module,  // N
                                   onionpir_phase1_results& result, const uint8_t* plaintext_db,
                                   const onionpir_expanded_query& qexp, uint64_t db_ncols);
TEST(onionpir, online_phase1) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  uint64_t db_ncols = 16;
  // generate a random db
  std::vector<uint8_t> plaintext_db(7 * ONIONPIR_db_nrows * db_ncols * ONIONPIR_N);
  std::uniform_int_distribution<uint8_t> randu8;
  for (uint8_t& b : plaintext_db) b = randu8(randgen());
  //
  uint64_t row = 7;
  onionpir_secret_key skey;
  onionpir_expanded_query qexp;
  onionpir_phase1_results res;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_queryexp_phase1(module.get(), row, qexp, skey);
  // exec phase 1
  uint64_t tbeg = std::chrono::steady_clock::now().time_since_epoch().count();
  onionpir_online_phase1(module.get(), res, plaintext_db.data(), qexp, db_ncols);
  uint64_t tend = std::chrono::steady_clock::now().time_since_epoch().count();
  double rtime = double(tend - tbeg) / 1e9;
  double plainbytes = 1024 * 7 * ONIONPIR_N * db_ncols;
  std::cerr << "online phase dotp time (s): " << rtime << std::endl;
  std::cerr << "db MB/s (strict mode): " << plainbytes / rtime / 1e6 << std::endl;
  std::cerr << "(x2 or x3 if preprocessing is allowed)" << std::endl;
  // test
  for (uint64_t i = 0; i < db_ncols; ++i) {
    message_and_noise decr = decrypt(module.get(), 56, skey,                                    //
                                     res.res.data() + i * ONIONPIR_results1_size * ONIONPIR_N,  //
                                     ONIONPIR_results1_size);
    std::cerr << "noise: " << log2(decr.noise_amplitude) << std::endl;
  }
}

void onionpir_prgsw_to_rgsw(const MODULE* module, const onionpir_cloud_key& ckey, int64_t* out_rgsw, uint64_t out_nrows,
                            uint64_t out_ncols, const int64_t* in_prgsw, uint64_t in_nrows, uint64_t in_ncols);

TEST(onionpir, prgsw_to_rgsw) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_secret_key skey;
  onionpir_cloud_key ckey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase2(module.get(), ckey, skey);
  // create a random rlwe matrix
  uint64_t in_nrows = 5;
  uint64_t in_ncols = 2 * ONIONPIR_rk_nrows;
  uint64_t out_nrows = 2 * in_nrows;
  uint64_t out_ncols = ONIONPIR_rk_ncols;
  std::vector<int64_t> in(in_nrows * in_ncols * ONIONPIR_N, 0);
  std::vector<int64_t> out(out_nrows * out_ncols * ONIONPIR_N);
  for (uint64_t i = 0; i < in_nrows; ++i) {
    for (uint64_t j = 0; j < in_ncols / 2; ++j) {
      in[(i * in_ncols + 2 * j + 1) * ONIONPIR_N] = 1;
    }
  }
  for (uint64_t i = 0; i < in_nrows; ++i) {
    onionpir_rlwe_trivial_encrypt_inplace(module.get(),                           //
                                          in.data() + i * in_ncols * ONIONPIR_N,  //
                                          in_ncols, skey);
  }
  onionpir_prgsw_to_rgsw(module.get(), ckey,                //
                         out.data(), out_nrows, out_ncols,  //
                         in.data(), in_nrows, in_ncols);
  for (uint64_t i = 0; i < out_nrows; ++i) {
    message_and_noise decr = decrypt(module.get(), 112, skey,                  //
                                     out.data() + i * out_ncols * ONIONPIR_N,  //
                                     out_ncols);
    uint64_t max_depth = ONIONPIR_N > 100 ? 5 : 6;
    if (i % 2 == 0) {
      // message should be -s on all rows
      for (uint64_t j = 0; j < max_depth; ++j) {
        for (uint64_t k = 0; k < ONIONPIR_N; ++k) {
          ASSERT_EQ(decr.decrypted[j * ONIONPIR_N + k], -skey.s[k]);
        }
      }
    } else {
      // message should be 1 on all rows
      for (uint64_t j = 0; j < max_depth; ++j) {
        for (uint64_t k = 0; k < ONIONPIR_N; ++k) {
          ASSERT_EQ(decr.decrypted[j * ONIONPIR_N + k], k == 0);
        }
      }
    }
  }
}

void onionpir_cmux_eval(const MODULE* module, int64_t* out_rlwe, uint64_t out_size, const int64_t* c1_rlwe,
                        const int64_t* c0_rlwe, const VMP_PMAT* rgsw);

TEST(onionpir, onionpir_cmux_eval) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_secret_key skey;
  for (uint64_t sbit : {0, 1}) {
    onionpir_generate_secret_key(module.get(), skey);
    // create a random rlwe matrix
    constexpr uint64_t _2P = ONIONPIR_query2exp_ncols;
    constexpr uint64_t P = _2P / 2;
    constexpr uint64_t _2S = ONIONPIR_query2exp_nrows;
    // constexpr uint64_t S = _2S/2;
    constexpr uint64_t N = ONIONPIR_N;
    std::vector<int64_t> message_raw(4 * N, 0);
    std::vector<int64_t> c0_raw(_2S * N, 0);
    std::vector<int64_t> c1_raw(_2S * N, 0);
    std::vector<int64_t> out_raw(_2S * N, 0);
    std::vector<int64_t> rgsw_raw(_2P * _2P * N, 0);
    int64_t(*rgsw)[_2P][N] = (int64_t(*)[_2P][N])rgsw_raw.data();
    int64_t(*c0)[N] = (int64_t(*)[N])c0_raw.data();
    int64_t(*c1)[N] = (int64_t(*)[N])c1_raw.data();
    int64_t(*out)[N] = (int64_t(*)[N])out_raw.data();
    int64_t(*message)[N] = (int64_t(*)[N])message_raw.data();
    for (uint64_t i = 0; i < _2P; ++i) {
      rgsw[i][i][0] = sbit;
      random_log2bound_symmetric(N, 5, rgsw[i][2 * P - 1]);
      onionpir_rlwe_trivial_encrypt_inplace(module.get(), *rgsw[i], _2P, skey, false);
    }
    for (uint64_t i = 0; i < 4; ++i) {
      random_centered_reduced(N, ONIONPIR_K, c0[2 * i + 1]);
      random_centered_reduced(N, ONIONPIR_K, c1[2 * i + 1]);
      if (sbit) {
        memcpy(message[i], c1[2 * i + 1], N * sizeof(int64_t));
      } else {
        memcpy(message[i], c0[2 * i + 1], N * sizeof(int64_t));
      }
    }
    random_log2bound_symmetric(N, 5, c0[_2S - 1]);
    random_log2bound_symmetric(N, 5, c1[_2S - 1]);
    onionpir_rlwe_trivial_encrypt_inplace(module.get(), *c0, _2S, skey);
    onionpir_rlwe_trivial_encrypt_inplace(module.get(), *c1, _2S, skey);
    VMP_PMAT_UNIPTR pmat(new_vmp_pmat(module.get(), _2S, _2P));
    uint8_t* tmp_bytes = get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module.get(), _2P, _2P));
    vmp_prepare_contiguous(module.get(), pmat.get(),  //
                           **rgsw, _2S, _2P, tmp_bytes);
    onionpir_cmux_eval(module.get(), *out, _2S, *c1, *c0, pmat.get());
    message_and_noise decr = decrypt(module.get(), 64, skey,  //
                                     *out, _2S);
    std::cerr << "noise: " << log2(decr.noise_amplitude) << std::endl;
    for (uint64_t i = 0; i < 4; ++i) {
      for (uint64_t j = 0; j < N; ++j) {
        ASSERT_EQ(decr.decrypted[i * ONIONPIR_N + j], message[i][j]);
      }
    }
  }
}

TEST(onionpir, onionpir_generate_query_phase2) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_secret_key skey;
  onionpir_input_query query;
  uint64_t col = 7;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_query_phase2(module.get(), col, query, skey);

  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t K = ONIONPIR_K;
  constexpr uint64_t ELLT = ONIONPIR_query2exp_ellt;
  constexpr uint64_t out_size = ONIONPIR_query2_ncols;
  ASSERT_EQ(query.phase2.size(), N * out_size);
  message_and_noise decr = decrypt(module.get(), 7 + K * ELLT, skey,  //
                                   query.phase2.data(), out_size);
  std::cerr << "noise: " << log2(decr.noise_amplitude) << std::endl;

  const uint64_t db_ncols_max_bits = 3;
  for (uint64_t k = 0; k < db_ncols_max_bits; ++k) {
    const uint64_t colbit = (col >> k) & 1;
    for (uint64_t i = 0; i < ELLT; ++i) {
      int64_t& v = decr.decrypted[(i + 1) * N + k * ELLT + i];
      ASSERT_EQ(v, colbit * (1L << K) / ONIONPIR_query2_denom) << k << " " << i << " " << colbit;
      v = 0;
    }
  }
  check_2D_zero(decr.decrypted.data(), ELLT + 1, N);
}

TEST(onionpir, onionpir_generate_queryexp_phase2) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_secret_key skey;
  onionpir_expanded_query qexp;
  uint64_t col = 5;
  onionpir_generate_secret_key(module.get(), skey);
  /*
  onionpir_generate_queryexp_phase2(module.get(), col, qexp, skey);

  constexpr uint64_t N = ONIONPIR_N;
#ifdef TEST_MODE
  constexpr uint64_t K = ONIONPIR_K;
  constexpr uint64_t ELLT = ONIONPIR_query2exp_ellt;
  constexpr uint64_t out_ncols = ONIONPIR_query2exp_ncols;
  constexpr uint64_t out_nrows = ONIONPIR_query2exp_nrows;
#endif
  // create two constants
  std::vector<int64_t> ONE(N, 0);
  std::vector<int64_t> MINUS_S(N);
  ONE[0] = 1;
  vec_znx_negate(module.get(),          //
                 MINUS_S.data(), 1, N,  //
                 skey.s.data(), 1, N);
  for (uint64_t idx = 0; idx < ONIONPIR_query2exp_nb; ++idx) {
    // check that we have a prgsw of colbit
    ASSERT_NE(qexp.query_exp_phase2[idx], nullptr);
#ifdef TEST_MODE
    uint64_t colbit = (col >> idx) & 1;
    int64_t(*q)[out_ncols][N] = (int64_t(*)[out_ncols][N])qexp.query_exp_phase2_raw[idx].data();
    for (uint64_t i = 0; i < out_nrows; ++i) {
      message_and_noise decr = decrypt(module.get(), K * ELLT, skey,  //
                                       *q[i], out_ncols);
      ASSERT_LE(decr.noise_amplitude, pow(2., -80));
      if (colbit == 0) {
        check_2D_zero(decr.decrypted.data(), ELLT, ONIONPIR_N);
      } else {
        for (uint64_t j = 0; j < N; ++j) {
          int64_t& vj = decr.decrypted[(i / 2) * N + j];
          if (i % 2 == 1) {
            ASSERT_EQ(vj, ONE[j]);
          } else {
            ASSERT_EQ(vj, MINUS_S[j]);
          }
          vj = 0;
        }
        check_2D_zero(decr.decrypted.data(), ELLT, ONIONPIR_N);
      }
    }
#endif
  }
  */
}

TEST(onionpir, online_phase2) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_secret_key skey;
  onionpir_expanded_query qexp;
  onionpir_phase1_results rin;
  onionpir_phase2_results rout;
  uint64_t db_ncols = ONIONPIR_db_ncols;
  uint64_t col = 5;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_queryexp_phase2(module.get(), col, qexp, skey);
  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t RLWE_size = ONIONPIR_results1_size;
  constexpr uint64_t b_size = RLWE_size / 2;

  rin.res.resize(db_ncols * RLWE_size * N, 0);
  int64_t(*in)[RLWE_size][N] = (int64_t(*)[RLWE_size][N])rin.res.data();
  std::vector<uint8_t> cleartext(7 * N);
  std::vector<int64_t> plaintext(4 * N);
  std::vector<int64_t> golden_plaintext(4 * N);
  std::uniform_int_distribution<uint8_t> uniform_u8;
  for (uint64_t c = 0; c < db_ncols; ++c) {
    // encrypt a random message
    for (uint8_t& b : cleartext) b = uniform_u8(randgen());
    onionpir_encode_plaintext(module.get(), plaintext.data(), cleartext.data());
    vec_znx_copy(module.get(),               //
                 *in[c] + N, b_size, 2 * N,  //
                 plaintext.data(), 4, N);
    // add noise (env 70 bits)
    random_log2bound_symmetric(N, 5, in[c][2 * b_size - 1]);
    // encrypt
    onionpir_rlwe_trivial_encrypt_inplace(module.get(),  //
                                          *in[c], RLWE_size, skey);
    if (c == col) golden_plaintext = plaintext;
  }
  // apply online phase
  double t0 = std::chrono::steady_clock::now().time_since_epoch().count();
  onionpir_online_phase2(module.get(), rout, rin, qexp, db_ncols);
  double t1 = std::chrono::steady_clock::now().time_since_epoch().count();
  std::cerr << "eval phase2 on " << db_ncols << " columns: " << (t1 - t0) / 1e9 << " seconds" << std::endl;
  // memset(rout.res.data()+9*N, 0, 8*N);
  // check result
  message_and_noise decr = decrypt(module.get(), 56, skey, rout.res.data(), ONIONPIR_results2_size);
  std::cerr << "decr noise: " << log2(decr.noise_amplitude) << std::endl;
  for (uint64_t i = 0; i < 4; ++i) {
    const int64_t* expect = golden_plaintext.data() + i * N;
    const int64_t* actual = decr.decrypted.data() + i * N;
    for (uint64_t j = 0; j < N; ++j) {
      ASSERT_EQ(expect[j], actual[j]);
    }
  }
}

TEST(onionpir, final_decrypt) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_secret_key skey;
  onionpir_phase2_results rout;
  onionpir_generate_secret_key(module.get(), skey);
  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t RLWE_size = ONIONPIR_results2_size;
  constexpr uint64_t b_size = RLWE_size / 2;
  std::vector<uint8_t> cleartext(7 * N);
  std::vector<int64_t> plaintext(4 * N);
  std::vector<int64_t> rlwe_raw(RLWE_size * N);
  int64_t(*rlwe)[N] = (int64_t(*)[N])rlwe_raw.data();
  std::uniform_int_distribution<uint8_t> uniform_u8;
  // encrypt a random message
  for (uint8_t& b : cleartext) b = uniform_u8(randgen());
  onionpir_encode_plaintext(module.get(), plaintext.data(), cleartext.data());
  // add noise (env 63 bits)
  random_log2bound_symmetric(N, 2, rlwe[2 * b_size - 1]);
  // add plaintext
  vec_znx_add(module.get(),              //
              *rlwe + N, b_size, 2 * N,  //
              *rlwe + N, b_size, 2 * N,  //
              plaintext.data(), 4, N);
  // encrypt
  onionpir_rlwe_trivial_encrypt_inplace(module.get(),  //
                                        *rlwe, RLWE_size, skey);
  rout.res = rlwe_raw;
  // apply online phase
  std::vector<uint8_t> decrypt(7 * N);
  onionpir_final_decrypt(module.get(), decrypt.data(), skey, rout);
  // check result
  for (uint64_t i = 0; i < 7 * N; i++) {
    ASSERT_EQ(decrypt[i], cleartext[i]);
  }
}

/** returns centermod and carry */
std::pair<int64_t, int64_t> centermod_and_carry(const int64_t inp) {
  int64_t carry = (inp + (1l << (ONIONPIR_K - 1))) >> ONIONPIR_K;
  return {inp - (carry << ONIONPIR_K), carry};
}

void rlwe_trace_expand_1step_check(int64_t* inp, int64_t* res, uint64_t size, bool odd_coeff) {
  // even coefficients are 2x input even or odd coefficients
  for (uint64_t j = 0; j < ONIONPIR_N; j += 2) {
    int64_t ci = 0;
    for (int64_t i = size - 1; i >= 0; --i) {
      const auto exp_co = centermod_and_carry(ci + inp[i * ONIONPIR_N + j + (int64_t)odd_coeff] * 2);
      ci = exp_co.second;
      int64_t act = res[i * ONIONPIR_N + j];
      ASSERT_EQ(act, exp_co.first);
    }
  }

  // odd coefficients are zero
  for (uint64_t i = 0; i < size; ++i) {
    for (uint64_t j = 1; j < ONIONPIR_N; j += 2) {
      ASSERT_EQ(res[i * ONIONPIR_N + j], 0);
    }
  }
}

TEST(tinyfhe, rlwe_trace_expand_1step) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_cloud_key ckey;
  onionpir_secret_key skey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);

  int64_t autom_p = ONIONPIR_N + 1;
  int64_t rotate_p = -1;

  uint64_t mu_size = 6;
  uint64_t a_size = 8;
  uint64_t b_size = 7;
  uint64_t rlwe_size = a_size + b_size;
  uint64_t res_rlwe_size = 16;

  // generate some random plaintext
  std::vector<int64_t> message(mu_size * ONIONPIR_N, 0);
  random_centered_reduced(message.size(), ONIONPIR_K, message.data());

  // encrypt into RLWE
  std::vector<int64_t> rlwe(rlwe_size * ONIONPIR_N, 0);
  for (uint64_t i = 0; i < mu_size; ++i) {
    memcpy(rlwe.data() + (2 * i + 1) * ONIONPIR_N, message.data() + i * ONIONPIR_N, ONIONPIR_N * sizeof(int64_t));
  }
  random_log2bound_symmetric(ONIONPIR_N, 2, rlwe.data() + (2 * b_size - 1) * ONIONPIR_N);
  onionpir_rlwe_trivial_encrypt_inplace(module.get(), rlwe.data(), rlwe_size, skey);

  std::vector<int64_t> res0_rlwe(res_rlwe_size * ONIONPIR_N, 0);
  std::vector<int64_t> res1_rlwe(res_rlwe_size * ONIONPIR_N, 0);
  rlwe_trace_expand_1step(module.get(), autom_p, rotate_p, ONIONPIR_K, res0_rlwe.data(), res_rlwe_size,
                          res1_rlwe.data(), res_rlwe_size, rlwe.data(), rlwe_size, ckey.autom.at(autom_p).get(),
                          ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols);
  {
    // decrypt and check res0 = in + in(X^p)
    message_and_noise decr = decrypt(module.get(), mu_size * ONIONPIR_K, skey, res0_rlwe.data(), res_rlwe_size);
    ASSERT_GT(decr.noise_amplitude, 0.);  // noiseless?
    ASSERT_LE(decr.noise_amplitude, pow(2., -96));
    // std::cerr << "output noise: " << log2(decr.noise_amplitude) << std::endl;

    // check result is reduced
    check_normalized(decr.decrypted.data(), mu_size * ONIONPIR_N, ONIONPIR_K);

    // check even coefficients are 2x input even coefficients and odd coefficients are zero
    rlwe_trace_expand_1step_check(message.data(), decr.decrypted.data(), mu_size, false);
  }

  {
    // decrypt and check res1 = rotate(in - in(X^p))
    message_and_noise decr = decrypt(module.get(), mu_size * ONIONPIR_K, skey, res1_rlwe.data(), res_rlwe_size);
    ASSERT_GT(decr.noise_amplitude, 0.);  // noiseless?
    ASSERT_LE(decr.noise_amplitude, pow(2., -96));
    // std::cerr << "output noise: " << log2(decr.noise_amplitude) << std::endl;

    // check result is reduced
    check_normalized(decr.decrypted.data(), mu_size * ONIONPIR_N, ONIONPIR_K);

    // check even coefficients are 2x input odd coefficients and odd coefficients are zero
    rlwe_trace_expand_1step_check(message.data(), decr.decrypted.data(), mu_size, true);
  }
}

TEST(tinyfhe, rlwe_trace_expand_1step_no_rotation) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_cloud_key ckey;
  onionpir_secret_key skey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);

  int64_t autom_p = ONIONPIR_N + 1;

  uint64_t mu_size = 6;
  uint64_t a_size = 8;
  uint64_t b_size = 7;
  uint64_t rlwe_size = a_size + b_size;
  uint64_t res_rlwe_size = 16;

  // generate some random plaintext
  std::vector<int64_t> message(mu_size * ONIONPIR_N, 0);
  random_centered_reduced(message.size(), ONIONPIR_K, message.data());

  // encrypt into RLWE
  std::vector<int64_t> rlwe(rlwe_size * ONIONPIR_N, 0);
  for (uint64_t i = 0; i < mu_size; ++i) {
    memcpy(rlwe.data() + (2 * i + 1) * ONIONPIR_N, message.data() + i * ONIONPIR_N, ONIONPIR_N * sizeof(int64_t));
  }
  random_log2bound_symmetric(ONIONPIR_N, 2, rlwe.data() + (2 * b_size - 1) * ONIONPIR_N);
  onionpir_rlwe_trivial_encrypt_inplace(module.get(), rlwe.data(), rlwe_size, skey);

  std::vector<int64_t> res0_rlwe(res_rlwe_size * ONIONPIR_N, 0);
  rlwe_trace_expand_1step_no_rotation(module.get(), autom_p, ONIONPIR_K, res0_rlwe.data(), res_rlwe_size, rlwe.data(),
                                      rlwe_size, ckey.autom.at(autom_p).get(), ONIONPIR_automkey_nrows,
                                      ONIONPIR_automkey_ncols);
  {
    // decrypt and check res0 = in + in(X^p)
    message_and_noise decr = decrypt(module.get(), mu_size * ONIONPIR_K, skey, res0_rlwe.data(), res_rlwe_size);
    ASSERT_GT(decr.noise_amplitude, 0.);  // noiseless?
    ASSERT_LE(decr.noise_amplitude, pow(2., -96));
    // std::cerr << "output noise: " << log2(decr.noise_amplitude) << std::endl;

    // check result is reduced
    check_normalized(decr.decrypted.data(), mu_size * ONIONPIR_N, ONIONPIR_K);

    // check even coefficients are 2x input even coefficients and odd coefficients are zero
    rlwe_trace_expand_1step_check(message.data(), decr.decrypted.data(), mu_size, false);
  }
}

TEST(onionpir, onionpir_trace_expand) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  onionpir_cloud_key ckey;
  onionpir_secret_key skey;
  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);

  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t K = ONIONPIR_K;
  constexpr uint64_t mu_size = 5;
  constexpr uint64_t a_size = 8;
  constexpr uint64_t b_size = 7;
  constexpr uint64_t rlwe_size = a_size + b_size;
  constexpr uint64_t res_nrows = 13;
  constexpr uint64_t res_ncols = 16;

  // output is scaled by 2^#steps
  const int64_t steps = ceil(log2(res_nrows));
  const int64_t res_mult_coef = 1l << steps;

  // generate some random plaintext in the first res_rlwe_nrows positions
  std::vector<int64_t> message(mu_size * ONIONPIR_N, 0);
  for (uint64_t i = 0; i < mu_size; ++i) {
    random_centered_reduced(res_nrows, ONIONPIR_K,  //
                            message.data() + i * ONIONPIR_N);
  }

  // encrypt into RLWE
  std::vector<int64_t> rlwe(rlwe_size * ONIONPIR_N, 0);
  vec_znx_copy(module.get(),                    //
               rlwe.data() + N, b_size, 2 * N,  //
               message.data(), mu_size, N);
  random_log2bound_symmetric(ONIONPIR_N, 2, rlwe.data() + (2 * b_size - 1) * ONIONPIR_N);
  onionpir_rlwe_trivial_encrypt_inplace(module.get(), rlwe.data(), rlwe_size, skey);

  // multiplied message
  std::vector<int64_t> mulmessage = message;
  for (int64_t& x : mulmessage) x *= res_mult_coef;
  vec_znx_normalize_base2k(module.get(), K, mulmessage.data(), mu_size, N, mulmessage.data(), mu_size, N,
                           get_tmp_space(vec_znx_normalize_base2k_tmp_bytes(module.get())));

  // call tested function
  std::vector<int64_t> res_rlwes(res_nrows * res_ncols * ONIONPIR_N, 0);
  onionpir_trace_expand(module.get(), res_rlwes.data(), res_nrows, res_ncols, rlwe.data(), rlwe_size, ckey);

  for (uint64_t i = 0; i < res_nrows; ++i) {
    message_and_noise decr = decrypt(module.get(), mu_size * ONIONPIR_K - steps, skey,  //
                                     res_rlwes.data() + i * res_ncols * ONIONPIR_N, res_ncols);
    ASSERT_GT(decr.noise_amplitude, 0.);  // noiseless?
    ASSERT_LE(decr.noise_amplitude, pow(2., -90));
    std::cerr << "output noise: " << log2(decr.noise_amplitude) << std::endl;

    // check result is reduced
    for (uint64_t j=0; j<mu_size; ++j) {
      int64_t& v = decr.decrypted[j*N];
      ASSERT_EQ(v, mulmessage[j*N+i]);
      v=0;
    }
    check_2D_zero(decr.decrypted.data(), mu_size, ONIONPIR_N);
  }
}

TEST(onionpir, onionpir_query_expand_phase1) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  uint64_t row = 11;
  onionpir_input_query query_inp;
  onionpir_expanded_query query;
  onionpir_secret_key skey;
  onionpir_cloud_key ckey;

  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);

  // generate query
  onionpir_generate_query_phase1(module.get(), row, query_inp, skey);

  // expand query
  double t0 = std::chrono::steady_clock::now().time_since_epoch().count();
  onionpir_query_expand_phase1(module.get(), query, ckey, query_inp);
  double t1 = std::chrono::steady_clock::now().time_since_epoch().count();
  std::cerr << "query_phase1_expand time (s): " << (t1 - t0) / 1e9 << std::endl;

  // verify that the phase of the query is correct
  ASSERT_NE(query.query_exp_phase1, nullptr);

#ifdef TEST_MODE
  double noise_bound = pow(2., -87);  // TODO
  ASSERT_EQ(query.query_exp_phase1_raw.size(), ONIONPIR_query1exp_nrows * ONIONPIR_query1exp_ncols * ONIONPIR_N);
  int64_t* q = query.query_exp_phase1_raw.data();
  for (uint64_t r = 0; r < ONIONPIR_query1exp_nrows; ++r) {
    message_and_noise decr =
        decrypt(module.get(), 64, skey, q + r * ONIONPIR_query1exp_ncols * ONIONPIR_N, ONIONPIR_query1exp_ncols);
    ASSERT_GT(decr.noise_amplitude, 0.);
    ASSERT_LE(decr.noise_amplitude, noise_bound);
    for (uint64_t i = 0; i < 4; ++i) {
      if (r == 4 * row + i) {
        int64_t& v = decr.decrypted[i * ONIONPIR_N + 0];
        ASSERT_EQ(v, 1);
        v = 0;
      }
    }
    // the rest should be zero
    check_2D_zero(decr.decrypted.data(), 4, ONIONPIR_N);
  }
#endif
}

TEST(onionpir, onionpir_query_expand_phase2) {
  MODULE_UNIPTR module(new_module_info(ONIONPIR_N, FFT64));
  uint64_t col = 7;
  onionpir_input_query query;
  onionpir_expanded_query qexp;
  onionpir_secret_key skey;
  onionpir_cloud_key ckey;

  onionpir_generate_secret_key(module.get(), skey);
  onionpir_generate_cloud_key_phase1(module.get(), ckey, skey);
  onionpir_generate_cloud_key_phase2(module.get(), ckey, skey);

  // generate query
  onionpir_generate_query_phase2(module.get(), col, query, skey);

  // expand query
  double t0 = std::chrono::steady_clock::now().time_since_epoch().count();
  onionpir_query_expand_phase2(module.get(), qexp, ckey, query);
  double t1 = std::chrono::steady_clock::now().time_since_epoch().count();
  std::cerr << "query_phase2_expand time (s): " << (t1 - t0) / 1e9 << std::endl;

  //   constexpr uint64_t N = ONIONPIR_N;
  // #ifdef TEST_MODE
  //   constexpr uint64_t K = ONIONPIR_K;
  //   constexpr uint64_t ELLT = ONIONPIR_query2exp_ellt;
  //   constexpr uint64_t out_ncols = ONIONPIR_query2exp_ncols;
  //   constexpr uint64_t out_nrows = ONIONPIR_query2exp_nrows;
  // #endif
  //   // create two constants
  //   std::vector<int64_t> ONE(N, 0);
  //   std::vector<int64_t> MINUS_S(N);
  //   ONE[0] = 1;
  //   vec_znx_negate(module.get(),          //
  //                  MINUS_S.data(), 1, N,  //
  //                  skey.s.data(), 1, N);
  //   for (uint64_t idx = 0; idx < ONIONPIR_query2exp_nb; ++idx) {
  //     // check that we have a prgsw of colbit
  //     ASSERT_NE(qexp.query_exp_phase2[idx], nullptr);
  // #ifdef TEST_MODE
  //     uint64_t colbit = (col >> idx) & 1;
  //     int64_t(*q)[out_ncols][N] = (int64_t(*)[out_ncols][N])qexp.query_exp_phase2_raw[idx].data();
  //     for (uint64_t i = 0; i < out_nrows; ++i) {
  //       message_and_noise decr = decrypt(module.get(), K * ELLT, skey,  //
  //                                        *q[i], out_ncols);
  //       ASSERT_LE(decr.noise_amplitude, pow(2., -80));
  //       if (colbit == 0) {
  //         check_2D_zero(decr.decrypted.data(), ELLT, ONIONPIR_N);
  //       } else {
  //         for (uint64_t j = 0; j < N; ++j) {
  //           int64_t& vj = decr.decrypted[(i / 2) * N + j];
  //           if (i % 2 == 1) {
  //             ASSERT_EQ(vj, ONE[j]);
  //           } else {
  //             ASSERT_EQ(vj, MINUS_S[j]);
  //           }
  //           vj = 0;
  //         }
  //         check_2D_zero(decr.decrypted.data(), ELLT, ONIONPIR_N);
  //       }
  //     }
  // #endif
  //   }
}
