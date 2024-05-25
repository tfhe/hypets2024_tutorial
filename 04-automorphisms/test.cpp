#include <gtest/gtest.h>

#include "automorphisms.h"

TEST(automorphism, encrypt_decrypt) {
  constexpr uint64_t lwe_size = autom_ncols;
  MODULE* module = new_module_info(N, FFT64);
  SVP_PPOL* skey = svp_ppol_alloc(module);
  Int64VecN mu[message_limbs];
  Int64VecN a[lwe_size];
  Int64VecN b[lwe_size];
  Int64VecN decrypted[message_limbs];

  // generate a secret key
  std::vector<int64_t> skey_raw(N);
  random_binary(N, skey_raw.data());
  svp_prepare(module, skey, skey_raw.data());

  // generate a random message
  for (uint64_t i = 0; i < message_limbs; ++i) {
    random_centered_reduced(N, K, mu[i]);
  }

  // encrypt it (noise level = K*ell)
  rlwe_encrypt(module, K, *a, *b, lwe_size, *mu, message_limbs, skey);
  double noise_log2 = rlwe_decrypt(module, K, *decrypted, message_limbs, *a, *b, lwe_size, skey);
  ASSERT_LE(noise_log2, -double(K * lwe_size));
  ASSERT_TRUE(memcmp(*mu, *decrypted, message_limbs * N * sizeof(int64_t)) == 0);

  free(skey);
  delete_module_info(module);
}

TEST(automorphism, automorphism) {
  MODULE* module = new_module_info(N, FFT64);

  Int64VecN skey_raw;
  SVP_PPOL* skey = svp_ppol_alloc(module);

  Int64VecN mu[message_limbs];

  int64_t p = 3;

  VMP_PMAT* ks_a = vmp_pmat_alloc(module, autom_nrows, autom_ncols);
  VMP_PMAT* ks_b = vmp_pmat_alloc(module, autom_nrows, autom_ncols);

  Int64VecN a[ell];
  Int64VecN b[ell];

  Int64VecN autom_a[ell];
  Int64VecN autom_b[ell];

  Int64VecN decrypted[message_limbs];

  Int64VecN autom_mu[message_limbs];

  // generate a secret key
  random_binary(N, skey_raw);
  svp_prepare(module, skey, skey_raw);

  // generate the autom ks
  create_keyswitch(module, p, K,  //
                   ks_a, ks_b,    //
                   skey_raw, skey);

  // generate a random message
  for (uint64_t i = 0; i < message_limbs; ++i) {
    random_centered_reduced(N, K, mu[i]);
  }

  // encrypt it (noise level = K*ell)
  rlwe_encrypt(module, K, *a, *b, ell, *mu, message_limbs, skey);

  // apply automorphism on ciphertext
  apply_automorphism(module, p, K,        //
                     *autom_a, *autom_b,  //
                     *a, *b,              //
                     ks_a, ks_b);         //

  // apply the automorphism in plaintext
  apply_automorphism_on_plaintext(module, p, K, *autom_mu, *mu);

  double noise_log2 = rlwe_decrypt(module, K, *decrypted, message_limbs, *autom_a, *autom_b, ell, skey);
  // ASSERT_LE(noise_log2, -double(K * ell));
  for (uint64_t i = 0; i < message_limbs; ++i) {
    for (uint64_t j = 0; j < N; j++) {
      ASSERT_EQ(autom_mu[i][j], decrypted[i][j]);
    }
  }

  free(skey);
  free(ks_b);
  free(ks_a);
  delete_module_info(module);
}
