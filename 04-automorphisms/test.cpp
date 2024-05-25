#include <gtest/gtest.h>

#include "automorphisms.h"

TEST(automorphism, encrypt_decrypt) {
  MODULE* module = new_module_info(N, FFT64);
  SVP_PPOL* skey = svp_ppol_alloc(module);
  Int64VecN mu[message_limbs];
  Int64VecN a[ell];
  Int64VecN b[ell];
  Int64VecN decrypted[message_limbs];

  // generate a secret key
  std::vector<int64_t> skey_raw(N);
  random_binary(N, skey_raw.data());
  svp_prepare(module, skey, skey_raw.data());

  // generate a random message
  for (uint64_t i=0; i<message_limbs; ++i) {
    random_centered_reduced(N, K, mu[i]);
  }

  // encrypt it (noise level = K*ell)
  rlwe_encrypt(module, K, *a, *b, ell, *mu, skey);
  double noise_log2 = rlwe_decrypt(module, K, *decrypted, *a, *b, ell, skey);
  ASSERT_LE(noise_log2, -double(K*ell));
  ASSERT_TRUE(memcmp(*mu,*decrypted, message_limbs*N*sizeof(int64_t))==0);

  free(skey);
  delete_module_info(module);
}
