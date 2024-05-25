#ifndef SPQLIOS_HYPETS_SAMPLES_AUTOMORPHISMS_H
#define SPQLIOS_HYPETS_SAMPLES_AUTOMORPHISMS_H

#include <cstdlib>
#include <iostream>
#include <random>

#include "spqlios/arithmetic/vec_znx_arithmetic.h"

constexpr uint64_t N = 64;             // fixed ring dimension
constexpr uint64_t K = 16;             // fixed limb size (2^16 bits)
constexpr uint64_t message_limbs = 4;  // number of message limbs
constexpr uint64_t ell = 5;            // size of input and output ringlwe
constexpr uint64_t autom_nrows = 5;    // nrows of the autom matrix
constexpr uint64_t autom_ncols = 6;    // nrows of the autom matrix

typedef int64_t Int64VecN[N];

uint8_t* get_tmp_space(uint64_t bytes);

/**
 * decrypts the input ciphertext, and outputs the noise level detected
 * plaintext is base 2^k normalized, and  */
double rlwe_decrypt(const MODULE* module, uint64_t k,                     //
                    int64_t* mu,                                          //
                    const int64_t* a, const int64_t* b, uint64_t nlimbs,  //
                    const SVP_PPOL* skey);

void rlwe_encrypt(const MODULE* module, uint64_t k,         //
                  int64_t* a, int64_t* b, uint64_t nlimbs,  //
                  const int64_t* mu, const SVP_PPOL* skey);

void create_keyswitch(const MODULE* module, int64_t p, uint64_t k,  //
                      VMP_PMAT* autom_ks_a, VMP_PMAT* autom_ks_b,   //
                      const int64_t* skey, const SVP_PPOL* skey_dft);

void apply_automorphism(const MODULE* module, int64_t p, uint64_t k,  //
                        int64_t* res_a, int64_t* res_b,               //
                        const int64_t* in_a, const int64_t* in_b,     //
                        const VMP_PMAT* autom_ks_a, const VMP_PMAT* autom_ks_b);

void apply_automorphism_on_plaintext(const MODULE* module, int64_t p, uint64_t k,  //
                                     int64_t* res_mu,                              //
                                     const int64_t* in_mu);

/** @brief macro that crashes if the condition are not met */
#define REQUIRE_DRAMATICALLY(req_contition, error_msg)                                                        \
  do {                                                                                                        \
    if (!(req_contition)) {                                                                                   \
      std::cerr << "REQUIREMENT FAILED at " << __FILE__ << ":" << __LINE__ << ": " << error_msg << std::endl; \
      abort();                                                                                                \
    }                                                                                                         \
  } while (0)

typedef std::default_random_engine rng;
/** @brief reference to the default test rng */
rng& randgen();

/**
 * @brief      Uniformly random base2k limb of size n in
 *             [-2^(log2_base2k-1),2^(log2_base2k-1)-1]
 */
void random_centered_reduced(uint64_t n, uint64_t log2_base2k, int64_t* res);

/**
 * @brief      Uniformly random limb of size n in [-2^log2bound,2^log2bound]
 */
void random_log2bound_symmetric(uint64_t n, uint64_t log2bound, int64_t* res);

/**
 * @brief      Normal(0,stddev) noise lwe limb of size n in
 *             [-2^log2bound,2^log2bound] (right bound included)
 */
void random_normal(uint64_t n, double stddev, int64_t* res);

/**
 * @brief      Uniformly random limb of size n in {0, 1}
 */
void random_binary(uint64_t n, int64_t* res);

/**
 * rounds a k-normalized input to the nearest multiple of 2^-mess_basebit
 * and returns the remainder in [-2^-(mess_basebit-1),2^-(mess_basebit-1)[
 */
void round_polynomial_noise(uint64_t nn,              // ring dimension
                            int64_t k,                // message is k-normalized
                            int64_t message_basebit,  //
                            double* rem, int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a,
                            uint64_t a_size, uint64_t a_sl);

#define DECLARE_3D_ZERO_MATRIX(varname, vartype, vardim1, vardim2, vardim3) \
  std::vector<vartype> varname##_raw((vardim1) * (vardim2) * (vardim3), 0); \
  vartype(*varname)[vardim2][vardim3] = (vartype(*)[vardim2][vardim3])varname##_raw.data()

#endif  // SPQLIOS_HYPETS_SAMPLES_AUTOMORPHISMS_H
