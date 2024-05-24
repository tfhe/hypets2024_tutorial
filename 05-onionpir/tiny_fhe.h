#ifndef SPQLIOS_SAMPLES_TINY_FHE_H
#define SPQLIOS_SAMPLES_TINY_FHE_H

#include <cassert>
#include <cstdint>
#include <iostream>

#include "spqlios/arithmetic/vec_znx_arithmetic.h"

/** @brief macro that crashes if the condition are not met */
#define REQUIRE_DRAMATICALLY(req_contition, error_msg)                                                        \
  do {                                                                                                        \
    if (!(req_contition)) {                                                                                   \
      std::cerr << "REQUIREMENT FAILED at " << __FILE__ << ":" << __LINE__ << ": " << error_msg << std::endl; \
      abort();                                                                                                \
    }                                                                                                         \
  } while (0)

#include <cstdint>
#include <random>

uint8_t* get_tmp_space(const uint64_t bytes);

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

/**
 * returns the offset to the nereast exact multiple of 2^-mess_basebit
 */
void noise_values(uint64_t nn,               // ring dimension
                  uint64_t k,                // message is k-normalized
                  uint64_t message_basebit,  //
                  double* res, int64_t* a, uint64_t a_size, uint64_t a_sl);

/**
 * returns the double value of the coefficients, in [-1/2,1/2[
 */
void torus_coeffs_values(uint64_t nn,  // ring dimension
                         uint64_t k,   // message is k-normalized
                         double* res, const int64_t* a, uint64_t a_size, uint64_t a_sl);

/**
 * rounds a k-normalized input to the nearest multiple of 2^-mess_basebit
 * and returns the remainder in [-2^-(mess_basebit-1),2^-(mess_basebit-1)[
 */
void round_polynomial(uint64_t nn,               // ring dimension
                      uint64_t k,                // message is k-normalized
                      uint64_t message_basebit,  //
                      int64_t* res, uint64_t res_size, uint64_t res_sl, const int64_t* a, uint64_t a_size,
                      uint64_t a_sl);

// res = in(X^p)
// res_size must be max(in_size, ncols) in this function
// the output does not need to be fully normalized
void apply_automorphism(const MODULE* module,                      // N
                        int64_t p,                                 // power of automophism
                        uint64_t k,                                //
                        int64_t* res_rlwe, uint64_t res_size,      //
                        const int64_t* in_rlwe, uint64_t in_size,  //
                        const VMP_PMAT* autom_ks, uint64_t autom_key_nrows, uint64_t autom_key_ncols);

void rlwe_encrypt_base2k(const MODULE* module,                                    //
                         uint64_t log2_base2k,                                    // output base 2^K
                         int64_t* new_b, uint64_t new_b_size, uint64_t new_b_sl,  // b part of rlwe
                         const SVP_PPOL* s,                                       // secret key
                         const int64_t* a, uint64_t a_size, uint64_t a_sl,        // a part of rlwe
                         const int64_t* phi, uint64_t phi_size, uint64_t phi_sl   // message + noise
);

void rlwe_phase_base2k(const MODULE* module,                              //
                       uint64_t log2_base2k,                              // output base 2^K
                       int64_t* phi, uint64_t phi_size, uint64_t phi_sl,  // decrypted phase
                       const SVP_PPOL* s,                                 // secret key
                       const int64_t* a, uint64_t a_size, uint64_t a_sl,  // a part of rlwe
                       const int64_t* b, uint64_t b_size, uint64_t b_sl   // message + noise
);

// res0 = in + in(X^p)
// res1 = rotate(in - in(X^p))
// the output is normalized
// function should work perfectly in place too
void rlwe_trace_expand_1step(const MODULE* module,                      //
                             int64_t autom_p, int64_t rotate_p,         //
                             uint64_t log2_base2k,                      //
                             int64_t* res0_rlwe, uint64_t res0_size,    //
                             int64_t* res1_rlwe, uint64_t res1_size,    //
                             const int64_t* in_rlwe, uint64_t in_size,  //
                             const VMP_PMAT* autom_ks, uint64_t autom_key_nrows, uint64_t autom_key_ncols);

// res0 = in + in(X^p)
// the output is normalized
// function should work perfectly in place too
void rlwe_trace_expand_1step_no_rotation(const MODULE* module,                      //
                                         int64_t autom_p,                           //
                                         uint64_t log2_base2k,                      //
                                         int64_t* res0_rlwe, uint64_t res0_size,    //
                                         const int64_t* in_rlwe, uint64_t in_size,  //
                                         const VMP_PMAT* autom_ks, uint64_t autom_key_nrows, uint64_t autom_key_ncols);

#endif  // SPQLIOS_SAMPLES_TINY_FHE_H
