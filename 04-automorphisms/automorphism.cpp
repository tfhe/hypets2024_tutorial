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
                        const VMP_PMAT* autom_ks_a, const VMP_PMAT* autom_ks_b);
