#ifndef SPQLIOS_ONIONPIR_H
#define SPQLIOS_ONIONPIR_H

#include <cstring>
#include <map>
#include <memory>
#include <vector>

#define NO_COPY(typenam)            \
  typenam(const typenam&) = delete; \
  void operator=(const typenam&) = delete

#include "spqlios/arithmetic/vec_znx_arithmetic.h"

template <>
struct std::default_delete<MODULE> {
  void operator()(MODULE* ptr) const { delete_module_info(ptr); }
};
template <>
struct std::default_delete<VMP_PMAT> {
  void operator()(VMP_PMAT* ptr) const { free(ptr); }
};
template <>
struct std::default_delete<SVP_PPOL> {
  void operator()(VMP_PMAT* ptr) const { free(ptr); }
};

// some unique pointers definition that autodelete objects
typedef std::unique_ptr<MODULE> MODULE_UNIPTR;
typedef std::unique_ptr<VMP_PMAT> VMP_PMAT_UNIPTR;
typedef std::unique_ptr<SVP_PPOL> SVP_PPOL_UNIPTR;

// some dimensions
#ifndef BIG_N_MODE
#define TEST_MODE
#endif  // BIG_N_MODE
#ifdef TEST_MODE
constexpr uint64_t ONIONPIR_N = 64;  // can use 64 during dev/tests

constexpr uint64_t ONIONPIR_K = 16;
// automorphism keys
constexpr uint64_t ONIONPIR_automkey_nrows = 8;
constexpr uint64_t ONIONPIR_automkey_ncols = 16;
constexpr uint64_t ONIONPIR_rk_nrows = 8;
constexpr uint64_t ONIONPIR_rk_ncols = 16;
// input query
constexpr uint64_t ONIONPIR_query1_ell = 8;
constexpr uint64_t ONIONPIR_query1_ncols = 16;
constexpr uint64_t ONIONPIR_query2_ncols = 16;
constexpr uint64_t ONIONPIR_query2_denom = 16;  // shall be >= log2(dbcols)*q2_ellt
// expanded query
constexpr uint64_t ONIONPIR_query1exp_ell = 6;
constexpr uint64_t ONIONPIR_query1exp_ellt = 4;    // plaintexts have 4 limbs
constexpr uint64_t ONIONPIR_query1exp_nrows = 64;  // N/ellt prgsw stacked
constexpr uint64_t ONIONPIR_query1exp_ncols = 12;  // 2*ell
// database
constexpr uint64_t ONIONPIR_db_nrows = 16;  // N/ellt
constexpr uint64_t ONIONPIR_db_ncols = 8;   // 16M/nrows
// results
constexpr uint64_t ONIONPIR_query2exp_nb = 3;  // log2(ONIONPIR_db_ncols)
constexpr uint64_t ONIONPIR_query2exp_ellt = 5;
constexpr uint64_t ONIONPIR_query2exp_nrows = 10;
constexpr uint64_t ONIONPIR_query2exp_ncols = 12;

constexpr uint64_t ONIONPIR_results1_size = 10;
constexpr uint64_t ONIONPIR_results2_size = 9;

#else
constexpr uint64_t ONIONPIR_N = 4096;  // can use 64 during dev/tests

constexpr uint64_t ONIONPIR_K = 16;
// automorphism keys
constexpr uint64_t ONIONPIR_automkey_nrows = 8;
constexpr uint64_t ONIONPIR_automkey_ncols = 16;
constexpr uint64_t ONIONPIR_rk_nrows = 8;
constexpr uint64_t ONIONPIR_rk_ncols = 16;
// input query
constexpr uint64_t ONIONPIR_query1_ell = 8;
constexpr uint64_t ONIONPIR_query1_ncols = 16;
constexpr uint64_t ONIONPIR_query2_ncols = 16;
constexpr uint64_t ONIONPIR_query2_denom = 128;  // shall be >= log2(dbcols)*q2_ellt
// expanded query
constexpr uint64_t ONIONPIR_query1exp_ell = 6;
constexpr uint64_t ONIONPIR_query1exp_ellt = 4;      // plaintexts have 4 limbs
constexpr uint64_t ONIONPIR_query1exp_nrows = 4096;  // N/ellt prgsw stacked
constexpr uint64_t ONIONPIR_query1exp_ncols = 12;    // 2*ell
// database
constexpr uint64_t ONIONPIR_db_nrows = 1024;   // N/ellt
constexpr uint64_t ONIONPIR_db_ncols = 16384;  // 16M/nrows
// results
constexpr uint64_t ONIONPIR_query2exp_nb = 7;  // log2(ONIONPIR_query2_denom)
constexpr uint64_t ONIONPIR_query2exp_ellt = 5;
constexpr uint64_t ONIONPIR_query2exp_nrows = 10;
constexpr uint64_t ONIONPIR_query2exp_ncols = 12;

constexpr uint64_t ONIONPIR_results1_size = 10;
constexpr uint64_t ONIONPIR_results2_size = 9;

#endif

struct onionpir_secret_key {
  // main binary key s of size N
  std::vector<int64_t> s;
  SVP_PPOL* ppol_s;

  inline onionpir_secret_key() : ppol_s(nullptr) {}
  inline ~onionpir_secret_key() { free(ppol_s); }
  NO_COPY(onionpir_secret_key);
};

struct onionpir_cloud_key {
  // 12 automorphisms for the trace expansion (phases 1 and 2)
  std::map<int64_t, VMP_PMAT_UNIPTR> autom;
  // rgsw key to multiply by s
  VMP_PMAT_UNIPTR rk_s;
#ifdef TEST_MODE
  std::map<int64_t, std::vector<int64_t>> autom_raw;
  std::vector<int64_t> rk_s_raw;
#endif

  onionpir_cloud_key() = default;
  NO_COPY(onionpir_cloud_key);
};
struct onionpir_input_query {
  // phase1 query is one rlwe (ell=8)
  std::vector<int64_t> phase1;
  // phase2 query is one rlwe (ell=8)
  std::vector<int64_t> phase2;

  onionpir_input_query() = default;
  ~onionpir_input_query() = default;
  NO_COPY(onionpir_input_query);
};
struct onionpir_expanded_query {
  // phase1 expanded query is one 4096 x 12 pmat (1024 prgsw stacked together)
  VMP_PMAT* query_exp_phase1;
  // phase2 expanded query are <=14 rgsw pmat
  VMP_PMAT* query_exp_phase2[14];
#ifdef TEST_MODE
  std::vector<int64_t> query_exp_phase1_raw;
  std::vector<int64_t> query_exp_phase2_raw[14];
#endif  // TEST_MODE

  onionpir_expanded_query() {
    query_exp_phase1 = nullptr;
    for (uint64_t i = 0; i < 14; ++i) query_exp_phase2[i] = nullptr;
  }
  ~onionpir_expanded_query() {
    free(query_exp_phase1);
    for (uint64_t i = 0; i < 14; ++i) free(query_exp_phase2[i]);
  }
  NO_COPY(onionpir_expanded_query);
};

struct onionpir_phase1_results {
  std::vector<int64_t> res;

  onionpir_phase1_results() = default;
  NO_COPY(onionpir_phase1_results);
};

struct onionpir_phase2_results {
  std::vector<int64_t> res;

  onionpir_phase2_results() = default;
  NO_COPY(onionpir_phase2_results);
};

/** encode 56 bits cleartext into a plaintext of 4x coeffs in [-2^15;2^15[ */
void onionpir_encode_plaintext(const MODULE* module,        // N
                               int64_t* encoded_plaintext,  // small int64[4*N] K-reduced
                               const uint8_t* cleartext     // bytes[7*N]
);

/** decode plaintext of 4x coeffs in [-2^15;2^15[ back to 56 bits cleartext */
void onionpir_decode_plaintext(const MODULE* module,             // N
                               uint8_t* decoded_cleartext,       // bytes[7*N]
                               const int64_t* encoded_plaintext  // small int64[4*N] K-reduced
);

/** generate the onion pir secret key */
void onionpir_generate_secret_key(const MODULE* module, onionpir_secret_key& skey);

/** generate the onion pir cloud key (phase1 only) */
void onionpir_generate_cloud_key_phase1(const MODULE* module, onionpir_cloud_key& ckey,
                                        const onionpir_secret_key& skey);

/** generate the onion pir cloud key (phase2 only) */
void onionpir_generate_cloud_key_phase2(const MODULE* module, onionpir_cloud_key& ckey,
                                        const onionpir_secret_key& skey);

/** generate the onion pir query (phase1 only) */
void onionpir_generate_query_phase1(const MODULE* module,            // N
                                    uint64_t row,                    // row index in [0,1023]
                                    onionpir_input_query& qin,       // compressed query
                                    const onionpir_secret_key& skey  // secret key
);

/** generate the onion pir query (phase2 only) */
void onionpir_generate_query_phase2(const MODULE* module,                  //
                                    uint64_t col,                          // col to be queried
                                    onionpir_input_query& res,             //
                                    const onionpir_secret_key& skey,       //
                                    uint64_t db_ncols = ONIONPIR_db_ncols  // number of DB query columns 1.. 16384
);

/** generate the expanded onion pir query from the secret key (phase1 only) */
void onionpir_generate_queryexp_phase1(const MODULE* module,            // N
                                       uint64_t row,                    // row index in [0,1023]
                                       onionpir_expanded_query& qexp,   // expanded query
                                       const onionpir_secret_key& skey  // secret key
);

/** expand the phase 1 query with the cloud key (cloud latency) */
void onionpir_query_expand_phase1(const MODULE* module,                 //
                                  onionpir_expanded_query& query_exp,   //
                                  const onionpir_cloud_key& cloud_key,  //
                                  const onionpir_input_query& input_query);

/** expand the phase 2 query with the cloud key (cloud latency) */
void onionpir_query_expand_phase2(const MODULE* module,                 //
                                  onionpir_expanded_query& query_exp,   //
                                  const onionpir_cloud_key& cloud_key,  //
                                  const onionpir_input_query& input_query);

/** generate the expanded onion pir query from the secret key (phase2 only) */
void onionpir_generate_queryexp_phase2(const MODULE* module,                  // N
                                       uint64_t col,                          // col index in DB_cols
                                       onionpir_expanded_query& qexp,         // expanded query
                                       const onionpir_secret_key& skey,       // secret key
                                       uint64_t db_ncols = ONIONPIR_db_ncols  // number of DB query columns 1.. 16384
);

/** evaluate the server computation (phase1 only) */
EXPORT void onionpir_online_phase1(const MODULE* module,  // N
                                   onionpir_phase1_results& result, const uint8_t* plaintext_db,
                                   const onionpir_expanded_query& qexp, uint64_t db_ncols = ONIONPIR_db_ncols);

/** evaluate the server computation (phase2 only) */
EXPORT void onionpir_online_phase2(const MODULE* module,                 // N
                                   onionpir_phase2_results& result,      // final result
                                   onionpir_phase1_results& input,       // result of phase 1
                                   const onionpir_expanded_query& qexp,  // expanded query
                                   uint64_t db_ncols = ONIONPIR_db_ncols);

/** final answer decryption */
EXPORT void onionpir_final_decrypt(const MODULE* module,             // N
                                   uint8_t* final_cleartext,         // final result 7*N bytes
                                   const onionpir_secret_key& skey,  // secret key
                                   const onionpir_phase2_results& c);

/** generate the expanded onion pir query from the secret key (phase2 only) */
void onionpir_rlwe_trivial_encrypt_inplace(const MODULE* module,               // N
                                           int64_t* rlwe, uint64_t rlwe_size,  // interleaved rlwe
                                           const onionpir_secret_key& key,     //
                                           bool a_is_zero = true);

// out[2i]  encodes -S * phase(in[i])
// out[2i+1]  is a (truncated) copy of in[i])
// function will be called out of place
void onionpir_prgsw_to_rgsw(const MODULE* module, const onionpir_cloud_key& ckey, int64_t* out_rgsw, uint64_t out_nrows,
                            uint64_t out_ncols, const int64_t* in_prgsw, uint64_t in_nrows, uint64_t in_ncols);

// res[i] = encodes nrows*coeff[i] of the input phase
// the output is fully normalized
void onionpir_trace_expand(const MODULE* module,                                        //
                           int64_t* res_rlwes, uint64_t res_nrows, uint64_t res_ncols,  //
                           const int64_t* in_rlwe, uint64_t in_size,                    //
                           const onionpir_cloud_key& ckey                               //
);

#endif  // SPQLIOS_ONIONPIR_H
