#include "onionpir.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "spqlios/arithmetic/vec_znx_arithmetic.h"
#include "tiny_fhe.h"

// magic macros to define (large) 2d or 3d matrices initialized to zero
// when the last dimensions are static constant only
// the matrix dissapears at the end of the block

// for all practical purposes, it is equivalent to:
// vartype varname[vardim1][vardim2];
// but less prone to stack overflows...

#define DECLARE_2D_ZERO_MATRIX(varname, vartype, vardim1, vardim2) \
  std::vector<vartype> varname##_raw((vardim1) * (vardim2), 0);    \
  vartype(*varname)[vardim2] = (vartype(*)[vardim2])varname##_raw.data()

#define DECLARE_3D_ZERO_MATRIX(varname, vartype, vardim1, vardim2, vardim3) \
  std::vector<vartype> varname##_raw((vardim1) * (vardim2) * (vardim3), 0); \
  vartype(*varname)[vardim2][vardim3] = (vartype(*)[vardim2][vardim3])varname##_raw.data()

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

void onionpir_encode_plaintext(const MODULE* module,        // N
                               int64_t* encoded_plaintext,  // small int64[4*N] K-reduced
                               const uint8_t* cleartext     // bytes[7*N]
) {
  const uint64_t nn = ONIONPIR_N;

  // fill least-significant znx with payload
  for (uint64_t i = 0, k = 0; i < nn; ++i, k += 7) {
    encoded_plaintext[3 * nn + i] = *(int8_t*)(cleartext + k) << 8;
    encoded_plaintext[2 * nn + i] = *(int16_t*)(cleartext + k + 1);
    encoded_plaintext[nn + i] = *(int16_t*)(cleartext + k + 3);
    encoded_plaintext[i] = *(int16_t*)(cleartext + k + 5);
  }
}

void onionpir_decode_plaintext(const MODULE* module,             // N
                               uint8_t* decoded_cleartext,       // bytes[7*N]
                               const int64_t* encoded_plaintext  // small int64[4*N] K-reduced
) {
  const uint64_t nn = ONIONPIR_N;

  for (uint64_t i = 0, k = 0; i < nn; ++i, k += 7) {
    *(int8_t*)(decoded_cleartext + k) = (int16_t)encoded_plaintext[3 * nn + i] >> 8;
    *(int16_t*)(decoded_cleartext + k + 1) = (int16_t)encoded_plaintext[2 * nn + i];
    *(int16_t*)(decoded_cleartext + k + 3) = (int16_t)encoded_plaintext[nn + i];
    *(int16_t*)(decoded_cleartext + k + 5) = (int16_t)encoded_plaintext[i];
  }
}

/** encrypt a trivial (that already contains noise)  */
void onionpir_rlwe_trivial_encrypt_inplace(const MODULE* module,               // N
                                           int64_t* rlwe, uint64_t rlwe_size,  // interleaved rlwe
                                           const onionpir_secret_key& key,     //
                                           bool a_is_zero) {
  const uint64_t a_size = (rlwe_size + 1) / 2;
  const uint64_t b_size = rlwe_size / 2;
  int64_t(*c)[ONIONPIR_N] = (int64_t(*)[ONIONPIR_N])rlwe;
  std::vector<int64_t> save_a;
  // randomize a
  if (!a_is_zero) {
    save_a.resize(a_size * ONIONPIR_N);
    // vec_znx_copy(module, //
    //              save_a.data(), a_size, ONIONPIR_N,//
    //              *c, a_size, 2*ONIONPIR_N);
    for (uint64_t i = 0; i < a_size; ++i) {
      memcpy(save_a.data() + i * ONIONPIR_N, c[2 * i], ONIONPIR_N * sizeof(double));
    }
  }
  for (uint64_t i = 0; i < a_size; ++i) {
    random_centered_reduced(ONIONPIR_N, ONIONPIR_K, c[2 * i]);
  }
  // lock the ciphertext
  rlwe_encrypt_base2k(module, ONIONPIR_K,             //
                      c[1], b_size, 2 * ONIONPIR_N,   // b
                      key.ppol_s,                     // key
                      c[0], a_size, 2 * ONIONPIR_N,   // a
                      c[1], b_size, 2 * ONIONPIR_N);  // b
  // readd_a_save
  if (!a_is_zero) {
    vec_znx_add(module,                        //
                c[0], a_size, 2 * ONIONPIR_N,  //
                c[0], a_size, 2 * ONIONPIR_N,  //
                save_a.data(), a_size, ONIONPIR_N);
    vec_znx_normalize_base2k(module, ONIONPIR_K,            //
                             c[0], a_size, 2 * ONIONPIR_N,  //
                             c[0], a_size, 2 * ONIONPIR_N,  //
                             get_tmp_space(vec_znx_big_normalize_base2k_tmp_bytes(module)));
  }
}

/** decrypts the b part  */
void onionpir_rlwe_decrypt_inplace(const MODULE* module,               // N
                                   int64_t* rlwe, uint64_t rlwe_size,  // interleaved rlwe
                                   const onionpir_secret_key& key) {
  const uint64_t a_size = (rlwe_size + 1) / 2;
  const uint64_t b_size = rlwe_size / 2;
  int64_t(*c)[ONIONPIR_N] = (int64_t(*)[ONIONPIR_N])rlwe;
  // unlock the ciphertext
  rlwe_phase_base2k(module, ONIONPIR_K,            //
                    c[1], b_size, 2 * ONIONPIR_N,  // b
                    key.ppol_s,                    // key
                    c[0], a_size, 2 * ONIONPIR_N,  // a
                    c[1], b_size, 2 * ONIONPIR_N   // b
  );
}

void onionpir_generate_secret_key(const MODULE* module, onionpir_secret_key& skey) {
  skey.s.resize(ONIONPIR_N);
  random_binary(ONIONPIR_N, skey.s.data());
  skey.ppol_s = svp_ppol_alloc(module);
  svp_prepare(module, skey.ppol_s, skey.s.data());
}

/** shortcut to generate the phase1 input query directly */
void onionpir_generate_query_phase1(const MODULE* module,            // N
                                    uint64_t row,                    // row index in [0,1023]
                                    onionpir_input_query& qin,       // compressed query
                                    const onionpir_secret_key& skey  // secret key
) {
  const uint64_t nn = ONIONPIR_N;
  DECLARE_2D_ZERO_MATRIX(q, int64_t, ONIONPIR_query1_ncols, ONIONPIR_N);
  // the trivial message encoded by the query is:
  int64_t _2K_OVER_N = (1L << ONIONPIR_K) / nn;
  q[3][4 * row] = _2K_OVER_N;      // b1
  q[5][4 * row + 1] = _2K_OVER_N;  // b2
  q[7][4 * row + 2] = _2K_OVER_N;  // b3
  q[9][4 * row + 3] = _2K_OVER_N;  // b4
  // add small noise (last limb b7 at ONIONPIR_query1exp_ncols-1)
  // TODO vary the amount of noise
  random_log2bound_symmetric(nn, 8, q[ONIONPIR_query1_ncols - 1]);
  // generate encrypted b
  onionpir_rlwe_trivial_encrypt_inplace(module, q[0], ONIONPIR_query1_ncols, skey);
  qin.phase1 = std::move(q_raw);
}

/** shortcut function to generate the phase1 expanded query directly */
void onionpir_generate_queryexp_phase1(const MODULE* module,            // N
                                       uint64_t row,                    // row index in [0,1023]
                                       onionpir_expanded_query& qexp,   // expanded query
                                       const onionpir_secret_key& skey  // secret key
) {
  const uint64_t nn = ONIONPIR_N;
  DECLARE_3D_ZERO_MATRIX(q, int64_t, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols, ONIONPIR_N);
  // the trivial message encoded by the query is:
  q[4 * row][1][0] = 1;      // b0
  q[4 * row + 1][3][0] = 1;  // b1
  q[4 * row + 2][5][0] = 1;  // b2
  q[4 * row + 3][7][0] = 1;  // b3
  // add small noise (last limb b5 at ONIONPIR_query1exp_ncols-1)
  // TODO vary the amount of noise
  const uint64_t b_size = ONIONPIR_query1exp_ncols / 2;
  for (uint64_t r = 0; r < ONIONPIR_query1exp_nrows; ++r) {
    random_log2bound_symmetric(nn, 2, q[r][2 * b_size - 1]);
  }
  // generate encrypted b
  for (uint64_t r = 0; r < ONIONPIR_query1exp_nrows; ++r) {
    onionpir_rlwe_trivial_encrypt_inplace(module, *q[r], ONIONPIR_query1exp_ncols, skey);
  }
  // generate the pmat
  qexp.query_exp_phase1 = vmp_pmat_alloc(module, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols);
  VMP_PMAT* qexpptr = qexp.query_exp_phase1;
  uint8_t* tmp_space =
      get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols));
  vmp_prepare_contiguous(module, qexpptr, **q, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols, tmp_space);
#ifdef TEST_MODE
  qexp.query_exp_phase1_raw = std::move(q_raw);
#endif
}

void onionpir_generate_cloud_key_phase1(const MODULE* module, onionpir_cloud_key& ckey,
                                        const onionpir_secret_key& skey) {
  ckey.autom.clear();
#ifdef TEST_MODE
  ckey.autom_raw.clear();
#endif
  uint8_t* tmp_space =
      get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols));
  int64_t pm = 4096;
  for (uint64_t i = 0; i < 12; ++i) {
    int64_t p = pm + 1;
    std::vector<int64_t> autom_s(ONIONPIR_N);
    for (uint64_t j = 0; j < ONIONPIR_N; ++j) {
      autom_s[j] = -skey.s[j];
    }
    vec_znx_automorphism(module, p,                      //
                         autom_s.data(), 1, ONIONPIR_N,  //
                         autom_s.data(), 1, ONIONPIR_N);
    DECLARE_3D_ZERO_MATRIX(q, int64_t, ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols, ONIONPIR_N);
    uint64_t b_size = ONIONPIR_automkey_ncols / 2;
    for (uint64_t j = 0; j < ONIONPIR_automkey_nrows; ++j) {
      memcpy(q[j][2 * j + 1], autom_s.data(), ONIONPIR_N * sizeof(int64_t));
      random_log2bound_symmetric(ONIONPIR_N, 8, q[j][2 * b_size - 1]);
      // onionpir_rlwe_trivial_encrypt_inplace(module, *q[j], ONIONPIR_automkey_ncols, skey);
    }
    VMP_PMAT_UNIPTR autom_pmat(vmp_pmat_alloc(module, ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols));
    vmp_prepare_contiguous(module, autom_pmat.get(), **q, ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols, tmp_space);
    ckey.autom[p] = std::move(autom_pmat);
#ifdef TEST_MODE
    ckey.autom_raw[p] = std::move(q_raw);
#endif
    pm >>= 1;
  }
}

void onionpir_generate_cloud_key_phase2(const MODULE* module, onionpir_cloud_key& ckey,
                                        const onionpir_secret_key& skey) {
  uint8_t* tmp_space = get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, ONIONPIR_rk_nrows, ONIONPIR_rk_ncols));
  DECLARE_3D_ZERO_MATRIX(q, int64_t, ONIONPIR_rk_nrows, ONIONPIR_rk_ncols, ONIONPIR_N);
  std::vector<int64_t> minus_s(ONIONPIR_N);
  for (uint64_t i = 0; i < ONIONPIR_N; ++i) {
    minus_s[i] = -skey.s[i];
  }
  const uint64_t b_size = ONIONPIR_automkey_ncols / 2;
  for (uint64_t j = 0; j < ONIONPIR_automkey_nrows; ++j) {
    memcpy(q[j][2 * j], minus_s.data(), ONIONPIR_N * sizeof(int64_t));
    random_log2bound_symmetric(ONIONPIR_N, 8, q[j][2 * b_size - 1]);
    onionpir_rlwe_trivial_encrypt_inplace(module, *q[j], ONIONPIR_rk_ncols, skey, false);
  }
  VMP_PMAT_UNIPTR rk_pmat(vmp_pmat_alloc(module, ONIONPIR_rk_nrows, ONIONPIR_rk_ncols));
  vmp_prepare_contiguous(module, rk_pmat.get(), **q, ONIONPIR_rk_nrows, ONIONPIR_rk_ncols, tmp_space);
  ckey.rk_s = std::move(rk_pmat);
#ifdef TEST_MODE
  ckey.rk_s_raw = std::move(q_raw);
#endif
}

EXPORT void onionpir_online_phase1(const MODULE* module,  // N
                                   onionpir_phase1_results& result, const uint8_t* plaintext_db,
                                   const onionpir_expanded_query& qexp, uint64_t db_ncols) {
  const uint64_t result_a_size = (ONIONPIR_results1_size + 1) / 2;
  const uint64_t result_b_size = (ONIONPIR_results1_size) / 2;
  uint64_t tmp_bytes_size = vmp_apply_dft_tmp_bytes(  //
      module,                                         //
      ONIONPIR_query1exp_ncols,                       // res
      ONIONPIR_query1exp_nrows,                       // in
      ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols);
  uint64_t t2 = vec_znx_big_range_normalize_base2k_tmp_bytes(module);
  if (t2 > tmp_bytes_size) tmp_bytes_size = t2;
  uint8_t* tmp_bytes = get_tmp_space(tmp_bytes_size);  // prealloc
  VEC_ZNX_DFT* col_res = vec_znx_dft_alloc(module, ONIONPIR_query1exp_ncols);
  VEC_ZNX_BIG* col_res_big = (VEC_ZNX_BIG*)col_res;
  std::vector<int64_t> encoded_plaintext(ONIONPIR_db_nrows * ONIONPIR_N * 4);
  result.res.resize(ONIONPIR_results1_size * db_ncols * ONIONPIR_N);
  const uint8_t* src = plaintext_db;  // stream pointer

  for (uint64_t db_col = 0; db_col < db_ncols; ++db_col) {
    // encode the plaintext
    int64_t* dest = encoded_plaintext.data();
    for (uint64_t i = 0; i < ONIONPIR_db_nrows; ++i) {
      onionpir_encode_plaintext(module, dest, src);
      dest += 4 * ONIONPIR_N;
      src += 7 * ONIONPIR_N;  // can make it mod something in test mode
    }
    // apply the dft
    vmp_apply_dft(module,                                                                     //
                  col_res, ONIONPIR_query1exp_ncols,                                          // res
                  encoded_plaintext.data(), ONIONPIR_query1exp_nrows, ONIONPIR_N,             // in
                  qexp.query_exp_phase1, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols,  //
                  tmp_bytes);
    // idft the result
    vec_znx_idft_tmp_a(module,                                 //
                       col_res_big, ONIONPIR_query1exp_ncols,  //
                       col_res, ONIONPIR_query1exp_ncols);
    // big-normalize to destination
    int64_t* dest_ptr = result.res.data() + db_col * ONIONPIR_results1_size * ONIONPIR_N;
    vec_znx_big_range_normalize_base2k(               // a part
        module, ONIONPIR_K,                           //
        dest_ptr, result_a_size, 2 * ONIONPIR_N,      //
        col_res_big, 0, ONIONPIR_query1exp_ncols, 2,  //
        tmp_bytes);
    vec_znx_big_range_normalize_base2k(                        // b part
        module, ONIONPIR_K,                                    //
        dest_ptr + ONIONPIR_N, result_b_size, 2 * ONIONPIR_N,  //
        col_res_big, 1, ONIONPIR_query1exp_ncols, 2,           // odd pos
        tmp_bytes);
  }
  free(col_res);
}

void onionpir_prgsw_to_rgsw(const MODULE* module, const onionpir_cloud_key& ckey, int64_t* out_rgsw, uint64_t out_nrows,
                            uint64_t out_ncols, const int64_t* in_prgsw, uint64_t in_nrows, uint64_t in_ncols) {
  REQUIRE_DRAMATICALLY(out_nrows == 2 * in_nrows, "nrows NOT_SUPPORTED");
  REQUIRE_DRAMATICALLY(out_rgsw != in_prgsw, "in place NOT_SUPPORTED");
  int64_t(*out)[ONIONPIR_N] = (int64_t(*)[ONIONPIR_N])out_rgsw;
  const int64_t(*in)[ONIONPIR_N] = (const int64_t(*)[ONIONPIR_N])in_prgsw;

  const uint64_t in_a_size = (in_ncols + 1) / 2;
  const uint64_t in_b_size = in_ncols / 2;
  const uint64_t out_a_size = (out_ncols + 1) / 2;
  const uint64_t out_b_size = out_ncols / 2;
  std::vector<int64_t> tmp_a(ONIONPIR_N * in_a_size);
  VEC_ZNX_DFT* tmp_as2 = vec_znx_dft_alloc(module, ONIONPIR_rk_ncols);
  VEC_ZNX_BIG* tmp_as2_big = (VEC_ZNX_BIG*)tmp_as2;
  uint64_t tb = 0;
  tb |= vmp_apply_dft_tmp_bytes(module,             //
                                ONIONPIR_rk_ncols,  //
                                in_a_size,          //
                                ONIONPIR_rk_nrows, ONIONPIR_rk_ncols);
  tb |= vec_znx_big_range_normalize_base2k_tmp_bytes(module);
  uint8_t* tmp_bytes = get_tmp_space(tb);

  for (uint64_t i = 0; i < in_nrows; ++i) {
    vec_znx_copy(module,                                               //
                 out[(2 * i + 1) * out_ncols], out_ncols, ONIONPIR_N,  //
                 in[(i)*in_ncols], in_ncols, ONIONPIR_N);
    // tmp_a = extract the a part of in[i] (it must be contiguous)
    vec_znx_copy(module,                               //
                 tmp_a.data(), in_a_size, ONIONPIR_N,  //
                 in[i * in_ncols], in_a_size, 2 * ONIONPIR_N);
    // tmp_as2 = vmp_apply(tmp_a, ckey.rk_s) -- apply idft
    vmp_apply_dft(module,                               //
                  tmp_as2, ONIONPIR_rk_ncols,           //
                  tmp_a.data(), in_a_size, ONIONPIR_N,  //
                  ckey.rk_s.get(), ONIONPIR_rk_nrows, ONIONPIR_rk_ncols, tmp_bytes);
    vec_znx_idft_tmp_a(module,                          //
                       tmp_as2_big, ONIONPIR_rk_ncols,  //
                       tmp_as2, ONIONPIR_rk_ncols);
    // out[2i] =  big_normalize(tmp_as2)
    vec_znx_big_range_normalize_base2k(                        //
        module, ONIONPIR_K,                                    //
        out[(2 * i) * out_ncols], out_a_size, 2 * ONIONPIR_N,  // a part
        tmp_as2_big, 0, ONIONPIR_rk_ncols, 2,                  // even range
        tmp_bytes);
    vec_znx_big_range_normalize_base2k(                                     //
        module, ONIONPIR_K,                                                 //
        out[(2 * i) * out_ncols] + ONIONPIR_N, out_b_size, 2 * ONIONPIR_N,  // b part
        tmp_as2_big, 1, ONIONPIR_rk_ncols, 2,                               // odd range
        tmp_bytes);
    // add the b part of in[i] to the a part of out[2i]
    vec_znx_add(module,                                                     //
                out[(2 * i) * out_ncols], out_a_size, 2 * ONIONPIR_N,       // out_a
                out[(2 * i) * out_ncols], out_a_size, 2 * ONIONPIR_N,       // out_a
                in[i * in_ncols] + ONIONPIR_N, in_b_size, 2 * ONIONPIR_N);  // in_b
    // and renormalize out[2i] - a part only.
    vec_znx_normalize_base2k(module, ONIONPIR_K,                                    //
                             out[(2 * i) * out_ncols], out_a_size, 2 * ONIONPIR_N,  // out_a
                             out[(2 * i) * out_ncols], out_a_size, 2 * ONIONPIR_N,  // out_a
                             tmp_bytes);
  }
  free(tmp_as2);
}

/** cmux. input size is ONIONPIR_result1_size */
void onionpir_cmux_eval(const MODULE* module, int64_t* out_rlwe, uint64_t out_size, const int64_t* c1_rlwe,
                        const int64_t* c0_rlwe, const VMP_PMAT* rgsw) {
  static constexpr uint64_t Sin = ONIONPIR_query2exp_nrows;
  static constexpr uint64_t Sout = ONIONPIR_query2exp_ncols;
  std::vector<int64_t> tmp64(ONIONPIR_N * Sout);
  VEC_ZNX_DFT* tmp = vec_znx_dft_alloc(module, Sout);
  VEC_ZNX_BIG* tmp_big = (VEC_ZNX_BIG*)tmp;
  uint64_t tb = 0;
  tb |= vmp_apply_dft_tmp_bytes(module,  //
                                Sout,    //
                                Sin,     //
                                Sin, Sout);
  tb |= vec_znx_big_range_normalize_base2k_tmp_bytes(module);
  uint8_t* tmp_bytes = get_tmp_space(tb);

  // diff = c1 - c0
  vec_znx_sub(module,                         //
              tmp64.data(), Sin, ONIONPIR_N,  //
              c1_rlwe, Sin, ONIONPIR_N,       //
              c0_rlwe, Sin, ONIONPIR_N);
  // tmp = vmp_apply(diff, pmat)
  vmp_apply_dft(module,                         //
                tmp, Sout,                      //
                tmp64.data(), Sin, ONIONPIR_N,  //
                rgsw, Sin, Sout,                //
                tmp_bytes);
  vec_znx_idft_tmp_a(module,         //
                     tmp_big, Sout,  //
                     tmp, Sout);     //

  // add c0 (in big space)
  vec_znx_big_add_small(module,         //
                        tmp_big, Sout,  //
                        tmp_big, Sout,  //
                        c0_rlwe, Sin, ONIONPIR_N);

  // normalize
  uint64_t out_a_size = (Sin + 1) / 2;
  uint64_t out_b_size = Sin / 2;
  vec_znx_big_range_normalize_base2k(module, ONIONPIR_K,                    // a
                                     out_rlwe, out_a_size, 2 * ONIONPIR_N,  //
                                     tmp_big, 0, Sout, 2,                   //
                                     tmp_bytes);
  vec_znx_big_range_normalize_base2k(module, ONIONPIR_K,                                 // b
                                     out_rlwe + ONIONPIR_N, out_b_size, 2 * ONIONPIR_N,  //
                                     tmp_big, 1, Sout, 2,                                //
                                     tmp_bytes);
  free(tmp);
}

void onionpir_generate_query_phase2(const MODULE* module,             //
                                    uint64_t col,                     //
                                    onionpir_input_query& res,        //
                                    const onionpir_secret_key& skey,  //
                                    uint64_t db_ncols                 // number of DB query columns 1.. 16384
) {
  REQUIRE_DRAMATICALLY(db_ncols > 0 and db_ncols <= ONIONPIR_db_ncols, "db_ncols must be in [1,16384]");
  REQUIRE_DRAMATICALLY(col < db_ncols, "col must be in [0,db_ncols)");

  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t ELLT = ONIONPIR_query2exp_ellt;
  constexpr uint64_t out_size = ONIONPIR_query2_ncols;
  constexpr uint64_t b_size = ONIONPIR_query2_ncols / 2;
  const uint64_t query_col_nbits = uint64_t(ceil(log2(db_ncols)));

  // get a pointer to the out rlwe
  std::vector<int64_t>& r = res.phase2;
  r.resize(ONIONPIR_query2_ncols * ONIONPIR_N, 0);
  int64_t(*q)[N] = (int64_t(*)[N])r.data();
  // set the message
  const int64_t POW2_K_OVER_DENOM = (1L << ONIONPIR_K) / ONIONPIR_query2_denom;
  for (uint64_t k = 0; k < query_col_nbits; ++k) {
    if ((col >> k) & 1) {
      for (uint64_t i = 0; i < ELLT; ++i) {
        q[2 * (i + 1) + 1][k * ELLT + i] = POW2_K_OVER_DENOM;
      }
    }
  }
  // add noise
  random_log2bound_symmetric(N, 1, q[2 * b_size - 1]);  // TODO noise levels
  // encrypt
  onionpir_rlwe_trivial_encrypt_inplace(module, *q, out_size, skey);
}

void onionpir_generate_queryexp_phase2(const MODULE* module,             //
                                       uint64_t col,                     //
                                       onionpir_expanded_query& res,     //
                                       const onionpir_secret_key& skey,  //
                                       uint64_t db_ncols                 // number of DB query columns 1.. 16384
) {
  REQUIRE_DRAMATICALLY(db_ncols > 0 and db_ncols <= ONIONPIR_db_ncols, "db_ncols must be in [1,16384]");
  REQUIRE_DRAMATICALLY(col < db_ncols, "col must be in [0,db_ncols)");

  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t out_nrows = ONIONPIR_query2exp_nrows;
  constexpr uint64_t out_ncols = ONIONPIR_query2exp_ncols;
  constexpr uint64_t out_b_size = out_ncols / 2;
  const uint64_t query_col_nbits = uint64_t(ceil(log2(db_ncols)));

  uint64_t tb = vmp_prepare_contiguous_tmp_bytes(module, out_nrows, out_ncols);
  uint8_t* tmp_bytes = get_tmp_space(tb);

  // encode all positions of col
  for (uint64_t i = 0, _2i = 1; i < query_col_nbits; ++i, _2i <<= 1) {
    int64_t colbit = (col & _2i) ? 1 : 0;
    // generate a prgsw(colbit)
    DECLARE_3D_ZERO_MATRIX(q, int64_t, out_nrows, out_ncols, N);
    for (uint64_t j = 0; j < out_nrows; ++j) {
      q[j][j][0] = colbit;                                         // message
      random_log2bound_symmetric(N, 3, q[j][2 * out_b_size - 1]);  // noise
      onionpir_rlwe_trivial_encrypt_inplace(module, *q[j], out_ncols, skey, false);
    }
    // save it
    VMP_PMAT* pmat = vmp_pmat_alloc(module, out_nrows, out_ncols);
    vmp_prepare_contiguous(module, pmat, **q, out_nrows, out_ncols, tmp_bytes);
    res.query_exp_phase2[i] = pmat;
#ifdef TEST_MODE
    res.query_exp_phase2_raw[i] = std::move(q_raw);
#endif
  }
}

EXPORT void onionpir_online_phase2(const MODULE* module,                 // N
                                   onionpir_phase2_results& result,      // final result
                                   onionpir_phase1_results& input,       // result of phase 1
                                   const onionpir_expanded_query& qexp,  // expanded query
                                   uint64_t db_ncols) {
  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t RLWE_size = ONIONPIR_query2exp_nrows;
  uint64_t current_dim = db_ncols;
  DECLARE_3D_ZERO_MATRIX(tmp, int64_t, ((db_ncols + 1) / 2), RLWE_size, N);
  int64_t(*in)[RLWE_size][N] = (int64_t(*)[RLWE_size][N])input.res.data();
  uint64_t pmax = ceil(log2(db_ncols));
  for (uint64_t p = 0; p < pmax; p++) {
    uint64_t next_dim = current_dim / 2;
    for (uint64_t i = 0; i < next_dim; i++) {
      onionpir_cmux_eval(module,                      //
                         *tmp[i], RLWE_size,          //
                         *in[2 * i + 1], *in[2 * i],  //
                         qexp.query_exp_phase2[p]);
    }
    // deal with the final position if dim is odd
    if (current_dim % 2 == 1) {
      vec_znx_copy(module, *tmp[next_dim], RLWE_size, N, *in[current_dim - 1], RLWE_size, N);
      ++next_dim;
    }
    current_dim = next_dim;
    in = tmp;  // all next iters
  }
  // truncate and copy the final result
  REQUIRE_DRAMATICALLY(current_dim == 1, "bug");
  result.res.resize(ONIONPIR_results2_size * N);
  vec_znx_copy(module,                                        //
               result.res.data(), ONIONPIR_results2_size, N,  //
               *tmp[0], RLWE_size, N);
}

EXPORT void onionpir_final_decrypt(const MODULE* module,             // N
                                   uint8_t* final_cleartext,         // final result 7*N bytes
                                   const onionpir_secret_key& skey,  // secret key
                                   const onionpir_phase2_results& c) {
  constexpr uint64_t N = ONIONPIR_N;
  constexpr uint64_t K = ONIONPIR_K;
  constexpr uint64_t RLWE_size = ONIONPIR_results2_size;
  constexpr uint64_t a_size = (RLWE_size + 1) / 2;
  constexpr uint64_t b_size = RLWE_size / 2;
  constexpr uint64_t message_basebit = 56;
  constexpr uint64_t final_size = (message_basebit + K - 1) / K;
  std::vector<int64_t> decrypted(a_size * ONIONPIR_N);
  const int64_t* rlwe = c.res.data();
  rlwe_phase_base2k((MODULE*)module, K,           //
                    decrypted.data(), a_size, N,  //
                    skey.ppol_s,                  //
                    rlwe, a_size, 2 * N,          //
                    rlwe + N, b_size, 2 * N       //
  );
  round_polynomial_noise(ONIONPIR_N, ONIONPIR_K, message_basebit,  //
                         nullptr,                                  //
                         decrypted.data(), final_size, N,          //
                         decrypted.data(), a_size, N);
  onionpir_decode_plaintext(module, final_cleartext, decrypted.data());
}

void onionpir_trace_expand(const MODULE* module,                                        //
                           int64_t* res_rlwes, uint64_t res_nrows, uint64_t res_ncols,  //
                           const int64_t* in_rlwe, uint64_t in_size,                    //
                           const onionpir_cloud_key& ckey                               //
) {
  const uint64_t nn = module_get_n(module);
  REQUIRE_DRAMATICALLY(res_nrows > 1, "res_nrows must be larger than 2");
  REQUIRE_DRAMATICALLY(res_nrows <= nn, "res_nrows must be smaller than ring size");

  // res[i] is the i-th output rlwe
  std::vector<int64_t*> res(res_nrows);
  for (uint64_t i = 0; i < res_nrows; ++i) {
    res[i] = res_rlwes + i * nn * res_ncols;
  }

  if (res_nrows == 2) {
    // 1 step, expand into output directly
    rlwe_trace_expand_1step(module, nn + 1, -1, ONIONPIR_K, res.at(0), res_ncols, res.at(1), res_ncols, in_rlwe,
                            in_size, ckey.autom.at(nn + 1).get(), ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols);
    return;
  }

  // largest power of 2 smaller than res_nrows
  const uint64_t working_nrows = 1ul << (uint64_t)std::floor(std::log2(res_nrows - 1));
  uint64_t working_ncols = std::max(in_size, ONIONPIR_automkey_ncols);

  // tmp[i] is the i-th temporary rlwe of size working_ncols
  std::vector<int64_t> tmp_space(working_nrows * working_ncols * nn);
  std::vector<int64_t*> tmp(working_nrows);
  for (uint64_t i = 0; i < working_nrows; ++i) {
    tmp[i] = tmp_space.data() + i * nn * working_ncols;
  }

  // first step: X -> X^(N+1) and rotate by -1
  rlwe_trace_expand_1step(module, nn + 1, -1, ONIONPIR_K, tmp.at(0), working_ncols, tmp.at(1), working_ncols, in_rlwe,
                          in_size, ckey.autom.at(nn + 1).get(), ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols);

  for (uint64_t p = 2; p < working_nrows; p *= 2) {
    const int64_t rotate_p = -p;
    const int64_t autom_p = nn / p + 1;

    for (uint64_t k = 0; k < p; ++k) {
      rlwe_trace_expand_1step(module, autom_p, rotate_p, ONIONPIR_K, tmp.at(k), working_ncols, tmp.at(k + p),
                              working_ncols, tmp.at(k), working_ncols, ckey.autom.at(autom_p).get(),
                              ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols);
    }
  }

  // last step output res_nrows rlwe with res_ncols limbs only
  uint64_t p = working_nrows;
  const int64_t rotate_p = -p;
  const int64_t autom_p = nn / p + 1;

  // first (res_nrows - working_nrows) of tmp will generate 2 output rlwe
  // the rest of tmp will generate 1 output rlwe only
  uint64_t k = 0;
  for (; k < res_nrows - working_nrows; ++k) {
    rlwe_trace_expand_1step(module, autom_p, rotate_p, ONIONPIR_K, res.at(k), res_ncols, res.at(k + p), res_ncols,
                            tmp.at(k), working_ncols, ckey.autom.at(autom_p).get(), ONIONPIR_automkey_nrows,
                            ONIONPIR_automkey_ncols);
  }

  for (; k < working_nrows; ++k) {
    rlwe_trace_expand_1step_no_rotation(module, autom_p, ONIONPIR_K, res.at(k), res_ncols, tmp.at(k), working_ncols,
                                        ckey.autom.at(autom_p).get(), ONIONPIR_automkey_nrows, ONIONPIR_automkey_ncols);
  }
}

void onionpir_query_expand_phase1(const MODULE* module,                //
                                  onionpir_expanded_query& query_exp,  //
                                  const onionpir_cloud_key& ckey,      //
                                  const onionpir_input_query& input_query) {
  const uint64_t nn = module_get_n(module);

#ifdef TEST_MODE
  query_exp.query_exp_phase1_raw.clear();
#endif

  // extract each query coefficient as RLWE
  std::vector<int64_t> q_rlwes_raw(ONIONPIR_query1exp_nrows * ONIONPIR_query1exp_ncols * nn);
  int64_t* q_rlwes = q_rlwes_raw.data();
  onionpir_trace_expand(module, q_rlwes, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols, input_query.phase1.data(),
                        ONIONPIR_query1_ncols, ckey);

  // to PRGSW
  if (query_exp.query_exp_phase1 != nullptr) free(query_exp.query_exp_phase1);
  query_exp.query_exp_phase1 = vmp_pmat_alloc(module, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols);

  uint8_t* tmp_space =
      get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, ONIONPIR_query1exp_nrows, ONIONPIR_query1exp_ncols));
  vmp_prepare_contiguous(module, query_exp.query_exp_phase1, q_rlwes, ONIONPIR_query1exp_nrows,
                         ONIONPIR_query1exp_ncols, tmp_space);

#ifdef TEST_MODE
  query_exp.query_exp_phase1_raw = std::move(q_rlwes_raw);
#endif
}

void onionpir_query_expand_phase2(const MODULE* module,                //
                                  onionpir_expanded_query& query_exp,  //
                                  const onionpir_cloud_key& ckey,      //
                                  const onionpir_input_query& input_query) {
  const uint64_t nn = module_get_n(module);

#ifdef TEST_MODE
  for (auto& v : query_exp.query_exp_phase2_raw) {
    v.clear();
  }
#endif

  // extract each query column bit as RLWE
  const uint64_t query2bits = uint64_t(ceil(log2(ONIONPIR_db_ncols)));

  std::vector<int64_t> q_rlwes_raw(query2bits * ONIONPIR_query2exp_ellt * ONIONPIR_query2_ncols * nn);
  int64_t* q_rlwes = q_rlwes_raw.data();
  onionpir_trace_expand(module, q_rlwes, query2bits * ONIONPIR_query2exp_ellt, ONIONPIR_query2_ncols,
                        input_query.phase2.data(), ONIONPIR_query2_ncols, ckey);

  uint8_t* tmp_space =
      get_tmp_space(vmp_prepare_contiguous_tmp_bytes(module, ONIONPIR_query2exp_nrows, ONIONPIR_query2exp_ncols));
  for (uint64_t i = 0; i < query2bits; ++i) {
    std::vector<int64_t> q_rgsw_raw(ONIONPIR_query2exp_nrows * ONIONPIR_query2exp_ncols * nn);
    int64_t* q_rgsw = q_rgsw_raw.data();

    onionpir_prgsw_to_rgsw(module, ckey, q_rgsw, ONIONPIR_query2exp_nrows, ONIONPIR_query2exp_ncols,
                           q_rlwes + i * ONIONPIR_query2exp_ellt * ONIONPIR_query2_ncols * nn, ONIONPIR_query2exp_ellt,
                           ONIONPIR_query2_ncols);

    if (query_exp.query_exp_phase2[i] != nullptr) free(query_exp.query_exp_phase2[i]);
    query_exp.query_exp_phase2[i] = vmp_pmat_alloc(module, ONIONPIR_query2exp_nrows, ONIONPIR_query2exp_ncols);
    vmp_prepare_contiguous(module, query_exp.query_exp_phase2[i], q_rgsw, ONIONPIR_query2exp_nrows,
                           ONIONPIR_query2exp_ncols, tmp_space);

#ifdef TEST_MODE
    query_exp.query_exp_phase2_raw[i] = std::move(q_rgsw_raw);
#endif
  }
}
