#include <spqlios/reim/reim_fft.h>
#include <vector>
#include <iostream>

int main() {
    // I solemnly swear that my entries are smaller than 2^50
    uint64_t Log2Bound = 50;

    uint64_t N = 16;   // N is the real dimension
    uint64_t m = N/2;  // m is the complex dimension

    // let's multiply those two integer polynomials
    std::vector<int64_t> a = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    std::vector<int64_t> b = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    // compute a_fft
    std::vector<double> a_fft(N);
    reim_from_znx64_simple(m, Log2Bound, a_fft.data(), a.data());
    reim_fft_simple(m, a_fft.data());

    // compute b_fft
    std::vector<double> b_fft(N);
    reim_from_znx64_simple(m, Log2Bound, b_fft.data(), b.data());
    reim_fft_simple(m, b_fft.data());

    // multiply in FFT space
    std::vector<double> product_fft(N);
    reim_fftvec_mul_simple(m, product_fft.data(), a_fft.data(), b_fft.data());

    // inverse FFT (not scaled by 1/m)
    reim_ifft_simple(m, product_fft.data());

    // scale by 1/m and round to int
    std::vector<int64_t> product(N);
    reim_to_znx64_simple(m, m, Log2Bound, product.data(), product_fft.data());

    for (uint64_t i=0; i<N; ++i) {
        std::cout << product[i] << ",";
    }
}
