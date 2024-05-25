#include <spqlios/reim/reim_fft.h>

#include <iostream>
#include <vector>

int main() {
  // I solemnly swear that my entries are smaller than 2^50
  uint64_t Log2Bound = 50;

  uint64_t N = 16;     // N is the real dimension
  uint64_t m = N / 2;  // m is the complex dimension

  // let's multiply those two integer polynomials
  std::vector<int64_t> a = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  std::vector<int64_t> b = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // compute a_fft
  std::vector<double> a_fft(N);
  // TODO: Convert a from ZnX to reim then apply fft. Store in a_fft

  // compute b_fft
  std::vector<double> b_fft(N);
  // TODO: Convert b from ZnX to reim then apply fft. Store in b_fft

  // multiply in FFT space
  std::vector<double> product_fft(N);
  // TODO: Compute the multiplication in FFT space. Can use the "simple" variant

  // inverse FFT (not scaled by 1/m)
  // TODO: Apply ifft

  // scale by 1/m and round to int
  std::vector<int64_t> product(N);
  // TODO: Convert from reim format back to ZnX

  for (uint64_t i = 0; i < N; ++i) {
    std::cout << product[i] << ",";
  }
}
