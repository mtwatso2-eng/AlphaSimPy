#include <armadillo>
#include <bitset>
#include <cmath>

std::bitset<8> toBits(unsigned char byte) {
  return std::bitset<8>(byte);
}

unsigned char toByte(std::bitset<8> bits) {
  return bits.to_ulong();
}

arma::uvec sampleInt(arma::uword n, arma::uword N) {
  arma::uvec output;
  output.set_size(n);
  if (n == 0) {
    return output;
  }
  double q, v, x, y1, y2;
  arma::uword threshold = 13 * n;
  arma::uword S, limit, top, bottom;
  arma::vec u(1, arma::fill::randu);
  v = std::exp(std::log(u(0)) / double(n));
  q = double(N - n + 1);
  while ((n > 1) & (threshold < N)) {
    while (true) {
      while (true) {
        x = double(N) * (1 - v);
        S = std::floor(x);
        if (double(S) < q) {
          break;
        }
        u.randu();
        v = std::exp(std::log(u(0)) / double(n));
      }
      u.randu();
      y1 = std::exp(std::log(u(0) * double(N) / q) / double(n - 1));
      v = y1 * (1 - x / double(N)) * (q / (q - double(S)));
      if (v <= 1) {
        break;
      }
      y2 = 1;
      top = N - 1;
      if ((n - 1) > S) {
        bottom = N - n;
        limit = N - S;
      } else {
        bottom = N - S - 1;
        limit = N - n + 1;
      }
      for (arma::uword i = N - 1; i >= limit; --i) {
        y2 *= double(top) / double(bottom);
      }
      u.randu();
      if ((double(N) / (double(N) - x)) >=
          (y1 * std::exp(std::log(y2) / double(n - 1)))) {
        v = std::exp(std::log(u(0)) / double(n - 1));
        break;
      }
      v = std::exp(std::log(u(0)) / double(n));
    }
    output(n - 1) = S + 1;
    N = N - S - 1;
    --n;
    q = double(N - n + 1);
    threshold -= 13;
  }
  if (n > 1) {
    top = N - n;
    while (n >= 2) {
      u.randu();
      S = 0;
      q = double(top) / double(N);
      while (q > u(0)) {
        ++S;
        --top;
        --N;
        q = (q * double(top)) / double(N);
      }
      output(n - 1) = S + 1;
      --N;
      --n;
    }
    u.randu();
    output(0) = std::floor(u(0) * N);
  } else {
    output(0) = std::floor(v * N);
  }
  return arma::cumsum(output);
}
