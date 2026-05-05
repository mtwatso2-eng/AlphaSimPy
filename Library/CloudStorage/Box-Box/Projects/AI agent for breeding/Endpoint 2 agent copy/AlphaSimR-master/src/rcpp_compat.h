#ifndef RCPP_COMPAT_H
#define RCPP_COMPAT_H

#include <stdexcept>
#include <iostream>

namespace Rcpp {
static std::ostream& Rcerr = std::cerr;

[[noreturn]] inline void stop(const char* msg) {
  throw std::runtime_error(msg);
}

[[noreturn]] inline void stop(const std::string& msg) {
  throw std::runtime_error(msg);
}
}  // namespace Rcpp

#endif
