/**
 * AlphaSimPy C++ bindings - PyBind11 wrapper for MaCS simulation
 * Calls AlphaSimR's runFromAlphaSimR and implements the MaCS loop
 * to convert output to Python/numpy types.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <bitset>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

// AlphaSimR headers - standalone Armadillo
#include <armadillo>
#include "simulator.h"
#include "misc.h"

// runFromAlphaSimR is in simulator.cpp but not declared in simulator.h
std::vector<AlphaSimRReturn> runFromAlphaSimR(std::string in);

py::dict macs(const std::string& command,
              py::array_t<uint32_t> max_sites,
              bool inbred,
              unsigned int ploidy,
              int n_threads,
              const std::vector<std::string>& seed) {
  py::buffer_info buf = max_sites.request();
  if (buf.ndim != 1) {
    throw std::runtime_error("max_sites must be a 1D array");
  }
  size_t n_chr = buf.shape[0];
  const uint32_t* max_sites_ptr = static_cast<const uint32_t*>(buf.ptr);

  arma::uvec maxSites(n_chr);
  for (size_t i = 0; i < n_chr; i++) {
    maxSites(i) = max_sites_ptr[i];
  }

  if (n_threads <= 0) {
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#else
    n_threads = 1;
#endif
  }

  arma::field<arma::Cube<unsigned char>> geno(n_chr);
  arma::field<arma::vec> genMap(n_chr);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (arma::uword chr = 0; chr < n_chr; chr++) {
    std::vector<AlphaSimRReturn> macsOutput = runFromAlphaSimR(command + seed[chr]);

    arma::uword nSites = macsOutput.size();
    arma::uword nHap = macsOutput[0].haplotypes.size();
    arma::uword nInd = inbred ? nHap : nHap / ploidy;

    arma::uvec selSites;
    if (maxSites(chr) > 0) {
      if (nSites < maxSites(chr)) {
        maxSites(chr) = nSites;
      }
      selSites = sampleInt(maxSites(chr), nSites);
      nSites = maxSites(chr);
    } else {
      selSites.set_size(nSites);
      for (arma::uword i = 0; i < nSites; ++i) {
        selSites(i) = i;
      }
    }

    genMap(chr).set_size(nSites);
    for (arma::uword site = 0; site < nSites; ++site) {
      genMap(chr).at(site) = macsOutput[selSites(site)].length;
    }

    arma::uword nBins = nSites / 8;
    if ((nSites % 8) > 0) ++nBins;

    geno(chr).set_size(nBins, ploidy, nInd);
    arma::uword grp = 0;
    arma::uword locus;
    std::bitset<8> workBits;

    for (arma::uword i = 0; i < nHap; ++i) {
      if (inbred) {
        for (grp = 0; grp < ploidy; ++grp) {
          locus = 0;
          if (nBins > 1) {
            for (arma::uword j = 0; j < (nBins - 1); ++j) {
              for (arma::uword k = 0; k < 8; ++k) {
                workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
                ++locus;
              }
              geno(chr).slice(i).col(grp).row(j) = toByte(workBits);
            }
          }
          for (arma::uword k = 0; k < 8; ++k) {
            workBits[k] = (locus < nSites) ? macsOutput[selSites(locus)].haplotypes[i] : 0;
            if (locus < nSites) ++locus;
          }
          geno(chr).slice(i).col(grp).row(nBins - 1) = toByte(workBits);
        }
      } else {
        locus = 0;
        if (nBins > 1) {
          for (arma::uword j = 0; j < (nBins - 1); ++j) {
            for (arma::uword k = 0; k < 8; ++k) {
              workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
              ++locus;
            }
            geno(chr).slice(i / ploidy).col(grp).row(j) = toByte(workBits);
          }
        }
        for (arma::uword k = 0; k < 8; ++k) {
          workBits[k] = (locus < nSites) ? macsOutput[selSites(locus)].haplotypes[i] : 0;
          if (locus < nSites) ++locus;
        }
        geno(chr).slice(i / ploidy).col(grp).row(nBins - 1) = toByte(workBits);
        ++grp;
        grp = grp % ploidy;
      }
    }
  }

  py::list geno_list;
  py::list gen_map_list;

  for (arma::uword chr = 0; chr < n_chr; chr++) {
    arma::Cube<unsigned char>& g = geno(chr);
    py::array_t<unsigned char> geno_arr({g.n_rows, g.n_cols, g.n_slices});
    py::buffer_info geno_buf = geno_arr.request();
    unsigned char* geno_ptr = static_cast<unsigned char*>(geno_buf.ptr);
    for (size_t i = 0; i < g.n_elem; i++) {
      geno_ptr[i] = g(i);
    }
    geno_list.append(geno_arr);

    arma::vec& gm = genMap(chr);
    py::array_t<double> gen_map_arr(gm.n_elem);
    py::buffer_info gm_buf = gen_map_arr.request();
    double* gm_ptr = static_cast<double*>(gm_buf.ptr);
    for (size_t i = 0; i < gm.n_elem; i++) {
      gm_ptr[i] = gm(i);
    }
    gen_map_list.append(gen_map_arr);
  }

  py::dict result;
  result["geno"] = geno_list;
  result["gen_map"] = gen_map_list;
  return result;
}

int get_num_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

py::array_t<uint32_t> sample_int(size_t n, size_t N) {
  if (n == 0) {
    std::vector<size_t> shape = {0};
    return py::array_t<uint32_t>(shape);
  }
  if (n >= N) {
    py::array_t<uint32_t> result(N);
    auto r = result.mutable_unchecked<1>();
    for (size_t i = 0; i < N; i++) r(i) = static_cast<uint32_t>(i);
    return result;
  }
  arma::uvec samples = sampleInt(static_cast<arma::uword>(n), static_cast<arma::uword>(N));
  py::array_t<uint32_t> result(n);
  auto r = result.mutable_unchecked<1>();
  for (size_t i = 0; i < n; i++) {
    r(i) = static_cast<uint32_t>(samples(i));
  }
  return result;
}

PYBIND11_MODULE(_alphasimpy_cpp, m) {
  m.doc() = "AlphaSimPy C++ bindings for MaCS simulation";
  m.def("macs", &macs,
        py::arg("command"), py::arg("max_sites"), py::arg("inbred"),
        py::arg("ploidy"), py::arg("n_threads"), py::arg("seed"));
  m.def("get_num_threads", &get_num_threads);
  m.def("sample_int", &sample_int, py::arg("n"), py::arg("N"));
}
