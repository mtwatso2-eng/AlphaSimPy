#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <armadillo>  // Use Armadillo directly instead of RcppArmadillo
#include <vector>
#include <string>
#include <bitset>
#include <stdexcept>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "simulator.h"
#include "misc.h"

// Forward declarations
std::vector<AlphaSimRReturn> runFromAlphaSimR(const std::string& in);
arma::uvec sampleInt(arma::uword n, arma::uword N);
int getNumThreads();

// Forward declaration for createDH2 - we'll need to implement a wrapper
// since the original uses Rcpp types
// For now, we'll implement makeDH in Python using a simplified approach

// Wrapper for runFromAlphaSimR that replaces Rcpp::stop with std::runtime_error
std::vector<AlphaSimRReturn> runFromAlphaSimR_python(const std::string& in) {
    std::vector<std::string> words;
    Simulator simulator;
    
    if (in == "" || in.empty()) {
        throw std::runtime_error("Not enough args for macs call");
    }
    
    boost::split(words, in, boost::is_any_of(", "), boost::token_compress_on);
    CommandArguments arguments;
    std::vector<std::string> subOption;
    
    // sample size
    subOption.emplace_back(words[0]);
    // seq length
    subOption.emplace_back(words[1]);
    arguments.push_back(subOption);
    subOption.clear();
    
    for (unsigned int i = 2; i < words.size(); ++i) {
        subOption.emplace_back(words[i]);
        if (i == words.size() - 1 || (words[i+1][0] == '-' && words[i+1][1] >= 65)) {
            arguments.push_back(subOption);
            subOption.clear();
        }
    }
    
    if (arguments.size() == 0) {
        throw std::runtime_error("Not enough args for macs call");
    }
    
    simulator.readInputParameters(arguments);
    std::vector<AlphaSimRReturn> test = simulator.beginSimulationMemory();
    return test;
}

namespace py = pybind11;

// Convert NumPy array to Armadillo uvec
arma::uvec numpy_to_uvec(py::array_t<unsigned int> arr) {
    auto buf = arr.request();
    arma::uvec result(static_cast<arma::uword>(buf.size));
    unsigned int* ptr = static_cast<unsigned int*>(buf.ptr);
    for (size_t i = 0; i < buf.size; i++) {
        result(i) = ptr[i];
    }
    return result;
}

// Convert Armadillo uvec to NumPy array
py::array_t<unsigned int> uvec_to_numpy(const arma::uvec& vec) {
    auto result = py::array_t<unsigned int>(vec.n_elem);
    auto buf = result.request();
    unsigned int* ptr = static_cast<unsigned int*>(buf.ptr);
    for (size_t i = 0; i < vec.n_elem; i++) {
        ptr[i] = vec(i);
    }
    return result;
}

// Convert NumPy array to Armadillo vec
arma::vec numpy_to_vec(py::array_t<double> arr) {
    auto buf = arr.request();
    arma::vec result(static_cast<arma::uword>(buf.size));
    double* ptr = static_cast<double*>(buf.ptr);
    for (size_t i = 0; i < buf.size; i++) {
        result(i) = ptr[i];
    }
    return result;
}

// Convert Armadillo vec to NumPy array
py::array_t<double> vec_to_numpy(const arma::vec& vec) {
    auto result = py::array_t<double>(vec.n_elem);
    auto buf = result.request();
    double* ptr = static_cast<double*>(buf.ptr);
    for (size_t i = 0; i < vec.n_elem; i++) {
        ptr[i] = vec(i);
    }
    return result;
}

// Convert NumPy 3D array to Armadillo Cube
arma::Cube<unsigned char> numpy_to_cube(py::array_t<unsigned char> arr) {
    auto buf = arr.request();
    if (buf.ndim != 3) {
        throw std::runtime_error("Input array must be 3-dimensional");
    }
    
    arma::uword n_bins = buf.shape[0];
    arma::uword ploidy = buf.shape[1];
    arma::uword n_ind = buf.shape[2];
    
    arma::Cube<unsigned char> result(n_bins, ploidy, n_ind);
    unsigned char* ptr = static_cast<unsigned char*>(buf.ptr);
    
    for (arma::uword i = 0; i < n_bins; i++) {
        for (arma::uword j = 0; j < ploidy; j++) {
            for (arma::uword k = 0; k < n_ind; k++) {
                size_t idx = i * ploidy * n_ind + j * n_ind + k;
                result(i, j, k) = ptr[idx];
            }
        }
    }
    return result;
}

// Convert Armadillo Cube to NumPy 3D array
py::array_t<unsigned char> cube_to_numpy(const arma::Cube<unsigned char>& cube) {
    std::vector<size_t> shape = {cube.n_rows, cube.n_cols, cube.n_slices};
    auto result = py::array_t<unsigned char>(shape);
    auto buf = result.request();
    unsigned char* ptr = static_cast<unsigned char*>(buf.ptr);
    
    for (arma::uword i = 0; i < cube.n_rows; i++) {
        for (arma::uword j = 0; j < cube.n_cols; j++) {
            for (arma::uword k = 0; k < cube.n_slices; k++) {
                size_t idx = i * cube.n_cols * cube.n_slices + j * cube.n_slices + k;
                ptr[idx] = cube(i, j, k);
            }
        }
    }
    return result;
}

// toByte is already defined in misc.cpp

// Python wrapper for MaCS function
py::dict macs_wrapper(const std::string& args,
                      py::array_t<unsigned int> max_sites,
                      bool inbred,
                      unsigned int ploidy,
                      int n_threads,
                      const std::vector<std::string>& seed) {
    // Convert inputs
    arma::uvec maxSites = numpy_to_uvec(max_sites);
    arma::uword ploidy_arma = static_cast<arma::uword>(ploidy);
    
    arma::uword nChr = maxSites.n_elem;
    arma::field<arma::Cube<unsigned char>> geno(nChr);
    arma::field<arma::vec> genMap(nChr);
    
    // Loop through chromosomes
    for (arma::uword chr = 0; chr < nChr; chr++) {
        std::string seed_str = (chr < seed.size()) ? seed[chr] : seed[0];
        std::string full_args = args + seed_str;
        
        // Run MaCS simulation
        std::vector<AlphaSimRReturn> macsOutput = runFromAlphaSimR_python(full_args);
        
        arma::uword nSites, nBins, nHap, nInd;
        nSites = macsOutput.size();
        if (nSites == 0) {
            throw std::runtime_error("MaCS simulation returned no sites");
        }
        nHap = macsOutput[0].haplotypes.size();
        if (inbred) {
            nInd = nHap;
        } else {
            nInd = nHap / ploidy_arma;
        }
        
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
        
        // Fill genMap
        genMap(chr).set_size(nSites);
        for (arma::uword site = 0; site < nSites; ++site) {
            genMap(chr).at(site) = macsOutput[selSites(site)].length;
        }
        
        // Fill Geno
        nBins = nSites / 8;
        if ((nSites % 8) > 0) {
            ++nBins;
        }
        
        geno(chr).set_size(nBins, ploidy_arma, nInd);
        arma::uword grp = 0;
        arma::uword locus;
        std::bitset<8> workBits;
        
        for (arma::uword i = 0; i < nHap; ++i) {
            if (inbred) {
                for (grp = 0; grp < ploidy_arma; ++grp) {
                    locus = 0;
                    // Fill in bins known to be complete
                    if (nBins > 1) {
                        for (arma::uword j = 0; j < (nBins - 1); ++j) {
                            for (arma::uword k = 0; k < 8; ++k) {
                                workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
                                ++locus;
                            }
                            geno(chr).slice(i).col(grp).row(j) = toByte(workBits);
                        }
                    }
                    // Fill in potentially incomplete bins
                    for (arma::uword k = 0; k < 8; ++k) {
                        if (locus < nSites) {
                            workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
                            ++locus;
                        } else {
                            workBits[k] = 0;
                        }
                    }
                    geno(chr).slice(i).col(grp).row(nBins - 1) = toByte(workBits);
                }
            } else {
                locus = 0;
                // Fill in bins known to be complete
                if (nBins > 1) {
                    for (arma::uword j = 0; j < (nBins - 1); ++j) {
                        for (arma::uword k = 0; k < 8; ++k) {
                            workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
                            ++locus;
                        }
                        geno(chr).slice(i / ploidy_arma).col(grp).row(j) = toByte(workBits);
                    }
                }
                // Fill in potentially incomplete bins
                for (arma::uword k = 0; k < 8; ++k) {
                    if (locus < nSites) {
                        workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
                        ++locus;
                    } else {
                        workBits[k] = 0;
                    }
                }
                geno(chr).slice(i / ploidy_arma).col(grp).row(nBins - 1) = toByte(workBits);
                ++grp;
                grp = grp % ploidy_arma;
            }
        }
    }
    
    // Convert results to Python dict
    py::dict result;
    py::list geno_list;
    py::list gen_map_list;
    
    for (arma::uword i = 0; i < nChr; i++) {
        geno_list.append(cube_to_numpy(geno(i)));
        gen_map_list.append(vec_to_numpy(genMap(i)));
    }
    
    result["geno"] = geno_list;
    result["gen_map"] = gen_map_list;
    
    return result;
}

PYBIND11_MODULE(_alphasimpy_cpp, m) {
    m.doc() = "AlphaSimPy C++ bindings via PyBind11";
    
    m.def("macs", &macs_wrapper,
          "Run MaCS simulation",
          py::arg("args"),
          py::arg("max_sites"),
          py::arg("inbred"),
          py::arg("ploidy"),
          py::arg("n_threads"),
          py::arg("seed"));
    
    m.def("get_num_threads", &getNumThreads,
          "Get number of available threads");
    
    m.def("sample_int", [](unsigned int n, unsigned int N) {
        arma::uvec result = sampleInt(n, N);
        return uvec_to_numpy(result);
    }, "Sample n integers without replacement from 0 to N-1",
    py::arg("n"), py::arg("N"));
}
