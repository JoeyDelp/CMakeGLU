#ifndef _SYMBOLIC_H_
#define _SYMBOLIC_H_
#include <iostream>
#include <unordered_set>
#include <vector>

#include "nicslu.h"
#include "type.h"

class Symbolic_Matrix {
 public:
  unsigned int n;
  unsigned int nnz;
  int num_lev;
  std::vector<unsigned> sym_c_ptr;
  std::vector<unsigned> sym_r_idx;
  std::vector<unsigned> csr_r_ptr;
  std::vector<unsigned> csr_c_idx;
  std::vector<unsigned> csr_diag_ptr;
  std::vector<REAL> val;
  std::vector<unsigned> l_col_ptr;  // Indices of diagonal elements
  std::vector<int> level_idx;
  std::vector<int> level_ptr;

  void predictLU(unsigned *, unsigned *, double *);
  void csr();
  void leveling();
  void fill_in(unsigned *, unsigned *);
  std::vector<REAL> solve(SNicsLU *, const std::vector<real__t> &);

#if GLU_DEBUG
  void PrintLevel();
  // ABFT
  std::vector<REAL> CCA;
  void ABFTCalculateCCA();
  void ABFTCheckResult();
#endif

  Symbolic_Matrix &operator=(const Symbolic_Matrix &other) {
    if (this != &other) {  // Self-assignment check
      n = other.n;
      nnz = other.nnz;
      num_lev = other.num_lev;
      sym_c_ptr = other.sym_c_ptr;
      sym_r_idx = other.sym_r_idx;
      csr_r_ptr = other.csr_r_ptr;
      csr_c_idx = other.csr_c_idx;
      csr_diag_ptr = other.csr_diag_ptr;
      val = other.val;
      l_col_ptr = other.l_col_ptr;
      level_idx = other.level_idx;
      level_ptr = other.level_ptr;
    }
    return *this;
  }

  Symbolic_Matrix(unsigned n, std::ostream &out, std::ostream &err)
      : n(n), num_lev(0), m_out(out), m_err(err){};

  Symbolic_Matrix(std::ostream &out, std::ostream &err)
      : m_out(out), m_err(err){};

  Symbolic_Matrix(){};

 private:
  std::ostream &m_out = std::cout;
  std::ostream &m_err = std::cerr;
};

#endif
