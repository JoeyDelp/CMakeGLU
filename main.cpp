#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include <nicslu.h>
#include "numeric.h"
#include "preprocess.h"
#include "symbolic.h"

using namespace std;

void help_message() {
  cout << endl;
  cout << "GLU program V3.0" << endl;
  cout << "Usage: ./lu_cmd -i inputfile" << endl;
  cout << "Additional usage: ./lu_cmd -i inputfile -p" << endl;
  cout << "-p to enable perturbation" << endl;
}

int main(int argc, char **argv) {
  double utime(0.0);
  SNicsLU *nicslu;

  char *matrixName = NULL;
  bool PERTURB = false;

  double *ax = NULL;
  unsigned int *ai = NULL, *ap = NULL;
  unsigned int n;

  if (argc < 3) {
    help_message();
    return -1;
  }

  for (int i = 1; i < argc;) {
    if (strcmp(argv[i], "-i") == 0) {
      if (i + 1 > argc) {
        help_message();
        return -1;
      }
      matrixName = argv[i + 1];
      i += 2;
    } else if (strcmp(argv[i], "-p") == 0) {
      PERTURB = true;
      i += 1;
    } else {
      help_message();
      return -1;
    }
  }

  nicslu = (SNicsLU *)malloc(sizeof(SNicsLU));

  int err = preprocess(matrixName, nicslu, &ax, &ai, &ap);
  if (err) {
    // cout << "Reading matrix error" << endl;
    exit(1);
  }

  n = nicslu->n;

  cout << "Matrix Row: " << n << endl;
  cout << "Original nonzero: " << nicslu->nnz << endl;

  Symbolic_Matrix A_sym(n, cout, cerr);
  A_sym.fill_in(ai, ap);
  cout << "Symbolic time: " << utime << " ms" << endl;

  A_sym.csr();
  cout << "CSR time: " << utime << " ms" << endl;

  A_sym.predictLU(ai, ap, ax);
  cout << "PredictLU time: " << utime << " ms" << endl;

  A_sym.leveling();
  cout << "Leveling time: " << utime << " ms" << endl;

#if GLU_DEBUG
  A_sym.ABFTCalculateCCA();
//    A_sym.PrintLevel();
#endif

  LUonDevice(A_sym, cout, cerr, PERTURB);

#if GLU_DEBUG
  A_sym.ABFTCheckResult();
#endif

  // solve Ax=b
  vector<REAL> b(n, 1.);
  vector<REAL> x = A_sym.solve(nicslu, b);
  {
    ofstream x_f("x.dat");
    for (double xx : x) x_f << xx << '\n';
  }
}
