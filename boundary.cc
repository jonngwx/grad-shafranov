/*!
 * @file boundary.cc
 * @author ???
 * @brief Implementation for the Boundary class.
 */
#include "include/boundary.h"
Boundary::Boundary(Grid* grid) : nr_(grid->nr_), nz_(grid->nz_) {}

Boundary::~Boundary() {}

//so this gets called in elliptic_test. The compiler gives a warning that there's no 'return' statement for this nonvoid function. Since this is not implemented, should I assume that our actual program will ever only call SlowBoundary's CalcB???
int Boundary::CalcB(Field* psi, Field* jphi) {}

/* L is numbered like:  This function will return:
 *
 * 6 5 4                0 1 2
 * 7   3                0   2
 * 0 1 2                0 1 2
 *
 * -helpful comment by JAS
 */
int Boundary::LtoI(int l) {
  int i = 0;

  if (l >= 0 && l <= (nr_ - 2)) {
    i = l;
  } else if (l >= (nr_ - 1) && l <= (nr_ + nz_ - 3)) {
    i = nr_ - 1;
  } else if (l >= (nr_ + nz_ - 2) && l <= (2 * nr_ + nz_ - 4)) {
    i = (2 * nr_ + nz_ - 3) - l;
  } else {
    i = 0;
  }
  return i;
}

/* L is numbered like:  This function will return:
 *
 * 6 5 4                2 2 2
 * 7   3                1   1
 * 0 1 2                0 0 0
 *
 * -helpful comment by JAS
 */
int Boundary::LtoJ(int l) {
  int j = 0;

  if (l >= 0 && l <= (nr_ - 2)) {
    j = 0;
  } else if (l >= (nr_ - 1) && l <= (nr_ + nz_ - 3)) {
    j = l - (nr_ - 1);
  } else if (l >= (nr_ + nz_ - 2) && l <= (2 * nr_ + nz_ - 4)) {
    j = nz_ - 1;
  } else {
    j = (2 * nr_ + 2 * nz_ - 4) - l;
  }
  return j;
}
