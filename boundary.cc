#include "include/boundary.h"
Boundary::Boundary(Grid* grid)
  : nr_(grid->nr_),
    nz_(grid->nz_) {}

Boundary::~Boundary()
{}

int Boundary::CalcB(Field* psi, Field* jphi) 
{}

int Boundary::LtoI(int l) {
  int i=0;
   
  if (l>=0 && l<=(nr_-2)) {
    i = l;
  }
  else if (l>=(nr_-1) && l<=(nr_+nz_-3)) {
    i = nr_ -1;
  }
  else if (l>=(nr_+nz_-2) && l<=(2*nr_+nz_-4)) {
    i = (2*nr_+nz_-3) - l;
  }
  else {
    i = 0;
  }
  return i;
}

int Boundary::LtoJ(int l) {
  int j=0;
 
  if (l>=0 && l<=(nr_-2)) {
    j = 0;
  }
  else if (l>=(nr_-1) && l<=(nr_+nz_-3)) {
    j = l - (nr_-1);
  }
  else if (l>=(nr_+nz_-2) && l<=(2*nr_+nz_-4)) {
    j = nz_-1;
  }
  else {
    j = (2*nr_ + 2*nz_ - 4) - l;
  }
  return j;
}
  
