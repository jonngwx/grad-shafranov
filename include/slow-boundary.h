#ifndef SLOWBOUNDARY_H_
#define SLOWBOUNDARY_H_

#include "boundary.h"
/* need to include header that defines TsvData */

class SlowBoundary : public Boundary {
  public:
    SlowBoundary(int nr, int nz, double dr, double dz, TsvData &cond_data); 
    ~SlowBoundary();
    int CalcB(Field &Psi);
  private:
    /* to be determined */

};

#endif // SLOWBOUNDARY_H_

