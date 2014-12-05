#ifndef SLOWBOUNDARY_H_
#define SLOWBOUNDARY_H_

#include "field.h"
#include "grid.h"
#include "tsv_reader.h"
#include "boundary.h"
/* need to include header that defines TsvData */

class SlowBoundary : public Boundary {
  public:
    SlowBoundary(Grid &grid, TsvData &cond_data); 
    ~SlowBoundary();
    int CalcB(Field &psi, Field &jphi);
  private:
    /* to be determined */

};

#endif // SLOWBOUNDARY_H_

