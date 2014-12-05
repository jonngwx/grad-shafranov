#include "field.h"
#include "grid.h"
#include "grad_output_txt.h"
#include <stdlib.h>
#include "rhs_func.h"
#include <stdio.h>

int main(){
	Field *f = new Field( 10, 11);
    PGData *pdata = new PGData();
    double bla[2] = {1.,2.};
    pdata->psis = bla;
    pdata->values = bla;

    RHSfunc *p = new RHSfunc("blah", pdata);
    Grid *g = new Grid(0, 0, 0 ,0 ,10, 11);
	Grad_Output_Txt out(f, g,p,p, "foo,bar");
	
	out.write_output("testing");
	delete f;
    delete pdata;
    delete p;
    delete g;
}