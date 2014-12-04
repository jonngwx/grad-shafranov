#include "field.h"
#include "grid.h"
#include "grad_output_txt.h"
#include <stdlib.h>

int main(){
	Field *f = new Field( 10, 11);
    Grid *g = new Grid(0, 0, 0 ,0 ,10, 11);
	Grad_Output_Txt out(f, g,f, "foo,bar");
	
	out.write_output("testing");
	delete f;
}