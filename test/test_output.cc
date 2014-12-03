#include "field.h"
#include "grad_output_txt.h"
#include <stdlib.h>

int main(){
	Field *f = new Field(.5, 1., -1., 1., 10, 10);
	Grad_Output_Txt g(f, f,f, "foo");
	
	g.write_output("testing");
	delete f;
}