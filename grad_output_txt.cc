#include "grad_output_txt.h"
#include "field.h"
#include <stdio.h>

Grad_Output_Txt::Grad_Output_Txt(Field* f0, Field* p0, Field* g0, const char* outputs) {
	f = f0;
	p = p0;
	g = g0;
	Grad_Output::parse_outputs(outputs);
}

void Grad_Output_Txt::write_output(const char* filename){
	// write the output to a text file
}