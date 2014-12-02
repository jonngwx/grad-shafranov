#include "grad_output_txt.h"
#include "field.h"
#include <stdio.h>

Grad_Output_Txt::Grad_Output_Txt(Field* f, Field* p, Field* g, const char* outputs) : f(f), p(p), g(g){
	Grad_Output::parse_outputs(outputs);
}

void Grad_Output_Txt::write_output(const char* filename){
	// write the output to a text file
}