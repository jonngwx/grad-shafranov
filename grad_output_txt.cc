#include "grad_output_txt.h"
#include "field.h"
#include <stdio.h>

Grad_Output_Txt::Grad_Output_Txt(Field* f0, Field* p0, Field* g0, const char* outputs) {
	f = f0;
	p = p0;
	g = g0;
	Grad_Output::parse_outputs(outputs);
}

Grad_Output_Txt::~Grad_Output_Txt(){
}

void Grad_Output_Txt::write_output(const char* filename){
	FILE * file;
	file = fopen(filename, "w");
	const int nr = f->nr;
	// write the output to a text file

	for (int i = 0; i < nr; i++){
		fprintf(file,"%15.8f ", f->R[i]);
	}
	
	fclose(file);
}

