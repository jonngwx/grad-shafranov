#ifndef GRAD_OUTPUT_TXT
#define GRAD_OUTPUT_TXT

#include <stdlib.h>
#include "grad_output.h"
#include "field.h"

class Grad_Output_Txt : public Grad_Output{
	public:
	Grad_Output_Txt(Field* f, Field* p, Field* g, const char* outputs);
	~Grad_Output_Txt();
	
	write_output(const char* filename);
	
};

#endif