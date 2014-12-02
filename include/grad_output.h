#ifndef GRAD_OUTPUT
#define GRAD_OUTPUT

#include <stdlib.h>
#include "field.h"

class Grad_Output{
	public:
	virtual ~Grad_Output(){};
	virtual void write_output(const char* filename)=0;
	
protected:
	Field *f;
	Field *p;
	Field *g;
	void parse_outputs(const char *outputs);
	// output_list
};

#endif