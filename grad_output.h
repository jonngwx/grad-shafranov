#ifndef GRAD_OUTPUT
#define GRAD_OUTPUT

#include <stdlib.h>

class Grad_Output{
	public:
	virtual ~Grad_Output(){};
	virtual write_output(const char* filename)=0;
};

#endif