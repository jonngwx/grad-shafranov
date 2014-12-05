#include "grad_output_txt.h"
#include "field.h"
#include "grid.h"
#include <stdio.h>

Grad_Output_Txt::Grad_Output_Txt(Field* f0, Grid* grid0, RHSfunc * p0, RHSfunc * g0, const char* outputs) {
    f = f0;
    p = p0;
    g = g0;
    grid = grid0;
    Grad_Output::parse_outputs(outputs);
}

Grad_Output_Txt::~Grad_Output_Txt(){
}

void Grad_Output_Txt::write_output(const char* filename){
	FILE * file;
	file = fopen(filename, "w");
	const int nr = grid->nr;
    const int nz = grid->nz;
	// write the output to a text file
	for (int i = 0; i < nr; i++){
		fprintf(file,"%15.8f ", grid->R[i]);
	}
    fprintf(file,"\n");
    for (int i = 0; i < nz; i++){
		fprintf(file,"%15.8f ", grid->z[i]);
	}
    fprintf(file,"\n");
    for (int i = 0; i < nr; i++){
        for (int j = 0; j < nz; j++){
            fprintf(file,"%15.8f ", f->f[i][j]);
        }
        fprintf(file,"\n");
    }
    
    for (int i = 0; i < nr; i++){
        for (int j = 0; j < nz; j++){
            fprintf(file,"%15.8f ", g->eval(f->f[i][j]));
        }
        fprintf(file,"\n");
    }
    
	// do all other fancy outputs
    // enums?
	fclose(file);
}

