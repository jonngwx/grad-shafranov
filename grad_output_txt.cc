#include "grad_output_txt.h"
#include "field.h"
#include "grid.h"
#include <stdio.h>
#include "util.h"

Grad_Output_Txt::Grad_Output_Txt(Field* f0, Field* jphi0, Grid* grid0, Field* p0, Field* g0, const char* outputs) {
    f = f0;
    jphi = jphi0;
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
	const int nr = grid->nr_;
        const int nz = grid->nz_;
	// write the output to a text file
        // format is name: data \n
    fprintf(file, "R: ");
    print1d(file,grid->R_,nr);
    fprintf(file,"\n");

    fprintf(file,"z: ");
    print1d(file,grid->z_,nz);
    fprintf(file,"\n");

    double psilo[2] = {f->f_l,f->f_0};
    fprintf(file,"psilo: ");
    print1d(file, psilo ,2);
    fprintf(file,"\n");

    fprintf(file,"psi: ");
    print2d(file,f->f_,nr,nz);
    fprintf(file,"\n");
    
    fprintf(file,"g: ");
    print2d(file,g->f_,nr,nz);
    fprintf(file,"\n");
    
    fprintf(file,"p: ");
    print2d(file,p->f_,nr,nz);
    fprintf(file,"\n");

	// do all other fancy outputs
    for (auto i : output_list){
        switch(i){
        case CURRENT:
            printf("current\n");
            fprintf(file,"j: ");
            for (int i = 0; i < nr; i++){
                for (int j = 0; j < nz; j++){
                    fprintf(file,"%15.8f ", jphi->f_[i][j]);
                }
            }
            fprintf(file,"\n");
            break;
        case TOROIDAL_FIELD:
            fprintf(file,"bt: ");
            for (int i = 0; i < nr; i++){
                for (int j = 0; j < nz; j++){
                    fprintf(file,"%15.8f ", g->f_[i][j]/grid->R_[i]);
                }
            }
            fprintf(file,"\n");
            break;
        }
    }
    fclose(file);
}

