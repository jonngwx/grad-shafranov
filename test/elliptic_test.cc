#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Elliptic
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/grid.h"
#include "../include/field.h"
#include "../include/elliptic_solver.h"
#include "../include/sor.h"
#include "../include/gauss_seidel.h"
#include <math.h>
#include "../include/boundary.h"

struct elliptictester {
  elliptictester() { R0= 0; }
  ~elliptictester() {}
  double R0;
};

BOOST_FIXTURE_TEST_SUITE(ell_suite, elliptictester)
/*!
 * Tests SOR elliptic solver by comparison to vacuum
 * solution constrained to be even about the Z = 0
 * midplane
 */
BOOST_AUTO_TEST_CASE (SOR_vacuum) {
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 10;
    double nz = 10;
    double max_iter = 100;
    double epsilon = 0.1;
    double perim = 2.0*(nr + nz -2.0);
    Grid *grid = new Grid(R0, Rend, Z0, Zend, nr, nz);
    Field *psi = new Field(*grid);
    Field *jphi = new Field(*grid);
    Boundary *psib = new Boundary(psi,grid);
    
    // Initialize psi and jphi
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            psi->f_[i][j] = 0;
            jphi->f_[i][j] = 0;
        }
    }
    
    // Vacuum solution multipole expansion
    Field *psi_sol = new Field(*grid);
    double Rc = 5;
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            double R = grid->R_[i];
            double Z = grid->z_[j];
            psi_sol->f_[i][j] = Rc*Rc + 0.5*(R*R - Rc*Rc);
            psi_sol->f_[i][j] += 1/(8*Rc*Rc)*((R*R-Rc*Rc)*(R*R-Rc*Rc) - 4*R*R*Z*Z);
            psi_sol->f_[i][j] += 1/(24*pow(Rc,4))*(pow((R*R-Rc*Rc),3) - 12*R*R*Z*Z*(R*R - Rc*Rc) + 8*R*R*pow(Z,4));
        }
    }
    
    // Initialize psi boundary using exact solution
    for (int l = 0; l < perim; ++l) {
        psi->f_[psib->LtoI(l)][psib->LtoJ(l)] = psi_sol->f_[psib->LtoI(l)][psib->LtoJ(l)];
    }
    
    /* 
    printf("\n SOLUTION \n");
    for (int i = 0; i < nr; ++i) {
        for (int j=0; j < nz; ++j) {
            printf("%f\t", psi_sol->f_[i][j]);
        }
        printf("\n");
    }
    */
    
    // Solve for psi given the proper boundary conditions
    EllipticSolver *solver = new SOR(*grid, *psi, epsilon);
    solver->coeff();
    for (int n = 0; n < max_iter; ++n) {
        if (n==0) solver->step_1(*jphi);
        else solver->step(*jphi);
        if (solver->norm() < epsilon) break;
    }
   
    //printf("\n FINAL \n");
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            //printf("%f\t", psi->f_[i][j]);
            BOOST_CHECK_CLOSE(psi_sol->f_[i][j], psi->f_[i][j], 4);
        }
        //printf("\n");
    }
}

/*!
 * Tests Gauss-Seidel elliptic solver by comparison to vacuum
 * solution constrained to be even about the Z = 0
 * midplane
 */

BOOST_AUTO_TEST_CASE (GS_vacuum) {
    double R0 = 0;
    double Rend = 10;
    double Z0 = -5;
    double Zend = 5;
    double nr = 10;
    double nz = 10;
    double max_iter = 100;
    double epsilon = 0.1;
    double error_ES = 0.1;
    double perim = 2.0*(nr + nz -2.0);
    Grid *grid = new Grid(R0, Rend, Z0, Zend, nr, nz);
    Field *psi = new Field(*grid);
    Field *jphi = new Field(*grid);
    Boundary *psib = new Boundary(psi, grid);
    
    // Initialize psi and jphi
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            psi->f_[i][j] = 0;
            jphi->f_[i][j] = 0;
        }
    }
    
    // Vacuum solution multipole expansion
    Field *psi_sol = new Field(*grid);
    double Rc = 5;
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            double R = grid->R_[i];
            double Z = grid->z_[j];
            psi_sol->f_[i][j] = Rc*Rc + 0.5*(R*R - Rc*Rc);
            psi_sol->f_[i][j] += 1/(8*Rc*Rc)*((R*R-Rc*Rc)*(R*R-Rc*Rc) - 4*R*R*Z*Z);
            psi_sol->f_[i][j] += 1/(24*pow(Rc,4))*(pow((R*R-Rc*Rc),3) - 12*R*R*Z*Z*(R*R - Rc*Rc) + 8*R*R*pow(Z,4));
        }
    }
    
    // Initialize psi boundary using exact solution
    for (int l = 0; l < perim; ++l) {
        psi->f_[psib->LtoI(l)][psib->LtoJ(l)] = psi_sol->f_[psib->LtoI(l)][psib->LtoJ(l)];
    }
    
    // Initialize psi boundary using exact solution
    for (int l = 0; l < perim; ++l) {
        psi->f_[psib->LtoI(l)][psib->LtoJ(l)] = psi_sol->f_[psib->LtoI(l)][psib->LtoJ(l)];
    }
   
    /* 
    printf("\n SOLUTION \n");
    for (int i = 0; i < nr; ++i) {
        for (int j=0; j < nz; ++j) {
            printf("%f\t", psi_sol->f_[i][j]);
        }
        printf("\n");
    }
    */
    
    // Solve for psi given the proper boundary conditions
    EllipticSolver *solver = new GaussSeidel(*grid, *psi, error_ES);
    solver->coeff();
    for (int n = 0; n < max_iter; ++n) {
        if (n==0) solver->step_1(*jphi);
        else solver->step(*jphi);
        if (solver->norm() < epsilon) break;
    }
    //printf("\n FINAL \n");
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            //printf("%f\t", psi->f_[i][j]);
            BOOST_CHECK_CLOSE(psi_sol->f_[i][j], psi->f_[i][j], 4);
        }
        //printf("\n");
    }
}


BOOST_AUTO_TEST_CASE (GS_shaf_solo) {
    /***********************************
     * Specify parameters for the grid
     **********************************/
    double Rmin = 5.0;
    double Rmax = 7.0;
    double zmin = -1.0;
    double zmax = 1.0;
    int nr = 20;
    int nz = 20;
  
    double max_iter = 1000;
    double epsilon = 1e-8;
    double error_ES = 0.01;
    double perim = 2.0*(nr + nz -2.0);
    Grid *grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    Field *psi = new Field(*grid);
    Field *psi_old = new Field(*grid);
    Field *jphi = new Field(*grid);
    Boundary *psib = new Boundary(psi, grid);
    /*********************************
     * Specify parameters for SS model
     *********************************/
    double B0 = 1.0; // B field in Teslas at mag axis
    double R0 = 6.0; // radial position of mag axis
    double kap0 = 1.1; // ellipticity of plasma
    double q0 = 2.01; // safety factor at mag axis
    double eps = 0.3; // inverse aspect ratio
    const double mu0 = 4*M_PI*1e-7; // mag perm of free space
    double psi_pv = eps*eps*kap0*R0*R0*B0/(2*q0); // defines psi at plas/vac interf.
    double R,Z;
    double coeff = B0/(R0*R0*kap0*q0);
    for (int i = 0; i < grid->nr_; ++i) {
      R = grid->R_[i];
      for (int j = 0; j < grid->nz_; ++j) {
	Z = grid->z_[j];
	psi_old->f_[i][j] = (coeff/2.0)*(R*R*Z*Z + (kap0*kap0/4.0)*pow((R*R-R0*R0),2));
	psi->f_[i][j] = psi_old->f_[i][j]*+1;
	if ((psi_old->f_[i][j] >= 0) && (psi_old->f_[i][j] <= psi_pv)) {
	  jphi->f_[i][j] = R*coeff*(kap0*kap0 + 1)/mu0;
	}
	else {
	  jphi->f_[i][j] = 0;
	} 
      }
    }
     
    // Initialize psi boundary using exact solution
    for (int l = 0; l < perim; ++l) {
        psi->f_[psib->LtoI(l)][psib->LtoJ(l)] = psi_old->f_[psib->LtoI(l)][psib->LtoJ(l)];
    }
       
    /* 
    printf("\n SOLUTION \n");
    for (int i = 0; i < nr; ++i) {
        for (int j=0; j < nz; ++j) {
            printf("%f\t", psi_sol->f_[i][j]);
        }
        printf("\n");
    }
    */
    
    // Solve for psi given the proper boundary conditions
    EllipticSolver *solver = new GaussSeidel(*grid, *psi, error_ES);
    solver->coeff();
    for (int n = 0; n < max_iter; ++n) {
        if (n==0) solver->step_1(*jphi);
        else solver->step(*jphi);
	//	printf("n = %d\n, norm = %f\n",n,solver->norm());
        if (solver->norm() < epsilon) {
	  break;

	}
    }
    //    printf("\n FINAL \n");
    for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nz; ++j) {
            //printf("%f\t", psi->f_[i][j]);
	       BOOST_CHECK_CLOSE(psi_old->f_[i][j], psi->f_[i][j], 1);
        }
        //printf("\n");
    }
}



BOOST_AUTO_TEST_SUITE_END()
