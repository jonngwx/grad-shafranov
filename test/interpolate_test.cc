#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Interpolate
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../include/grid.h"
#include "../include/field.h"
#include "../include/interpolate.h"
#include <math.h>

const int OutsideInterp = -2;

struct poly_test {
  poly_test() { BOOST_TEST_MESSAGE( "Setup fixture poly_test"); 
    nr = 50;
    nz = 50;
    Rmin = -5;
    Rmax = 5;
    zmin = -3;
    zmax = 3;

    grid = new Grid(Rmin, Rmax, zmin, zmax, nr, nz);
    psi = new Field(*grid);
    inter = new Interpolate(*grid, *psi);
    R0 = 2.301;
    z0 = 1.01;
  }

  ~poly_test(){ BOOST_TEST_MESSAGE( "Teardown fixture poly_test"); 
    delete psi;
    delete grid; 
  }

  int nr, nz;
  double Rmin, Rmax;
  double zmin, zmax;

  double R0, z0;
  Grid * grid;
  Field * psi;
  Interpolate * inter;
};

BOOST_FIXTURE_TEST_SUITE( suite1, poly_test)
  /*!
   * Tests Interpolate::F functions using a flat function = 0
   */
BOOST_AUTO_TEST_CASE (flat_interp) {
  for (int i=0; i < nr; ++i) {
    for (int j=0; j < nz; ++j) {
      psi->f_[i][j] = 0;
    }
  }

  inter->updateInterpolation(R0,z0);
  try {
    double psi_interp = inter->F(R0,z0);
    double psir_interp = inter->F_r(R0,z0);
    double psiz_interp = inter->F_z(R0,z0);
    double psizz_interp = inter->F_zz(R0,z0);
    double psirz_interp = inter->F_rz(R0,z0);
    double psirr_interp = inter->F_rr(R0,z0);
    BOOST_CHECK_EQUAL(0, psi_interp);
    BOOST_CHECK_EQUAL(0, psir_interp);
    BOOST_CHECK_EQUAL(0, psirr_interp);
    BOOST_CHECK_EQUAL(0, psiz_interp);
    BOOST_CHECK_EQUAL(0, psizz_interp);
    BOOST_CHECK_EQUAL(0, psirz_interp);
  }
  catch(int i) {
    if (i == OutsideInterp) {
      printf("Interpolation outside of current gridcell\n");
    }
  }
}

/*!
 * Tests Interpolate::F functions using a linear function z = r + z
 */
BOOST_AUTO_TEST_CASE (linear_interp) {
  for (int i=0; i < nr; ++i) {
    for (int j=0; j < nz; ++j) {
      double R = grid->R_[i];
      double z = grid->z_[j];
      psi->f_[i][j] = R + z;
    }
  }

  inter->updateInterpolation(R0,z0);
  try {
    double psi_interp = inter->F(R0,z0);
    double psir_interp = inter->F_r(R0,z0);
    double psiz_interp = inter->F_z(R0,z0);
    double psizz_interp = inter->F_zz(R0,z0);
    double psirz_interp = inter->F_rz(R0,z0);
    double psirr_interp = inter->F_rr(R0,z0);
    BOOST_CHECK_CLOSE(R0 + z0, psi_interp, 0.0001);
    BOOST_CHECK_CLOSE(1.0, psir_interp, 0.00001);
    BOOST_CHECK_SMALL(psirr_interp, 0.00001);
    BOOST_CHECK_CLOSE(1.0, psiz_interp, 0.00001);
    BOOST_CHECK_SMALL(psizz_interp, 0.00001);
    BOOST_CHECK_SMALL(psirz_interp, 0.00001);
  }
  catch(int i) {
    if (i == OutsideInterp) {
      printf("Interpolation outside of current gridcell\n");
    }
  }
}

/*!
 * Tests Interpolate::F functions using paraboloid
 */
BOOST_AUTO_TEST_CASE (Paraboloid_interp) {
  for (int i=0; i < nr; ++i) {
    for (int j=0; j < nz; ++j) {
      psi->f_[i][j] = grid->R_[i]*grid->R_[i] + grid->z_[j]*grid->z_[j];
    }
  }

  inter->updateInterpolation(0,0);

  try {
    double psi_interp = inter->F(0.05,0.05);
    double psi2_interp = inter->F(0.01,0.01);
    double psir2_interp = inter->F_r(0.01, 0.01);
    double psiz_interp = inter->F_z(0.05, 0.05);
    double psiz2_interp = inter->F_z(0.01, 0.01);
    double psizz_interp = inter->F_zz(0.01, 0.01);
    double psirr_interp = inter->F_rr(0.01, 0.01);
    BOOST_CHECK_CLOSE(0.005, psi_interp, .5);
    BOOST_CHECK_CLOSE(0.0002, psi2_interp, .5);
    BOOST_CHECK_CLOSE(0.02, psir2_interp, .5);
    BOOST_CHECK_CLOSE(0.1, psiz_interp, .5);
    BOOST_CHECK_CLOSE(0.02, psiz2_interp, .5);
    BOOST_CHECK_CLOSE(2, psizz_interp, .005);
    BOOST_CHECK_CLOSE(2, psirr_interp, .005);
  }
  catch(int i) {
    if (i == OutsideInterp) {
      printf("Interpolation outside of current gridcell\n");
    }
  }
}

/*!
 * Tests Interpolate::F functions using a bicubic function
 */
BOOST_AUTO_TEST_CASE (cubic_interp) {
  for (int i=0; i < nr; ++i) {
    for (int j=0; j < nz; ++j) {
      double R = grid->R_[i];
      double z = grid->z_[j];
      psi->f_[i][j] = R*R*R + 2*z*z*z + 3*R*z*z + 4*R*R*z*z;
    }
  }
  
  inter->updateInterpolation(R0,z0);
  try {
    double psi_interp = inter->F(R0,z0);
    double psir_interp = inter->F_r(R0,z0);
    double psiz_interp = inter->F_z(R0,z0);
    double psizz_interp = inter->F_zz(R0,z0);
    double psirz_interp = inter->F_rz(R0,z0);
    double psirr_interp = inter->F_rr(R0,z0);
    BOOST_CHECK_CLOSE(R0*R0*R0+2*z0*z0*z0+3*R0*z0*z0 + 4*R0*R0*z0*z0, psi_interp, .05);
    BOOST_CHECK_CLOSE(3*R0*R0 + 8*R0*z0*z0 + 3*z0*z0, psir_interp, .05);
    BOOST_CHECK_CLOSE(6*R0 + 8*z0*z0, psirr_interp, .05);
    BOOST_CHECK_CLOSE(6*z0*z0+ 8*R0*R0*z0 + 6*R0*z0, psiz_interp, .05);
    BOOST_CHECK_CLOSE(12*z0 + 8*R0*R0 + 6*R0, psizz_interp, .05);
    BOOST_CHECK_CLOSE( 16*R0*z0 + 6*z0, psirz_interp, .05);
  }
  catch(int i) {
    if (i == OutsideInterp) {
      printf("Interpolation outside of current gridcell\n");
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
