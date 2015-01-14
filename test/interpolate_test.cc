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
  poly_test() {
    BOOST_TEST_MESSAGE( "Setup fixture poly_test"); 
    nr = 11;
    nz = 11;
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

  ~poly_test() {
    BOOST_TEST_MESSAGE( "Teardown fixture poly_test"); 
    // inter->PrintAmnCoefficients();
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
 * Tests Interpolate::F functions using a general function paraboloid
 */
BOOST_AUTO_TEST_CASE (Paraboloid_interp) {
  for (int i=0; i < nr; ++i) {
    for (int j=0; j < nz; ++j) {
      double R = grid->R_[i];
      double z = grid->z_[j];
      psi->f_[i][j] = R*R + z*z - 2*R*z + 1;
    }
  }

  inter->updateInterpolation(R0,z0);

  try {
    double psi_interp   = inter->F   (R0,z0);
    double psir_interp  = inter->F_r (R0,z0);
    double psiz_interp  = inter->F_z (R0,z0);
    double psizz_interp = inter->F_zz(R0,z0);
    double psirr_interp = inter->F_rr(R0,z0);
    double psirz_interp = inter->F_rz(R0,z0);
    BOOST_CHECK_CLOSE(R0*R0 + z0*z0 - 2*R0*z0 + 1, psi_interp, 0.001);
    BOOST_CHECK_CLOSE(2*R0 - 2*z0, psir_interp, 0.0001);
    BOOST_CHECK_CLOSE(2*z0 - 2*R0, psiz_interp, 0.0001);
    BOOST_CHECK_CLOSE(2, psirr_interp, 0.0001);
    BOOST_CHECK_CLOSE(2, psizz_interp, 0.0001);
    BOOST_CHECK_CLOSE(-2, psirz_interp, 0.0001);
  }
  catch(int i) {
    if (i == OutsideInterp) {
      printf("Interpolation outside of current gridcell\n");
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
