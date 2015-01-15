/*!
 * @file critical.h
 * @author Elizabeth J. Paul
 * @brief Header file for the Critical class.
 */

#ifndef CRITICAL_H
#define CRITICAL_H
#include "field.h"
#include "grid.h"
#include "interpolate.h"
#include <vector>
#include "tsv_reader.h"

/*! 
 * @brief Contains methods to search for the critical points of Psi and update the
 * limiters and magnetic axis.
 */
class Critical {
public:
    Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, Table &limiters, double R_stag_up, double z_stag_up, double R_stag_down, double z_stag_down, double R0, double z0);
    ~Critical();
    /*!
     * @brief returns dr, dz to progress toward critical point search in psi
     * @param[in] r Current r coordinate in critical point search
     * @param[in] z Current z coordinate in critical point search
     * @param[out] dr Update in r direction to progress toward critical point
     * @param[out] dz Update in z direction to progress toward critical point
     */
    void Psi_search(double r, double z, double *dr, double *dz);
    /*!
     * @brief Perform search for critical points of limiter beginning
     * with initial guess r, z
     * @param[in] r Initial guess for magnetic axis r coordinate.
     * @param[in] z Initial guess for magnetic axis z coordinate.
     * @param rcrit R coordinate of determined magnetic axis point.
     * @param zcrit Z coordinate of determined magnetic axis point.
     * @param Psi_min Value of Psi at rcrit, zcrit.
     */
    void Psi_magnetic(double r, double z,double *rcrit, double *zcrit, double *Psi_min);
    /*!
     * @brief Performs search for limitet point using critical point search and 
     * comparison with physical limiter points
     * @return Value of Psi at innermost limiter point, either a critical point or value
     * at physical limiter.
     */
    double Psi_limiter();
    /*! 
     * @brief Performs critical point search; updates Psi_i and Psi_o
     */
    void update();
    /*!
     * @brief Finds saddle point using critical point search starting with r and z.
     *
     * Returns false if search lands outside of horizontal limiters or grid
     * boundaries.
     *
     * @param r Current r coordinate.
     * @param z Current z coordinate.
     * @return True if search critical point corresponds to saddle point, False otherwise or if outside boundaries
     */
    bool find_saddle(double &r, double &z);
private:
    Field &Psi_;
    const int max_iter; /** < maximum number of iterations for critical search */
    Grid &Grid_; /** < Underlying grid */
    Interpolate Inter_; /** < Instance of interpolation class */
    const double epsilon; /** < convergence criterion for critical point search */
    double Psi_stag_up; /** < Psi at first limiter */
    double Psi_stag_down; /** < Psi at second limiter */
    double R0; /** < R coordinate of magnetic axis */
    double z0; /** < Z coordinate of magnetic axis */
    const double R_stag_up_orig;
    const double z_stag_up_orig;
    double R_stag_up; /** < R coordinate of stagnation point */
    double z_stag_up; /** < Z coordinate of limiter */
    const double R_stag_down_orig;
    const double z_stag_down_orig;
    double R_stag_down;
    double z_stag_down;

    Table &limiters_;
    std::vector<double> Psi_phys_lim;
};

#endif
