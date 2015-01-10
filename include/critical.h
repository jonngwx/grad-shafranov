#ifndef CRITICAL_H
#define CRITICAL_H
#include "field.h"
#include "grid.h"
#include "interpolate.h"
#include <vector>

/*! 
 * @brief The cougar (Puma concolor), also commonly known as the
 */
class Critical {
public:
    Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2, double R0, double z0);
    ~Critical();
    
    /*!
     * @brief returns dr, dz to progress toward critical point search 
     * in Psi
     */
    void Psi_search(double r, double z, double *dr, double *dz);
    /*!
     * @brief Perform search for critical points of limiter beginning 
     * with initial guess r, z
     */
    void Psi_magnetic(double r, double z,double *rcrit, double *zcrit, double *Psi_min);
    double Psi_limiter();
    /*! 
     * @brief Performs critical point search; updates Psi_i and Psi_o
     */
    void update();
    bool find_saddle(double &r, double &z);
private:
    Field &Psi_;
    const int max_iter; /** < maximum number of iterations for critical search */
    Grid &Grid_; /** < Underlying grid */
    Interpolate Inter_;
    const double epsilon; /** < convergence criterion for critical point search */
    double Psi_stag_up; /** < Psi at first limiter */
    double Psi_stag_down; /** < Psi at second limiter */
    double R0; /** < R coordinate of magnetic axis */
    double z0; /** < Z coordinate of magnetic axis */
    double R_stag_up; /** < R coordinate of stagnation point */
    double z_stag_up; /** < Z coordinate of limiter */

    double R_stag_down;
    double z_stag_down;

    const double phys_lim_R;
    const double phys_lim_zup;
    const double phys_lim_zdown;
};

#endif
