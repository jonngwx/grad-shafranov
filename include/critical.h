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
    void Psi_limiter(double r, double z,double *rcrit, double *zcrit, double *Psi_min);
    /*! 
     * @brief Performs critical point search; updates Psi_i and Psi_o
     */
    void update();
private:
    Field &Psi_;
    const int max_iter; /** < maximum number of iterations for critical search */
    Grid &Grid_; /** < Underlying grid */
    Interpolate Inter_;
    const double epsilon; /** < convergence criterion for critical point search */
    double z_limiter1; /** < z coordinate of limiter 1 */
    double z_limiter2; /** < z coordinate of limiter 2 */
    double Psi_lim1; /** < Psi at first limiter */
    double Psi_lim2; /** < Psi at second limiter */
    double R0; /** < R coordinate of magnetic axis */
    double z0; /** < Z coordinate of magnetic axis */
    double Rl; /** < R coordinate of limiter */
    double zl; /** < Z coordinate of limiter */

    const double phys_lim_R;
    const double phys_lim_zup;
    const double phys_lim_zdown;
};

#endif
