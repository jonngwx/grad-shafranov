/*!
 * @file critical.cc
 * @brief Implementation of class Critical.
 */
#ifndef CRITICAL_H
#define CRITICAL_H
#include "field.h"
#include "grid.h"
#include "interpolate.h"
#include <vector>
#include "tsv_reader.h"

/*! 
 * @brief Finds 'critical points' like the magnetic axis, limiter points, and X points.
 *
 * These points are minima, maximima, and saddle points of \f$\psi\f$.
 */
class Critical {
public:
    Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, Table &limiters, double R_stag_up, double z_stag_up, double R_stag_down, double z_stag_down, double R0, double z0);
    ~Critical();
    
    /*!
     * @brief returns dr, dz to progress toward critical point search in psi
     * @param[in] r MEANS SOMETHING
     * @param[in] z MEANS SOMETHING
     * @param[out] dr MEANS SOMETHING
     * @param[out] dz MEANS SOMETHING
     * @param
     */
    void Psi_search(double r, double z, double *dr, double *dz);
    /*!
     * @brief Perform search for critical points of limiter beginning 
     * with initial guess r, z
     * @param[in] r MEANS SOMETHING
     * @param[in] z MEANS SOMETHING
     * @param rcrit MEANS SOMETHING
     * @param zcrit MEANS SOMETHING
     * @param Psi_min MEANS SOMETHING
     */
    void Psi_magnetic(double r, double z,double *rcrit, double *zcrit, double *Psi_min);

    /*!
     * @brief MEANS SOMETHING
     * @return MEANS SOMETHING
     */
    double Psi_limiter();
    /*! 
     * @brief Performs critical point search; updates Psi_i and Psi_o
     */
    void update();

    /*!
     * @brief MEANS SOMETHING (why are these &references and *pointers above in Psi_magnetic? ... sloppy)
     * @param r MEANS SOMETHING
     * @param z MEANS SOMETHING
     * @return MEANS SOMETHING
     */
    bool find_saddle(double &r, double &z);
private:
    Field &Psi_;
    const int max_iter; //!< maximum number of iterations for critical search
    Grid &Grid_; //!< Underlying grid
    Interpolate Inter_;
    const double epsilon; //!< convergence criterion for critical point search
    double Psi_stag_up; //!< Psi at first limiter
    double Psi_stag_down; //!< Psi at second limiter
    double R0; //!< R coordinate of magnetic axis
    double z0; //!< Z coordinate of magnetic axis
    
    const double R_stag_up_orig;
    const double z_stag_up_orig;
    double R_stag_up; //!< R coordinate of stagnation point
    double z_stag_up; //!< Z coordinate of limiter
    
    const double R_stag_down_orig;
    const double z_stag_down_orig;
    double R_stag_down;
    double z_stag_down;

    Table &limiters_;
    std::vector<double> Psi_phys_lim;
};

#endif
