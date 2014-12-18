#ifndef CRITICAL_H
#define CRITICAL_H
#include "field.h"
#include "grid.h"

class Critical {
public:
    Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2, double R0, double z0);
    ~Critical();
    /*!
     * @brief Calculates alpha, beta grid for interpolation
     */
    void interpolate();
    /*!
     * @brief Find alpha for given position
     */
    double cell_alpha(double r, double z);
    /*!
     * @brief Find beta for given position
     */
    double cell_beta(double r, double z);
    /*!
     * @brief Interpolated Psi defined for all r, z
     */
    double Psi_interp(double r, double z);
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
    const int max_iter;
    Grid &Grid_;
    const double epsilon;
    double z_limiter1;
    double z_limiter2;
    // Critical points - updated with every critical search
    double R0; // R coordinate of magnetic axis
    double z0; // Z coordinate of magnetic axis
    double Rl; // R coordinate of limiter
    double zl; // Z coordinate of limiter
    double **alpha;
    double **beta;
};

#endif
