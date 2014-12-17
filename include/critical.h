#ifndef CRITICAL_H
#define CRITICAL_H
#include "field.h"
#include "grid.h"

class Critical {
public:
    Critical(const Grid &GridS, const Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2);
    ~Critical();
    /*!
     * @brief Calculates alpha, beta grid for interpolation
     */
    void interpolate();
    /*!
     * @brief Interpolated Psi defined for all r, z
     */
    void Psi_interp(double r, double z);
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
    void Psi_interp(double r, double z);
    /*!
     * @brief returns dr, dz to progress toward critical point search 
     * in Psi
     */
    void Psi_search(double r, double z, double *dr, double *dz);
    /*!
     * @brief Perform search for critical points of limiter beginning 
     * with initial guess r, z
     */
    void Psi_limiter(double r, double z,double *rcrit, double *zcrit, double *Psi_min);
    /*! 
     * @brief Performs critical point search; updates Psi_i and Psi_o
     */
    void critical();
private:
    const double epsilon;
    const int max_iter;
    const double z_limiter1;
    const double z_limiter2;
    // Critical points - updated with every critical search
    double R0; // R coordinate of magnetic axis
    double z0; // Z coordinate of magnetic axis
    double Rl; // R coordinate of limiter
    double zl; // Z coordinate of limiter
    const Field &Psi_;
    const Grid &Grid_;
    double **alpha;
    double **beta;
};

#endif