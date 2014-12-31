#ifndef CRITICAL_H
#define CRITICAL_H
#include "field.h"
#include "grid.h"
#include <vector>

class Critical {
public:
    Critical(Grid &GridS, Field &Psi, int max_iter, double epsilon, double z_limiter1, double z_limiter2, double R0, double z0);
    ~Critical();
    
    double Psi_interp(double r, double z);
    double Psir_interp(double r, double z);
    double Psirr_interp(double r, double z);
    double Psirz_interp(double r, double z);
    double Psizz_interp(double r, double z);
    double Psiz_interp(double r, double z);
    
    /*!
     * @brief returns dr, dz to progress toward critical point search 
     * in Psi
     */
    void Psi_search(double r, double z, double *dr, double *dz);
    /*!
     * @brief Perform search for critical points of limiter beginning 
     * with initial guess r, z
     */
    void updateCoefficients();
    void updateP(double r, double z);
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
    // Psi at limiters
    double Psi_lim1;
    double Psi_lim2;
    // Critical points - updated with every critical search
    double R0; // R coordinate of magnetic axis
    double z0; // Z coordinate of magnetic axis
    double Rl; // R coordinate of limiter
    double zl; // Z coordinate of limiter
    double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33; // Coefficients for interpolation
    // Array of 4 points corresponding to surrounding used for interpolation
    std::vector<std::vector<double>> P;
};

#endif
