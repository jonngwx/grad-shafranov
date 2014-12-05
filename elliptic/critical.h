#ifndef CRITICAL_H
#define CRITICAL_H

class Critical {
public:
    Critical(const Grid &GridS, const Field &Psi, int max_iter, double epsilon);
    ~Critical();
    // Calculates alpha, beta grid for interpolation
    void interpolate();
    // Interpolated Psi defined for all r, z
    void Psi_interp(double r, double z);
    // Find alpha for given position
    double cell_alpha(double r, double z);
    // Find beta for given position
    double cell_beta(double r, double z);
    // Interpolated Psi defined for all r, z
    void Psi_interp(double r, double z);
    // returns dr, dz to progress toward critical point search in Psi
    void Psi_search(double r, double z, double *dr, double *dz);
    // Perform search for critical points beginning with initial
    // guess r, z
    void Psi_critical(double r, double z,double *rcrit, double *zcrit)
private:
    const double epsilon;
    const int max_iter;
    const Field &Psi_;
    const Grid &Grid_;
    double **alpha;
    double **beta;
};

#endif