#ifndef CRITICAL_H
#define CRITICAL_H

class Critical {
public:
    Critical(const Grid &GridS, const Field &Psi);
    ~Critical();
    // Calculates Psi_r, Psi_z, Psi_rr, Psi_zz, Psi_rz for interpolation
    void interpolate();
    // Interpolated Psi defined for all r, z
    void Psi_interp(double r, double z) {
private:
    const Field &Psi;
    const Grid &GridS;
    double **Psi_r;
    double **Psi_z;
    double **Psi_rr;
    double **Psi_zz;
    double **Psi_rz;
    double **alpha;
    double **beta;
};

#endif