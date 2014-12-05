#ifndef RHSFUNC_H_
#define RHSFUNC_H_

class RHSfunc {
  public:
    RHSfunc(string pgtype, PGData data); 
    ~RHSfunc();
    double eval(double var);
    double eval_prime(double var);
  private:
    /* to be determined */

};

#endif // RHSFUNC_H_

