#ifndef BOUNDARY_H_
#define BOUNDARY_H_

class Boundary {
    public:
        virtual ~Boundary() {}
        virtual int CalcB(Field Psi) = 0;
};

#endif // BOUNDARY_H_
