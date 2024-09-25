#ifndef GMPECC_H
#define GMPECC_H

#include <gmp.h>

struct Point
{
    mpz_t x;
    mpz_t y;
};

struct Elliptic_Curve
{
    mpz_t p;
    mpz_t n;
};

extern Elliptic_Curve EC;
extern Point G;
extern Point DoublingG[256];

void Point_Doubling(Point *P, Point *R);
void Point_Addition(Point *P, Point *Q, Point *R);
void Scalar_Multiplication(Point P, Point *R, mpz_t m);
void Point_Negation(Point *A, Point *S);
void init_doublingG(Point *P);

#endif // GMPECC_H