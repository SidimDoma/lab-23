#ifndef EXPONENT_H_INCLUDED
#define EXPONENT_H_INCLUDED

#include "matrix.h"
#include "operations_raws.h"
#include "operations.h"

struct matrix *matrix_power(const struct matrix *A, int power, struct matrix *res);

struct matrix* matrix_exponential(const struct matrix *A, double eps);


#endif // EXPONENT_H_INCLUDED
