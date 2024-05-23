#ifndef GAUSS_H_INCLUDED
#define GAUSS_H_INCLUDED

#include "matrix.h"
#include "operations_raws.h"
#include "operations.h"

struct matrix* add_multiply_row_by_scalar(struct matrix* m, size_t k, size_t t, double scalar);

struct matrix* inverse_matrix(struct matrix *m);

struct matrix* solve_system(struct matrix* A, struct matrix* B, struct matrix* x);

#endif // GAUSS_H_INCLUDED
