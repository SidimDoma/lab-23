#ifndef OPERATIONS_H_INCLUDED
#define OPERATIONS_H_INCLUDED

#include "matrix.h"

struct matrix *matrix_sum(const struct matrix *A, const struct matrix *B, struct matrix *res);

struct matrix *multiply_matrix_by_scalar(struct matrix *m, double scalar);

struct matrix *matrix_multiplication(const struct matrix *A, const struct matrix *B);

struct matrix *matrix_transpose_simple(struct matrix *matrix);

struct matrix* matrix_transpose(const struct matrix* matrix);

struct matrix *submatrix(struct matrix *m, size_t num_row, size_t num_col);

#endif // OPERATIONS_H_INCLUDED
