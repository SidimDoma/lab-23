#ifndef OPERATIONS_RAWS_H_INCLUDED
#define OPERATIONS_RAWS_H_INCLUDED

#include "matrix.h"
#include "operations.h"

struct matrix *multiply_row(struct matrix* matrix, double scalar, double row);

struct matrix *add_strings(struct matrix* matrix, size_t row1, size_t row2);

struct matrix *permute_string(struct matrix* matrix, size_t row1, size_t row2);

double matrix_norm(const struct matrix *m);

double matrix_determinant(struct matrix *m);

#endif // OPERATIONS_RAWS_H_INCLUDED
