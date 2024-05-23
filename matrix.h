#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED
#include <stddef.h>

struct matrix {
    double *data;
    size_t rows;
    size_t cols;
};

struct matrix* matrix_alloc(size_t rows, size_t cols);

void add_element(struct matrix *m, size_t i, size_t j, double value);

double *get_element(const struct matrix *m, size_t i, size_t j);

struct matrix* matrix_clone(const struct matrix* m);

struct matrix* matrix_resize(struct matrix* m, size_t rows, size_t cols);

void matrix_print(const struct matrix* matrix);

void matrix_free(struct matrix* matrix);

struct matrix* zero_matrix_set(struct matrix* matrix);

struct matrix* zero_matrix_alloc(size_t rows, size_t cols);

struct matrix *matrix_alloc_square(size_t h);

struct matrix *matrix_alloc_square_unit(size_t h);

struct matrix* matrix_assign(struct matrix* m1, struct matrix* m2);
#endif // MATRIX_H_INCLUDED
