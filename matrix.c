#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
/*
struct matrix {
    double *data;
    size_t rows;
    size_t cols;
};*/

struct matrix *matrix_alloc(size_t h, size_t w) {
    struct matrix *m = malloc(sizeof(struct matrix));
    double *new_data;
    new_data = malloc(h * w * sizeof(double));
    if (new_data && m) {
        m->rows = h;
        m->cols = w;
        m->data = new_data;
        return m;
    }
    free(new_data);
    free(m);
    return NULL;
}


void matrix_free(struct matrix *m) {
    free(m->data);                      // Freeing the memory allocated for the matrix data.
    free(m);                            // Freeing the memory allocated for the matrix structure.
}


void add_element(struct matrix *m, size_t i, size_t j, double value) {
    m->data[i * m->cols + j] = value;
}

double *get_element(const struct matrix *m, size_t i, size_t j) {             // Gets values from matrix sets
    if (i < m->rows && j < m->cols) {
        return &m->data[i * m->cols + j];
    }
    return 0;
}


struct matrix *matrix_clone(const struct matrix *m) {

    struct matrix *clone = matrix_alloc(m->rows, m->cols);

    if (clone) {
        for (size_t i = 0; i < m->rows; ++i) {
            for (size_t j = 0; j < m->cols; ++j) {
                add_element(clone, i, j, *get_element(m, i, j));
            }
        }
        return clone;
    }
    return NULL;
}


struct matrix* matrix_resize(struct matrix* m, size_t rows, size_t cols) {             // Changes size of matrix
    double* newData = (double*) realloc(m->data, rows * cols * sizeof(double));

    if (newData == NULL) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    m->rows = rows;
    m->cols = cols;
    m->data = newData;

    return m;
}


void print_matrix(const struct matrix* matrix) {                                 // Prints matrix
    for (size_t i = 0; i < matrix->rows; i++) {
        for (size_t j = 0; j < matrix->cols; j++) {
            printf("%-13lf", *get_element(matrix, i, j));
        }
        printf("\n");
    }
}


struct matrix* zero_matrix_set(struct matrix* matrix) {
    size_t i, j;

    // Setting each element of the matrix to zero.
    for (i = 0; i < matrix->rows; i++) {
        for (j = 0; j < matrix->cols; j++) {
            matrix->data[i * matrix->cols + j] = 0.0;
        }
    }

    return matrix;  // Returning the matrix with all elements set to zero.
}

struct matrix* zero_matrix_alloc(size_t rows, size_t cols) {
    struct matrix* matrix = matrix_alloc(rows, cols);  // Allocating memory for the matrix.
    zero_matrix_set(matrix);                    // Setting all elements to zero.
    return matrix;                              // Returning the allocated matrix.
}



struct matrix *matrix_alloc_square(size_t h) {
    return matrix_alloc(h, h);
}

struct matrix *matrix_alloc_square_unit(size_t h) {
    struct matrix *m = matrix_alloc_square(h);
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            if (i == j) {
                add_element(m, i, j, 1);
            } else {
                add_element(m, i, j, 0);
            }
        }
    }
    return m;
}

struct matrix* matrix_assign(struct matrix* m1, struct matrix* m2) {
    if (m1->rows != m2->rows || m1->cols != m2->cols) {
        return NULL;
    }
    for (size_t i = 0; i < m1->rows; ++i) {
        for (size_t j = 0; j < m1->cols; ++j) {
            double key = *get_element(m2, i, j);
            add_element(m1, i, j, key);
        }
    }
    return m1;
}
