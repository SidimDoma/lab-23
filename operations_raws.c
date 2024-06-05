#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "operations_raws.h"

struct matrix *multiply_row(struct matrix* matrix, double scalar, double row) {
    for (size_t i = 0; i < matrix->cols; i++) {
        double key = *get_element(matrix, row, i);
        add_element(matrix, row, i, key * scalar);
    }
    return matrix;
}

struct matrix *add_strings(struct matrix* matrix, size_t row1, size_t row2) { //прибавим к первому второй
    for (size_t i = 0; i < matrix->cols; i++) {
        double key = *get_element(matrix, row1, i) + *get_element(matrix, row2, i);
        add_element(matrix, row1, i, key);
    }
    return matrix;
}

struct matrix *permute_string(struct matrix* matrix, size_t row1, size_t row2) {
    for (size_t i = 0; i < matrix->cols; i++) {
        double key1 = *get_element(matrix, row1, i);
        double key2 = *get_element(matrix, row2, i);
        add_element(matrix, row1, i, key2);
        add_element(matrix, row2, i, key1);
    }
    return matrix;
}

double matrix_norm(const struct matrix *m) {
    double max = 0.0;
    for (size_t i = 0; i < m->rows; i++) {
        for (size_t j = 0; j < m->cols; j++) {
            double abs_value = fabs(*get_element(m, i, j));
            if (abs_value > max) {
                max = abs_value;
            }
        }
    }
    return max;
}

double matrix_determinant(struct matrix *m) {
    if (m->rows != m->cols) {
        return 0; // Not a square matrix
    }
    if (m->rows == 1) {
        return *get_element(m, 0, 0);
    }
    double det = 0;
    for (size_t j = 0; j < m->cols; ++j) {
        struct matrix *sub = submatrix(m, 0, j);
        if (!sub) {
            return 0; // Submatrix allocation failed
        }
        det += pow(-1, j+1) * *get_element(m, 0, j) * matrix_determinant(sub);
        matrix_free(sub);
    }
    return det;
}
