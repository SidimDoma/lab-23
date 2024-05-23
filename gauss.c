#include <stdio.h>
#include <stdlib.h>
#include "gauss.h"

struct matrix* add_multiply_row_by_scalar(struct matrix* m, size_t k, size_t t, double scalar){
    for (size_t i = 0; i < m->cols; ++i) {
        add_element(m, k, i, *get_element(m, k, i) + scalar * *get_element(m, t, i));
    }
    return m;
}

struct matrix* inverse_matrix(struct matrix* m) {
    if (m->rows != m->cols) {
        return NULL;
    }

    struct matrix* result = matrix_alloc_square_unit(m->rows);
    if (result == NULL) {
        return NULL;
    }

    struct matrix* temp = matrix_clone(m);
    if (temp == NULL) {
        matrix_free(result);
        return NULL;
    }

    for (size_t i = 0; i < temp->rows; ++i) {
        size_t row = i;
        while (row < temp->rows && *get_element(temp, row, i) == 0) {
            row++;
        }
        if (row == temp->rows) {
            matrix_free(result);
            matrix_free(temp);
            return NULL;
        }
        permute_string(temp, i, row);
        permute_string(result, i, row);

        double element = *get_element(temp, i, i);
        multiply_row(temp, 1 / element, i);
        multiply_row(result, 1 / element, i);

        for (size_t j = 0; j < temp->rows; ++j) {
            if (j != i) {
                double scalar = -*get_element(temp, j, i);

                for (size_t k = 0; k < temp->cols; ++k) {
                    add_element(temp, j, k, *get_element(temp, j, k) + scalar * *get_element(temp, i, k));
                    add_element(result, j, k, *get_element(result, j, k) + scalar * *get_element(result, i, k));
                }
            }
        }
    }

    matrix_free(temp);
    return result;
}

struct matrix* solve_system(struct matrix* A, struct matrix* B, struct matrix* x) {

    if (A->rows != A->cols) {
        return NULL;
    }
    struct matrix *A_inv = inverse_matrix(A);
    if (A_inv == NULL) {
        return NULL;
    }

    x = matrix_multiplication(A_inv, B);
    matrix_free(A_inv);
    return x;
}

