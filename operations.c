#include <stdio.h>
#include <stdlib.h>
#include "operations.h"

struct matrix *matrix_sum(const struct matrix *A, const struct matrix *B, struct matrix *res) {
    if (A->rows != B->rows || A->cols != B->cols) {
        return NULL;
    }
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            *get_element(res, i, j) = *get_element(A, i, j) + *get_element(B, i, j);
        }
    }
    return res;
}

struct matrix *multiply_matrix_by_scalar(struct matrix *m, double scalar) {
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            *get_element(m, i, j) *= scalar;
        }
    }
    return m;
}

struct matrix *matrix_multiplication(const struct matrix *A, const struct matrix *B) {
    if (A->cols != B->rows) {
        return NULL;
    }
    struct matrix *result = matrix_alloc(A->rows, B->cols);
    if (result == NULL) {
        return NULL; // Не удалось выделить память для результирующей матрицы
    }
    for (size_t i = 0; i < A->rows; ++i) {
        for (size_t j = 0; j < B->cols; ++j) {
            double sum = 0;
            for (size_t k = 0; k < A->cols; ++k) {
                sum += *get_element(A, i, k) * *get_element(B, k, j);
            }
            add_element(result, i, j, sum);
        }
    }
    return result;
}


struct matrix *matrix_transpose_simple(struct matrix *matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = i + 1; j < matrix->cols; j++) {
             double key1 = *get_element(matrix, i, j);
             double key2 = *get_element(matrix, j, i);
             add_element(matrix, i, j, key2);
             add_element(matrix, j, i, key1);
        }
    }
    return matrix;
}

struct matrix *matrix_transpose(const struct matrix *matrix) {
    struct matrix *result = matrix_alloc(matrix->cols, matrix->rows);

    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            double key = *get_element(matrix, j, i);
            add_element(result, i, j, key);;
        }
    }

    return result;
}

struct matrix *submatrix(struct matrix *m, size_t num_row, size_t num_col) {
    if (num_row >= m->rows || num_col >= m->cols) {
        return NULL; // Invalid row or column index
    }
    struct matrix *result = matrix_alloc(m->rows - 1, m->cols - 1);
    if (!result) {
        return NULL; // Allocation failed
    }
    size_t r = 0, c = 0;
    for (size_t i = 0; i < m->rows; ++i) {
        if (i == num_row) {
            continue;
        }
        for (size_t j = 0; j < m->cols; ++j) {
            if (j == num_col) {
                continue;
            }
            add_element(result, r, c, *get_element(m, i, j));
            c++;
        }
        r++;
        c = 0;
    }
    return result;
}
