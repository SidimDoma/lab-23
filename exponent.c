#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "exponent.h"

struct matrix* matrix_power(const struct matrix* A, int power, struct matrix* res) {
    if (power == 0) {
        matrix_free(res);
        return matrix_alloc_square_unit(A->rows);
    } else if (power == 1) {
        matrix_free(res);
        return matrix_clone(A);
    } else {
        struct matrix* temp = matrix_clone(A);
        struct matrix* temp_result = matrix_clone(A);
        matrix_free(res);
        res = matrix_clone(A);

        for (int i = 2; i <= power; i++) {
            matrix_free(temp_result);
            temp_result = matrix_multiplication(res, A);
            matrix_free(res);
            res = matrix_clone(temp_result);
        }

        matrix_free(temp);
        matrix_free(temp_result);
        return res;
    }
}

struct matrix* matrix_exponential(const struct matrix *A, double eps) {
    struct matrix *result = matrix_alloc_square_unit(A->rows);
    struct matrix *term = matrix_clone(A);
    struct matrix *temp;

    double n_factorial = 1;
    double norm_term = matrix_norm(term);

    for (int n = 1; norm_term >= eps; n++) {
        multiply_matrix_by_scalar(term, 1.0 / n);
        matrix_sum(result, term, result);

        n_factorial *= n;

        temp = matrix_multiplication(term, A);
        matrix_free(term);
        term = temp;

        norm_term = matrix_norm(term);
    }

    matrix_free(term);
    return result;
}
