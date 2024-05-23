#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "exponent.h"

struct matrix* matrix_power(const struct matrix* A, int power, struct matrix* res) {
    if (power == 0) {
        matrix_free(res);
        return matrix_alloc_square_unit(A->rows);
    }
    else if (power == 1) {
        matrix_free(res);
        return matrix_clone(A);
    }
    else {
        struct matrix* temp = matrix_clone(A);

        for (int i = 2; i <= power; i++) {
            matrix_free(res);
            res = matrix_multiplication(temp, A);
            matrix_assign(temp, res);
        }

        matrix_free(temp);
        return res;
    }
}

struct matrix* matrix_exponential(const struct matrix *A, double eps) {
    struct matrix *result = matrix_alloc_square_unit(A->rows);

    struct matrix *temp = matrix_clone(A);

    double n_factorial = 1;
    double norm_term = 0;
    for (int n = 1; ; n++) {
        matrix_power(A, n, temp);
        n_factorial *= n;
        multiply_matrix_by_scalar(temp, 1.0 / n_factorial);
        matrix_sum(result, temp, result);
        norm_term = matrix_norm(temp);
        if (norm_term < eps) {
            break;
        }
    }
    matrix_free(temp);
    return result;
}
