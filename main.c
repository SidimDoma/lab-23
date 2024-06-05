#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "matrix.h"
#include "operations.h"
#include "operations_raws.h"
#include "gauss.h"
#include "exponent.h"
#include <locale.h>

int main() {
    setlocale(LC_CTYPE, "Russian");
    struct matrix *m = matrix_alloc(3, 3);
    add_element(m, 0, 0, 1);
    add_element(m, 0, 1, 2);
    add_element(m, 0, 2, 3);

    add_element(m, 1, 0, 4);
    add_element(m, 1, 1, 5);
    add_element(m, 1, 2, 6);

    add_element(m, 2, 0, 7);
    add_element(m, 2, 1, 8);
    add_element(m, 2, 2, 9);

    printf("Matrix M1:\n");
    print_matrix(m);

    struct matrix *m2 = matrix_alloc(3, 3);
    add_element(m2, 0, 0, 2);
    add_element(m2, 0, 1, 4);
    add_element(m2, 0, 2, 6);

    add_element(m2, 1, 0, 8);
    add_element(m2, 1, 1, 15);
    add_element(m2, 1, 2, 12);

    add_element(m2, 2, 0, 14);
    add_element(m2, 2, 1, 16);
    add_element(m2, 2, 2, 20);

    printf("\nMatrix M2:\n");
    print_matrix(m2);

    printf("\nCloning the matrix M1:\n");
    struct matrix *m_clone1 = matrix_clone(m);
    struct matrix *m_clone2 = matrix_clone(m);
    struct matrix *m_clone3 = matrix_clone(m);
    struct matrix *m_clone4 = matrix_clone(m);
    struct matrix *m_clone5 = matrix_clone(m);
    struct matrix *m_clone6 = matrix_clone(m);
    struct matrix *m_clone7 = matrix_clone(m);

    print_matrix(m_clone1);
    matrix_free(m);
    matrix_free(m_clone1);

    struct matrix *m_clone1_2 = matrix_clone(m2);
    struct matrix *m_clone2_2 = matrix_clone(m2);
    struct matrix *m_clone3_2 = matrix_clone(m2);
    struct matrix *m_clone4_2 = matrix_clone(m2);
    struct matrix *m_clone5_2 = matrix_clone(m2);
    struct matrix *m_clone6_2 = matrix_clone(m2);
    struct matrix *m_clone7_2 = matrix_clone(m2);
    struct matrix *m_clone8_2 = matrix_clone(m2);

    printf("\nResizing the cloned matrix to 2x2:\n");
    struct matrix *m_res = matrix_resize(m_clone2, 2, 2);
    print_matrix(m_res);
    matrix_free(m_clone2);

    printf("\nCreating a zero matrix (4x3):\n");
    struct matrix *zero_matrix = zero_matrix_alloc(4, 3);
    print_matrix(zero_matrix);
    matrix_free(zero_matrix);

    printf("\nCreating a unit matrix (4x4):\n");
    struct matrix *unit_matrix = matrix_alloc_square_unit(4);
    print_matrix(unit_matrix);
    matrix_free(unit_matrix);

    printf("\nМатрица м1 со значениями м2:\n");
    struct matrix *assign_matrix = matrix_assign(m_clone3, m2);
    print_matrix(assign_matrix);
    matrix_free(assign_matrix);

    /**Operations*/
    printf("\n М1 matrix added to М2:\n");
    struct matrix *res = matrix_alloc(3, 3);
    struct matrix *sum_matrix = matrix_sum(m_clone4, m_clone1_2, res);
    print_matrix(sum_matrix);
    matrix_free(sum_matrix);

    printf("\nМ1 multiplyed with scalar 5 :\n");
    struct matrix *mult_scal_matrix = multiply_matrix_by_scalar(m_clone5, 5);
    print_matrix(mult_scal_matrix);
    matrix_free(mult_scal_matrix);

    printf("\nМ1 matrix multiplyed with М2 :\n");
    struct matrix *mult_matrix = matrix_multiplication(m_clone6, m_clone2_2);
    print_matrix(mult_matrix);
    matrix_free(mult_matrix);

    printf("\nМ1 matrix transponed :\n");
    struct matrix *tr_matrix = matrix_transpose_simple(m_clone7);
    print_matrix(tr_matrix);
    matrix_free(tr_matrix);

    printf("\nМ2 matrix transponed :\n");
    struct matrix *transpose_matrix = matrix_transpose(m_clone3_2);
    print_matrix(transpose_matrix);
    matrix_free(transpose_matrix);

    printf("\nМ2 matrix without 1 raw and 2 column :\n");
    struct matrix *sub_matrix = submatrix(m_clone3_2, 0, 1);
    print_matrix(sub_matrix);
    matrix_free(sub_matrix);

    printf("\nМ2 matrix where 1 raw multiplyed with 3 :\n");
    struct matrix *rowm_matrix = multiply_row(m_clone4_2, 3, 0);
    print_matrix(rowm_matrix);
    matrix_free(rowm_matrix);

    printf("\nМ2 matrix, with 1 raw added with 2:\n");
    struct matrix *add_matrix = add_strings(m_clone5_2, 0, 1);
    print_matrix(add_matrix);
    matrix_free(add_matrix);

    printf("\nМ2 matix, with 2 и 3 raws swapped:\n");
    struct matrix *per_matrix = permute_string(m_clone6_2, 1, 2);
    print_matrix(per_matrix);
    matrix_free(per_matrix);

    printf("\nМ2 matrix norm: %lf\n", matrix_norm(m_clone7_2));
    matrix_free(m_clone7_2);

    printf("\n Determinant of М2 matrix: %lf\n", matrix_determinant(m_clone8_2));
    matrix_free(m_clone8_2);

    // Create a matrix A
    printf("Enter 3х3 matrix А :\n");
            size_t h1 = 3;
            size_t w1 = 3;

            struct matrix* mA = matrix_alloc(h1, w1);
            if (mA == NULL) {
                return 1;
            }
            for (size_t i = 0; i < h1; i++) {
                for (size_t j = 0; j < w1; j++) {
                    double element;
                    scanf("%lf", &element);
                    add_element(mA, i, j, element);
                }
            }
    // Create matrix B (right-hand side of the system)
    printf("Enter 3х1 matrix B :\n");
            size_t h2 = 3;
            size_t w2 = 1;

            struct matrix* mB = matrix_alloc(h2, w2);
            if (mB == NULL) {
                return 1;
            }
            for (size_t i = 0; i < h2; i++) {
                for (size_t j = 0; j < w2; j++) {
                    double element;
                    scanf("%lf", &element);
                    add_element(mB, i, j, element);
                }
            }
    // Create a matrix for the solution 'x'
    struct matrix* x = matrix_alloc(3, 1);
    printf("\nMatrix A:\n");
    print_matrix(mA);
    printf("\nMatrix B:\n");
    print_matrix(mB);
    printf("\nSolving the system of equations:\n");
    x = solve_system(mA, mB, x);
    print_matrix(x);
    matrix_free(mA);
    matrix_free(mB);
    matrix_free(x);

    struct matrix *m1 = matrix_alloc(3, 3);
    add_element(m1, 0, 0, 1);
    add_element(m1, 0, 1, 0);
    add_element(m1, 0, 2, 0);

    add_element(m1, 1, 0, 0);
    add_element(m1, 1, 1, 2);
    add_element(m1, 1, 2, 0);

    add_element(m1, 2, 0, 0);
    add_element(m1, 2, 1, 0);
    add_element(m1, 2, 2, -1);

    print_matrix(m1);
    struct matrix *exp = matrix_exponential(m1, 0.00001);

    printf("\nМатричная экспонента:\n");
    print_matrix(exp);

    matrix_free(exp);

    return 0;
}
