#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EPS 1e-7
#define INCORRECT_MATRIX 1
#define CALCULATION_ERROR 2
#define CORRECT_MATRIX 0
#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double** matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t* result);
void s21_remove_matrix(matrix_t* A);
int s21_eq_matrix(matrix_t* A, matrix_t* B);
int s21_sum_matrix(matrix_t* A, matrix_t* B, matrix_t* result);
int s21_sub_matrix(matrix_t* A, matrix_t* B, matrix_t* result);
int s21_mult_number(matrix_t* A, double number, matrix_t* result);
int s21_mult_matrix(matrix_t* A, matrix_t* B, matrix_t* result);
int s21_transpose(matrix_t* A, matrix_t* result);
int s21_calc_complements(matrix_t* A, matrix_t* result);
int s21_determinant(matrix_t* A, double* result);
int s21_inverse_matrix(matrix_t* A, matrix_t* result);
matrix_t s21_minor(matrix_t* A, int rows, int columns);
int eq_calc(matrix_t* A, matrix_t* B);
int eq(matrix_t* A);
void zetloc(matrix_t* A);

#endif  // SRC_S21_MATRIX_H_
