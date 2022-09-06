#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t* result) {
  int report = 1;
  if (rows > 0 && columns > 0) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = (double**)calloc(rows, sizeof(double*));
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = (double*)calloc(columns, sizeof(double));
    }
    report = 0;
  } else {
    zetloc(result);
  }

  return report;
}
void s21_remove_matrix(matrix_t* A) {
  if (A->columns > 0 && A->rows > 0 && A->matrix) {
    for (int i = 0; i < A->rows; i++) free(A->matrix[i]);
    free(A->matrix);
    A->matrix = NULL;
    A->columns = 0;
    A->rows = 0;
  }
}
int s21_eq_matrix(matrix_t* A, matrix_t* B) {
  int report = FAILURE;
  if (eq_calc(A, B)) {
    report = SUCCESS;
    for (int i = 0; i < A->rows && report == SUCCESS; i++) {
      for (int j = 0; j < A->columns && report == SUCCESS; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) {
          report = FAILURE;
        }
      }
    }
  }
  return report;
}
int s21_sum_matrix(matrix_t* A, matrix_t* B, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  zetloc(result);
  if (eq(A) && eq(B)) {
    if (eq_calc(A, B)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
      report = CORRECT_MATRIX;
    } else {
      report = CALCULATION_ERROR;
    }
  }
  return report;
}
int s21_sub_matrix(matrix_t* A, matrix_t* B, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  zetloc(result);
  if (eq(A) && eq(B)) {
    if (eq_calc(A, B)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
      report = CORRECT_MATRIX;
    } else {
      report = CALCULATION_ERROR;
    }
  }
  return report;
}
int s21_mult_number(matrix_t* A, double number, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  zetloc(result);
  if (eq(A)) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
    report = CORRECT_MATRIX;
  }
  return report;
}
int s21_mult_matrix(matrix_t* A, matrix_t* B, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  zetloc(result);
  zetloc(result);
  if (A->columns > 0 && A->rows > 0 && B->columns > 0 && B->rows > 0) {
    if (A->columns != B->rows) {
      report = CALCULATION_ERROR;
    } else {
      s21_create_matrix(A->rows, B->columns, result);
      matrix_t tempB;
      s21_transpose(B, &tempB);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < tempB.rows; j++) {
          result->matrix[i][j] = 0;
          for (int k = 0; k < tempB.columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * tempB.matrix[j][k];
          }
        }
      }
      s21_remove_matrix(&tempB);
      report = CORRECT_MATRIX;
    }
  }
  return report;
}
int s21_transpose(matrix_t* A, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  zetloc(result);
  if (eq(A)) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
    report = CORRECT_MATRIX;
  }
  return report;
}
int s21_calc_complements(matrix_t* A, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  zetloc(result);
  if (eq(A)) {
    if (A->rows == A->columns) {
      s21_create_matrix(A->rows, A->columns, result);
      if (A->rows == 1) {
        result->matrix[0][0] = 1.0;
      } else {
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            matrix_t matrix_minor = s21_minor(A, i, j);
            double res_determ = 0.0;
            report = s21_determinant(&matrix_minor, &res_determ);
            if (!report) {
              result->matrix[i][j] = pow(-1, i + j) * res_determ;
              report = CORRECT_MATRIX;
            }
            s21_remove_matrix(&matrix_minor);
          }
        }
      }
    } else {
      report = CALCULATION_ERROR;
    }
  }
  return report;
}
int s21_determinant(matrix_t* A, double* result) {
  int report = INCORRECT_MATRIX;
  int inv = 1;
  *result = 0.0;
  if (eq(A)) {
    if (A->columns == A->rows) {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
        report = CORRECT_MATRIX;
      } else if (A->rows == 2) {
        *result = (A->matrix[0][0] * A->matrix[1][1]) -
                  (A->matrix[1][0] * A->matrix[0][1]);
        report = CORRECT_MATRIX;
      } else {
        double tmp = 0.0;
        for (int i = 0; i < A->columns; i++) {
          matrix_t matrix_minor = s21_minor(A, 0, i);
          report = s21_determinant(&matrix_minor, &tmp);
          if (!report) {
            *result += inv * A->matrix[0][i] * tmp;
            inv = -inv;
            s21_remove_matrix(&matrix_minor);
            report = CORRECT_MATRIX;
          }
        }
      }
    } else {
      report = CALCULATION_ERROR;
    }
  }
  return report;
}
int s21_inverse_matrix(matrix_t* A, matrix_t* result) {
  int report = INCORRECT_MATRIX;
  double res = 0.0;
  zetloc(result);
  if (eq(A)) {
    int report_determ = s21_determinant(A, &res);
    if ((!report_determ) && (A->columns == A->rows) && (fabs(res) >= EPS)) {
      matrix_t tmp1;
      matrix_t tmp2;
      s21_calc_complements(A, &tmp1);
      s21_transpose(&tmp1, &tmp2);
      res = s21_mult_number(&tmp2, 1.0 / res, result);
      s21_remove_matrix(&tmp1);
      s21_remove_matrix(&tmp2);
      report = CORRECT_MATRIX;
    } else {
      report = CALCULATION_ERROR;
    }
  }
  return report;
}
matrix_t s21_minor(matrix_t* A, int rows, int columns) {
  matrix_t result;
  s21_create_matrix(A->rows - 1, A->columns - 1, &result);
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      if (i >= rows && j >= columns) {
        result.matrix[i][j] = A->matrix[i + 1][j + 1];
      } else if (i >= rows) {
        result.matrix[i][j] = A->matrix[i + 1][j];
      } else if (j >= columns) {
        result.matrix[i][j] = A->matrix[i][j + 1];
      } else {
        result.matrix[i][j] = A->matrix[i][j];
      }
    }
  }
  return result;
}
int eq(matrix_t* A) {
  int result = FAILURE;
  if (A != NULL) {
    if (A->rows > 0 && A->columns > 0) result = SUCCESS;
  }
  return result;
}
int eq_calc(matrix_t* A, matrix_t* B) {
  int result = FAILURE;
  if (A->rows == B->rows && A->columns == B->columns) {
    result = SUCCESS;
  }
  return result;
}
void zetloc(matrix_t* A) {
  A->columns = 0;
  A->rows = 0;
}
