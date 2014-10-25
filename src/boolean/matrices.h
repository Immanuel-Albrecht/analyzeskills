/*
 * matrices.h, (c) 2014, Immanuel Albrecht; Dresden University of
 * Technology, Professur für die Psychologie des Lernen und Lehrens
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MATRICES_R20UTAJQ

#define MATRICES_R20UTAJQ

#include <stdio.h>
#include <inttypes.h>
#include <cprops/vector.h>
#include <cprops/avl.h>

typedef uint64_t chunk_t;

/**
 * @def OFFSET(x)
 *
 * get the offset of the x-th bit in an 64-bit integer vector
 */

#define OFFSET(x) ((unsigned)(x)>>6)
/**
 * @def BITNBR(x)
 *
 * get the remainder of the x-th bit in an 64-bit integer vector, i.e. 65=64+ @a 1
 */
#define BITNBR(x) (((unsigned)(x))&(63))

/**
 * @def WIDTH(x)
 *
 * determine the length of an 64-bit integer vector that can hold x bits
 */
#define WIDTH(x) ((((unsigned)(x))&(63))?((unsigned)(x)/64)+1:((unsigned)(x)/64))

/**
 * @def BITVALUE(x)
 *
 * gives the bit-value of the x-th bit. (Note that bit 0 is the most,
 * and bit 63 is the least significant bit)
 */

#define BITVALUE(x) ((1ULL<<(63-BITNBR(x))))

/**
 * @def CRIMPVALUE(x)
 *
 * gives an 64-bit integer that has set the bits 0 through x.
 *
 * CRIMPVALUE(0) == 0x8000000000000000
 * CRIMPVALUE(1) == 0xc000000000000000
 * CRIMPVALUE(2) == 0xe000000000000000
 * CRIMPVALUE(3) == 0xf000000000000000
 * CRIMPVALUE(4) == 0xf800000000000000
 * etc.
 *
 */

#define CRIMPVALUE(x) ((~(0ULL))>>(63-(BITNBR(x)))<<(63-BITNBR(x)))

/**
 * @def MASKVECTOR(v,x)
 *
 * set the unused attribute bits to zero. (i.e. attributes == 100 -> width == 2, BITNBR(99) == 35)
 * where v is a 64-bit integer vector, and x is the number used bits.
 */
#define MASKVECTOR(v,x) {if (BITNBR((x))) { *((v)+OFFSET((x)-1)) = ( (*((v)+OFFSET((x)-1))>>(63-BITNBR((x)-1))) ) << (63-BITNBR((x)-1));  }}


/**
 * @def CMPLASTMASK(x)
 *
 * x __ size of vector in bits
 *
 * gets the comparison bit-mask for the last element of a bit-vector, i.e.
 *   if BITNBR(x) == 0,  then CMPLASTMASK(x) == ~(0ULL),
 *   otherwise                CMPLASTMASK(x) == CRIMPVALUE(BITNBR(x-1))
 */

#define CMPLASTMASK(x) ( (BITNBR(x) ? CRIMPVALUE(BITNBR(x-1)) : ~(0ULL)) )

/**
 * @def  CROSSV(v,x)
 * crosses the x-th attribute of an attribute vector
 */
#define CROSSV(v,x) { *((v)+OFFSET(x)) |= BITVALUE(x); }

/**
 * @def  CLEARV(v,x)
 * clears the x-th attribute of an attribute vector
 */
#define CLEARV(v,x) { *((v)+OFFSET(x)) &= ~ (BITVALUE(x)); }

/**
 * @def INCIDESV(v,x)
 * checks whether the x-th attribute of an attribute vector is crossed
 */
#define INCIDESV(v,x) (  ( *((v)+OFFSET(x)) >> (63-BITNBR(x)) ) & 1  )

/**
 * @def ROW(g,I)
 * gives the attribute vector for a given object
 */
#define ROW(g,I) ((I)->incidence + ((I)->width * (g)))

typedef struct t_matrix {
    int rows,cols;
    int width; // = WIDTH(cols)
    chunk_t *incidence;
} s_matrix;
typedef s_matrix *matrix;


matrix matrix_alloc(int rows, int cols);
matrix matrix_dup(matrix M);
matrix matrix_dupT(matrix M);
matrix matrix_dupAxBT(matrix A, matrix B);
matrix matrix_dupAxBT_T(matrix A, matrix B);
matrix matrix_dupAimpliesB(matrix A, matrix B);
void matrix_clear(matrix M);
void matrix_ones(matrix M);
int matrix_cmp(void *A, void *B);
int matrix_row_cmp(matrix A, int i, matrix B, int j);
int matrix_isZero(matrix M);
int matrix_isRowAscending(matrix M);
void matrix_free(matrix x);
void matrix_fput( matrix  x, FILE *f);
void row_set( chunk_t *row, const char* coded );
void matrix_Inv(matrix M);
int matrix_increment(matrix M);
int matrix_increment_row(matrix M, int row);

cp_avltree * matrix_column_implications(matrix M, int n_grams);


typedef struct t_vector {
    int size;
    int width; // = WIDTH(size)
    chunk_t *incidence;
} s_vector;

typedef s_vector *vector;

vector vector_alloc(int length);
vector vector_dup(vector V);
vector vector_dup_matrix_row(matrix M, int row);
vector vector_from_cstr(const char* row);
int vector_isZero(vector V);
int vector_isOnes(vector V);
void vector_free(vector V);

void vector_clear(vector V);
void vector_ones(vector V);
void vector_Inv(vector V);

int vector_checksum(vector V);

vector vector_dup_AminusB(vector A, vector B);
vector vector_dup_AandB(vector A, vector B);
vector vector_dup_AorB(vector A, vector B);
vector vector_dup_AxorB(vector A, vector B);
void vector_fput( vector V, FILE *f);

vector vector_minusB(vector A, vector B);
vector vector_andB(vector A, vector B);
vector vector_orB(vector A, vector B);
vector vector_xorB(vector A, vector B);

vector vector_minusB_row(vector A, matrix B, int row);
vector vector_andB_row(vector A, matrix B, int row);
vector vector_orB_row(vector A, matrix B, int row);
vector vector_xorB_row(vector A, matrix B, int row);

int vector_cmp(void *V, void *W);

int vector_AsubsetB(vector A, vector B);
int vector_row_AsubsetB(matrix A, int row, vector B);
int vector_AsubsetB_row(vector A, matrix B, int row);

int vector_AdisjointB(vector A, vector B);

#define vector_count_crosses vector_checksum

/** @brief (int, vector)-pair
 *
 * 
 * This structure is used as a conjunctive-element in an DNF pseudo-matrix,
 * where i is the column that gets a mark if the input vector is a superset of V
 *
 */
typedef struct t_int_vector {
    int i;//!< output column that shall receive the cross.
    vector V;//!< produce a cross, if the input vector is a superset of this.
} s_int_vector;
typedef s_int_vector *int_vector;

int_vector int_vector_alloc(int cols);
void int_vector_free(int_vector x);
int_vector int_vector_alloc2(int i, vector V);

typedef struct t_positive_DNF_matrix {
    int rows,cols;
    cp_vector *entries; // array of int_vector
} s_positive_DNF_matrix;
typedef s_positive_DNF_matrix *positive_DNF_matrix;

positive_DNF_matrix positive_DNF_matrix_alloc(int rows, int cols);
void positive_DNF_matrix_free(positive_DNF_matrix x);
#define DNF_alloc positive_DNF_matrix_alloc
#define DNF_free positive_DNF_matrix_free

matrix matrix_dup_DNFxBT_T(positive_DNF_matrix DNF, matrix B);

void DNF_fput(const char* name,positive_DNF_matrix DNF, FILE* f);
void DNF_add_row_from_cstr(positive_DNF_matrix DNF, int row, const char* data);
void DNF_clear(positive_DNF_matrix DNF);

positive_DNF_matrix DNF_dup_matrix(matrix M);
void DNF_add_row_atom(positive_DNF_matrix DNF, int row, int col);

void DNF_clarify_join(positive_DNF_matrix DNF);

typedef struct t_cross_stats {
    int rows, cols;
    int *rows_X, *cols_X; //number of regular crosses per row/col
    int *rows_n, *cols_n; //number of overproduced crosses per row/col
} s_cross_stats;
typedef s_cross_stats *cross_stats;

cross_stats cross_stats_alloc();
void cross_stats_free(cross_stats x);
cross_stats cross_stats_alloc_DNFxBT_T(positive_DNF_matrix DNF, matrix B);

#endif /* end of include guard: MATRICES_R20UTAJQ */
