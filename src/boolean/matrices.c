/**
 * matrices.c, (c) 2014, Immanuel Albrecht; Dresden University of
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

#include "matrices.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <inttypes.h>

matrix matrix_alloc(int rows, int cols) {
    matrix x = malloc(sizeof(*x));
    assert(x);

    x->rows = rows;
    x->cols = cols;
    x->width= WIDTH(cols);

    x->incidence = calloc(x->width*x->rows,sizeof(chunk_t));
    assert(x->incidence);

    return x;
}

matrix matrix_dup(matrix M) {
    matrix x = matrix_alloc(M->rows, M->cols);

    memcpy(x->incidence, M->incidence, x->rows*x->width*sizeof(*x->incidence));

    return x;
}

matrix matrix_dupT(matrix M) {
    int i;
    int j;
    matrix x = matrix_alloc(M->cols, M->rows);

    for (i = 0; i < M->rows; i++) {
        for (j = 0; j< M->cols; ++j) {
            if (INCIDESV(ROW(i,M),j)) {
                CROSSV(ROW(j,x),i);
            }
        }
    }

    return x;
}

/**
 * returns matrix product of A and transposed(B).
 */
matrix matrix_dupAxBT(matrix A, matrix B) {
    assert(A->cols == B->cols);
    matrix x = matrix_alloc(A->rows, B->rows);
    int i;
    int j;
    int w;
    chunk_t *a,*b;

    for (i = 0; i < x->rows; i++) {
        for (j = 0; j < x->cols; j++) { 
            a = ROW(i,A);
            b = ROW(j,B);
            for (w = 0; w < A->width; ++w) {
                if (*a & *b)
                {
                    CROSSV(ROW(i,x),j);
                    break;
                }
                ++a; ++b;
            }
        }            
    }

    return x;
}
/**
 * returns the transposed matrix product of A and transposed(B).
 */
matrix matrix_dupAxBT_T(matrix A, matrix B) {
    assert(A->cols == B->cols);
    matrix x = matrix_alloc(B->rows, A->rows);
    int i;
    int j;
    int w;
    chunk_t *a,*b;

    for (i = 0; i < x->cols; i++) {
        for (j = 0; j < x->rows; j++) { 
            a = ROW(i,A);
            b = ROW(j,B);
            for (w = 0; w < A->width; ++w) {
                if (*a & *b)
                {
                    CROSSV(ROW(j,x),i);
                    break;
                }
                ++a; ++b;
            }
        }            
    }

    return x;
}


/**
 * returns X(i,j) = AND_x: A(i,x)--->B(j,x)
 */
matrix matrix_dupAimpliesB(matrix A, matrix B) {
    assert(A->cols == B->cols);
    matrix x = matrix_alloc(A->rows, B->rows);
    int i;
    int j;
    int w;
    int t;
    chunk_t *a,*b;

    for (i = 0; i < x->rows; i++) {
        for (j = 0; j < x->cols; j++) { 
            a = ROW(i,A);
            b = ROW(j,B);
            t = 1;
            for (w = 0; w < A->width; ++w) {
                if ((*a & (~ *b)) != 0)
                {
                    // printf("%"PRIx64" -/-> %"PRIx64" : %"PRIx64"\n",*a,*b,*a & (~*b));
                    t = 0;
                    break;
                }
                ++a; ++b;
            }
            if (t) {
                CROSSV(ROW(i,x),j);
            }
        }            
    }

    return x;
}


void matrix_clear(matrix M) {
    assert(M);

    memset(M->incidence,0,M->width*M->rows*sizeof(chunk_t));
}

/**
 * returns -1 if matrix A is less than matrix B,
 *          1 if matrix A is bigger than matrix B
 */
int matrix_cmp(void *A, void *B) {
    int i,n;
    assert(A);
    assert(B);

    matrix M = A;
    matrix N = B;

    if (M->cols < N->cols)
        return -1;
    if (M->cols > N->cols)
        return 1;
    if (M->cols < N->cols)
        return -1;
    if (M->cols > N->cols)
        return 1;
    
    assert(M->width == N->width);

    n = M->width*M->rows;

    for (i = 0; i < n; i++) {
        if (M->incidence[i] < N->incidence[i]) {
            return -1;
        }
        else if (M->incidence[i] > N->incidence[i])
            return 1;
    }

    return 0;
}
/**
 * returns -1 if matrix row A(i,.) is less than matrix row B(j,.),
 *          1 if matrix row A(i,.) is bigger than matrix row B(j,.)
 */
int matrix_row_cmp(matrix A, int i, matrix B, int j) {
    int w;
    assert(A);
    assert(B);
    assert(A->cols == B->cols);

    chunk_t *a,*b;

    a = ROW(i,A);
    b = ROW(j,B);

    for (w = 0; w<A->width; ++w) {
        if (*a < *b) 
            return -1;
        else if (*b < *a)
            return 1;

        ++a; ++b;
    }

    return 0;
}

int matrix_isZero(matrix M) {
    int i,N;
    N = M->rows*M->width;

    for (i = 0; i < N; i++) {
        if (M->incidence[i])
            return 0;
    }

    return 1;
}

int matrix_isRowAscending(matrix M) {
    int i;

    for (i = 0; i < M->rows-1; i++) {
        if (matrix_row_cmp(M,i,M,i+1)>0)
            return 0; /* row M(i+1,.) < M(i,.) */
    }

    return 1;
}

void matrix_free(matrix x) {
    free(x->incidence);
    free(x);
}

void matrix_fput( matrix  x, FILE *f) {
    int i,j;
    for (i = 0; i < x->rows; i++) {
        for (j = 0; j < x->cols; ++j) {
            if (INCIDESV(ROW(i,x),j))
                fputs("X",f);
            else
                fputs(".",f);
        }
        fputs("\n",f);
    }
}

void row_set( chunk_t *row, const char* coded ) {
    int i;

    i=0;
    while (*coded) {
        if (*coded == 'X') {
            CROSSV(row, i);
        }
        else {
            CLEARV(row, i);
        }
        ++i; ++coded;
    }
}

void matrix_Inv(matrix M) {
    int i;
    for (i = 0; i < M->rows*M->width; i++) {
        M->incidence[i] = ~M->incidence[i];
    }

    for (i = 0; i < M->rows; i++) {
        MASKVECTOR(ROW(i,M), M->cols);
    }
}
void matrix_ones(matrix M) {
    int i;
    for (i = 0; i < M->rows*M->width; i++) {
        M->incidence[i] = ~ 0ULL;
    }

    for (i = 0; i < M->rows; i++) {
        MASKVECTOR(ROW(i,M), M->cols);
    }
}

/**
 * iterate through all matrices of the same shape as M,
 * return 0 if the resulting matrix M is the zero matrix.
 */
int matrix_increment(matrix M) {
    int i,j;
    for (i=0;i<M->rows;++i) {
        for (j=0;j<M->width-1;++j) {
            // printf("v= "PRIu64"\n",ROW(i,M)[j]);
            ROW(i,M)[j] ++;
            if (ROW(i,M)[j])
                return 1;
        }
        ROW(i,M)[M->width-1]+=BITVALUE(BITNBR(M->cols-1));
        // printf("v= %"PRIu64"\n",ROW(i,M)[j]);
        MASKVECTOR(ROW(i,M),M->cols);
        if (ROW(i,M)[M->width-1])
            return 1;
    }

    return 0;
}

int matrix_increment_row(matrix M, int i) {
    int j;

    for (j=0;j<M->width-1;++j) {
        ROW(i,M)[j] ++;
        if (ROW(i,M)[j])
            return 1;
    }

    ROW(i,M)[M->width-1]+=BITVALUE(BITNBR(M->cols-1));

    MASKVECTOR(ROW(i,M),M->cols);
    if (ROW(i,M)[M->width-1])
        return 1;

    return 0;
}

vector vector_alloc(int length) {
    vector x = malloc(sizeof(*x));
    assert(x);
    x->size = length;
    x->width = WIDTH(length);
    
    x->incidence = calloc(x->width,sizeof(*x->incidence));
    assert( x->incidence );
    
    return x;    
}

vector vector_dup(vector V) {
    vector x = vector_alloc(V->size);

    memcpy(x->incidence, V->incidence, sizeof(*x->incidence)*V->width);
    return x;
}

void vector_free(vector V) {
    free(V->incidence);
    free(V);
}

vector vector_dup_matrix_row(matrix M, int row) {
    vector x = vector_alloc(M->cols);

    assert(M->width == x->width);

    memcpy(x->incidence, ROW(row,M), sizeof(*x->incidence)*x->width);

    return x;
}

int vector_cmp(void *a, void *b) {

    assert(a);
    assert(b);

    int i;
    vector V = a;
    vector W = b;

    if (V->size < W->size)
        return -1;
    if (V->size > W->size)
        return 1;
    
    for (i = 0; i < V->width; i++) {
        if (V->incidence[i] < W->incidence[i]) {
            return -1;
        }
        else if (V->incidence[i] > W->incidence[i])
            return 1;
    }

    return 0;
}

int vector_isZero(vector V) {
    int i;

    for (i = 0; i < V->width; i++) {
        if (V->incidence[i])
            return 0;
    }

    return 1;
}
int vector_isOnes(vector V) {
    int i;

    for (i = 0; i < V->width-1; i++) {
        if (~ V->incidence[i])
            return 0;
    }

    if ((~V->incidence[V->width-1]) & CMPLASTMASK(V->size))
        return 0;

    return 1;
}

void vector_fput(vector V, FILE* f) {
    int j;
    for (j = 0; j < V->size; ++j) {
        if (INCIDESV(V->incidence,j))
            fputs("X",f);
        else
            fputs(".",f);
    }
    fputs("\n",f);
}

vector vector_dup_AminusB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = vector_dup(A);

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] &= ~ B->incidence[i];
    }

    return x;
}

vector vector_dup_AandB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = vector_dup(A);

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] &= B->incidence[i];
    }

    return x;
}
vector vector_dup_AorB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = vector_dup(A);

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] |= B->incidence[i];
    }

    return x;
}
vector vector_dup_AxorB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = vector_dup(A);

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] ^= B->incidence[i];
    }

    return x;
}
vector vector_from_cstr(const char* row) {
    vector x = vector_alloc(strlen(row));
    int i;
    for (i = 0; i < x->size; i++) {
        switch (row[i]) {
            case 'X':
            case 'x':
            case '1':
            case '+':
                CROSSV(x->incidence,i);
        }
    }

    return x;
}
vector vector_minusB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] &= ~ B->incidence[i];
    }

    return x;
}

vector vector_minusB_row(vector A, matrix B, int row) {
    assert(A->size == B->cols);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] &= ~ ROW(row,B)[i];
    }

    return x;
}
vector vector_andB_row(vector A, matrix B, int row) {
    assert(A->size == B->cols);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] &= ROW(row,B)[i];
    }

    return x;
}
vector vector_orB_row(vector A, matrix B, int row) {
    assert(A->size == B->cols);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] |= ROW(row,B)[i];
    }

    return x;
}
vector vector_xorB_row(vector A, matrix B, int row) {
    assert(A->size == B->cols);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] ^= ROW(row,B)[i];
    }

    return x;
}
vector vector_andB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] &= B->incidence[i];
    }

    return x;
}
vector vector_orB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] |= B->incidence[i];
    }

    return x;
}
vector vector_xorB(vector A, vector B) {
    assert(A->size == B->size);

    vector x = A;

    int i;
    for (i = 0; i < x->width; i++) {
        x->incidence[i] ^= B->incidence[i];
    }

    return x;
}
int vector_row_AsubsetB(matrix A, int row, vector B) {
    assert(A->cols == B->size);

    int i;
    for (i = 0; i < A->width; i++) {
        if (ROW(row,A)[i] & (~ B->incidence[i]))
            return 0;
    }

    return 1;
}


int vector_AsubsetB(vector A, vector B) {
    assert(A->size == B->size);

    int i;
    for (i = 0; i < A->width; i++) {
        if (A->incidence[i] & (~ B->incidence[i]))
            return 0;
    }

    return 1;
}
int vector_AsubsetB_row(vector A, matrix B, int row) {
    assert(A->size == B->cols);

    int i;
    for (i = 0; i < A->width; i++) {
        if (A->incidence[i] & (~ ROW(row,B)[i]))
            return 0;
    }

    return 1;
}


void vector_clear(vector V) {
    memset(V->incidence,0,V->width*sizeof(*V->incidence));

}
void vector_ones(vector V) {
    int i;
    for (i = 0; i < V->width; i++) {
        V->incidence[i] = ~ 0ULL;
    }

    MASKVECTOR(V->incidence, V->size);

}
void vector_Inv(vector V) {
    int i;
    for (i = 0; i < V->width; i++) {
        V->incidence[i] = ~ V->incidence[i];
    }

    MASKVECTOR(V->incidence, V->size);

}

int_vector int_vector_alloc(int cols) {
    int_vector x = malloc(sizeof(s_int_vector));
    assert(x);

    x->V = vector_alloc(cols);

    return x;
}

void int_vector_free(int_vector x) {
    vector_free(x->V);
    free(x);
}

positive_DNF_matrix positive_DNF_matrix_alloc(int rows, int cols) {
    positive_DNF_matrix x = malloc(sizeof(s_positive_DNF_matrix));
    assert(x);
    x->rows = rows;
    x->cols = cols;

    x->entries = cp_vector_create(1);

    return x;
}

void positive_DNF_matrix_free(positive_DNF_matrix x) {
    cp_vector_destroy_custom(x->entries, int_vector_free);
    free(x);
}

/**
 *
 * returns the matrix C = (DNF)x(B^T)^T
 * (if DNF is a DSF matrix determining which compound-skills can
 *  be created by combining some elementary skills, and
 *  B is a matrix that contains the available elementary skills for some subjects,
 *  then the resulting matrix is the matrix of compound-skills available to
 *  the subjects)
 */
matrix matrix_dup_DNFxBT_T(positive_DNF_matrix DNF, matrix B) {
    assert(DNF->cols == B->cols);
    matrix C = matrix_alloc(B->rows, DNF->rows);

    int i,N,j;
    N = cp_vector_size(DNF->entries);

    int_vector iV;
    vector Q;

    for (j = 0; j < B->rows; j++) {
        Q = vector_dup_matrix_row(B,j);
        for (i = 0; i < N; i++) {
            iV = cp_vector_element_at(DNF->entries, i);

            if (vector_AsubsetB(iV->V,Q)) {
                CROSSV(ROW(j,C),iV->i);
            }
        }
        vector_free(Q);
    }

    return C;
}

void DNF_fput(const char* name,positive_DNF_matrix DNF, FILE* f) {
    int i,N;
    int empty;
    int j,M;
    int_vector iV;
    int l = strlen(name);
    int k;
 
    N = DNF->rows;
    M = cp_vector_size(DNF->entries);

    for (i = 0; i < N; i++) {
        empty = 1;
        for (j = 0; j < M; j++) {
            iV = cp_vector_element_at(DNF->entries, j);

            if (iV->i == i) {
                if (empty)
                    fprintf(f,"%s(%2d,.) = ",name,i);
                else {
                    for (k=0;k<l;++k)
                        fputs(" ",f);
                    fputs("      \\/ ",f);
                }

                vector_fput(iV->V,f);
                empty = 0;
            }
        }       
        if (empty)
           fprintf(f,"%s(%2d,.) = _|_\n",name,i);
    }
}

int_vector int_vector_alloc2(int i, vector V) {

    int_vector x = malloc(sizeof(s_int_vector));
    assert(x);

    x->i = i;
    x->V = V;

    return x;
}

void DNF_add_row_from_cstr(positive_DNF_matrix DNF, int row, const char* data) {
    assert(row >= 0);
    assert(row < DNF->rows);
    vector V = vector_from_cstr(data);
    assert(V->size == DNF->cols);
    
    cp_vector_add_element(DNF->entries, int_vector_alloc2(row, V));
}

positive_DNF_matrix DNF_dup_matrix(matrix M) {
    positive_DNF_matrix DNF = positive_DNF_matrix_alloc(M->rows, M->cols);
    
    int i,j;

    for (i = 0; i < M->rows; i++) {
        for (j = 0; j < M->cols; j++) {
            if (INCIDESV(ROW(i,M),j)) {
                DNF_add_row_atom(DNF, i, j);
            }            
        }
    }

    return DNF;
}

void DNF_add_row_atom(positive_DNF_matrix DNF, int row, int col) {
    assert(row >= 0);
    assert(row < DNF->rows);
    assert(col >= 0);
    assert(col < DNF->cols);

    vector V = vector_alloc(DNF->cols);
    CROSSV(V->incidence,col);
    cp_vector_add_element(DNF->entries, int_vector_alloc2(row,V));
}

cross_stats cross_stats_alloc(int rows, int cols) {
    cross_stats x = malloc(sizeof(s_cross_stats));
    assert(x);

    x->rows = rows;
    x->cols = cols;

    x->rows_X = calloc(rows,sizeof(int));
    x->rows_n = calloc(rows,sizeof(int));
    x->cols_X = calloc(cols,sizeof(int));
    x->cols_n = calloc(cols,sizeof(int));

    return x;
}

void cross_stats_free(cross_stats x) {
    free(x->rows_X);
    free(x->rows_n);
    free(x->cols_X);
    free(x->cols_n);
    free(x);
}

cross_stats cross_stats_alloc_DNFxBT_T(positive_DNF_matrix DNF, matrix B) {
    assert(DNF->cols == B->cols);
    matrix C = matrix_alloc(B->rows, DNF->rows);
    cross_stats S = cross_stats_alloc(B->rows, DNF->rows);


    int i,N,j;
    N = cp_vector_size(DNF->entries);

    int_vector iV;
    vector Q;

    for (j = 0; j < B->rows; j++) {
        Q = vector_dup_matrix_row(B,j);
        for (i = 0; i < N; i++) {
            iV = cp_vector_element_at(DNF->entries, i);

            if (vector_AsubsetB(iV->V,Q)) {
                if (! INCIDESV(ROW(j,C),iV->i)) {
                    CROSSV(ROW(j,C),iV->i);
                    S->rows_X[j] ++;
                    S->cols_X[iV->i] ++;
                } else {
                    S->rows_n[j] ++;
                    S->cols_n[iV->i] ++;
                }
            }
        }
        vector_free(Q);
    }

    matrix_free(C);

    return S;
}

int vector_count_crosses(vector V) {
    int count = 0;
    int i;
    for (i = 0; i < V->size; i++) {
        if (INCIDESV(V->incidence,i))
            count++;
    }

    return count;
}
/**
 * calculate the implications up to n_grams and return the
 * nontrivial ones as avltree,
 * keys are premise-vectors, values are conclusion-vectors.
 *
 * delete with
 *
 * cp_avltree_destroy_custom( ..., vector_free, vector_free )
 */
cp_avltree * matrix_column_implications(matrix M, int n_grams) {
    cp_avltree *implications = cp_avltree_create(vector_cmp);

    int* choices = calloc(n_grams+1,sizeof(int));
    int* indices = calloc(M->cols,sizeof(int));
    int i;
    int j,J,crosses;
    int l,k;
    int sum;
    vector row,r2;
    vector prem;
    vector concl;
    prem = vector_alloc(M->cols);
    for (i = 0; i < M->rows; i++) {
        row = vector_dup_matrix_row(M,i);

        j=0;
        for (l=0; l<M->cols;++l) {
            if (INCIDESV(row->incidence,l)) {
                indices[j] = l;
                ++j;
            }
        }

        J = vector_count_crosses(row);
        crosses = J;
        if (J>n_grams)
            J = n_grams;
        for (j = 1; j < J; j++) {
            for (l=0;l<j+1;++l) {
                choices[l] = 0;
            }

            do {
                vector_clear(prem);

                for (l=0;l<j;++l) {
                    if (INCIDESV(prem->incidence, indices[choices[l]]))
                    {
                        break;
                    } else CROSSV(prem->incidence, indices[choices[l]]);
                }

                if (vector_count_crosses(prem) == j) {
                    if (! cp_avltree_contains(implications, prem)) {

                        concl = vector_dup_matrix_row(M,i);
                        for (k = 0; k < M->rows; k++) {
                            if (k==i)
                                continue;

                            r2 = vector_dup_matrix_row(M,k);

                            if (vector_AsubsetB(prem,r2)) {
                                vector_andB(concl,r2);
                            }

                            vector_free(r2);
                        }

                        /*puts("IMPL ");
                        vector_fput(prem,stdout);
                        puts("---> ");
                        vector_fput(concl,stdout); */
                        cp_avltree_insert(implications, vector_dup(prem), concl);

                    }
                }

                l = 0;
                choices[l]++;
                while (choices[l] == crosses) {
                    choices[l] = 0;
                    l++;
                    choices[l] ++;
                }

                sum = 0;
                for (l=0;l<j;++l) {
                    sum += choices[l];
                }
            } while (sum);
        }
        vector_free(row);        
    }

    vector_free(prem);
    free(choices);
    free(indices);

    return implications;
}

void DNF_clear(positive_DNF_matrix DNF) {
    cp_vector_destroy_custom(DNF->entries, int_vector_free);
    DNF->entries = cp_vector_create(1);
}
/** 
 *
 * remove those rows who are supersets of other rows, i.e.
 * if
 *    (1, "XX..")
 *  and
 *    (1, "X...")
 *  are part of the DNF, then remove the first one entirely.
 */
void DNF_clarify_join(positive_DNF_matrix DNF) {
    int i,N,j;
    N = cp_vector_size(DNF->entries);
    int_vector A,B;
    vector obsolete = vector_alloc(DNF->entries);

    for (i = 0; i < N; i++) {
        if (INCIDESV(obsolete->incidence,i))
            continue;
        A = cp_vector_element_at(DNF->entries,i);

        for (j = 0; j<N; j++) {
            if (INCIDESV(obsolete->incidence,j))
                continue;
            if (i == j) 
                continue;
            B = cp_vector_element_at(DNF->entries,j);
            if (A->i != B->i)
                continue;
            if (vector_cmp(A->V,B->V) == 0) {
                if (i < j)
                {
                    CROSSV(obsolete->incidence,j);
                }
            } else if (vector_AsubsetB(A->V, B->V)) {
                    CROSSV(obsolete->incidence,j);
            }
        }
    }
    for (i = N-1; i >= 0; i--) {
        if (INCIDESV(obsolete->incidence, i)) {
            int_vector_free(cp_vector_element_at(DNF->entries,i));
            cp_vector_remove_element_at(DNF->entries,i);
        }
    }

    vector_free(obsolete);
}
/** returns number of bits that are set in V */
int vector_checksum(vector V) {
    int count;
    int i;

    count = 0;
    for (i = 0; i < V->size; i++) {
        if (INCIDESV(V->incidence, i))
            ++count;
    }

    return count;
}
