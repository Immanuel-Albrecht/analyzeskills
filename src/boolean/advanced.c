/*
 * advanced.c, (c) 2014, Immanuel Albrecht; Dresden University of
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

#include "advanced.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include <cprops/heap.h>
#include <cprops/util.h>
#include <cprops/hashtable.h>
#include <cprops/vector.h>
#include <cprops/multimap.h>
#include <cprops/trie.h>
#include <cprops/avl.h>
#include <cprops/rb.h>


/** @brief tries to find matrix A, such that A*transposed(B) = C
 *
 * @param A matrix, going to be changed.
 *          In case of a return value >0, A contains the desired matrix.
 *
 * @param B skill matrix, subjects=rows.
 * @param C response matrix, subjects=rows
 *
 *
 * Tries to find a valid matrix B, s.t. A{gprod}(B^T) = C^T.
 *
 * @returns 1, if the problem could be solved, 0 otherwise.
 */
int matrix_solveB_AimpliesBeqC(matrix A, matrix B, matrix C) {
    int i,j,k;
    assert(A->cols == B->cols);
    assert(A->rows == C->rows);
    assert(B->rows == C->cols);

    chunk_t *a,*b;
    matrix AimpliesB;

    matrix_clear(B);
    /*  puts("A=");
        matrix_fput(A,stdout);
        puts("C=");
        matrix_fput(C,stdout); */

    for (i = 0; i < A->rows; i++) {
        a = ROW(i,A);
        for (j = 0; j < B->rows; ++j) {
            if (INCIDESV(ROW(i,C),j)) {
                // printf("%d solves %d\n",j,i);
                b = ROW(j,B);
                for (k = 0; k < B->cols; ++k) {
                    if (INCIDESV(a,k))
                    {
                        // printf("%d has %d\n",j,k);
                        CROSSV(b,k);
                    }
                }
            }                    
        }
    }

    /*  puts("B=");
        matrix_fput(B,stdout); */

    AimpliesB = matrix_dupAimpliesB(A,B);
    if (matrix_cmp(AimpliesB,C)) {
        /*  puts("A-->B=");
            matrix_fput(AimpliesB,stdout);
            printf("cmp=%d\n",matrix_cmp(C,AimpliesB)); */
        matrix_free(AimpliesB);
        return 0;
    }

    matrix_free(AimpliesB);
    return 1;
}
/**
 * @struct t_fjb_data
 *
 * @brief struct private to fjb_callback .
 **/
typedef struct t_fjb_data {
    cp_vector *base;

} s_fjb_data;
typedef s_fjb_data *fjb_data;

static int fjb_callback(void* entry, void* user) {
    cp_avlnode *node=entry;
    vector V= node->key;
    fjb_data data = user;


    if (vector_isZero(V))
        return 0;

    /* check whether V is the join of the relevant base vectors */

    vector J = vector_alloc(V->size);
    vector diff = vector_dup(V);
    vector x,y;

    int i,N;
    int j;

    N = cp_vector_size(data->base);

    for (i = 0; i < N; i++) {
        x = cp_vector_element_at(data->base,i);
        if (vector_AsubsetB(x,V)) {
            vector_orB(J,x);
            vector_minusB(diff,x);
        }
    }

    if (vector_cmp(V,J)) {
        /*
         * V is not covered by the base! 
         *
         * if (diff) was part of the base, we would have a 
         * representation of V.
         *
         **/
        
        /* remove the parts that are covered by diff from the
         * other base vectors that*/
        for (i = 0; i < N; i++) {
            x = cp_vector_element_at(data->base,i);
            if (vector_AsubsetB(diff,x))
                vector_minusB(x,diff);  
        }

        /* this might have created vectors that can be represented
         * by other vectors
         */
        for (i = 0; i < N; i++) {
            x = cp_vector_element_at(data->base,i);
            vector_clear(J);
            for (j = 0; j < N; j++) {
                if (i==j)
                    continue;

                y = cp_vector_element_at(data->base,j);
                if (vector_AsubsetB(y,x)) {
                    vector_orB(J,y);
                }
            }

            if (vector_cmp(J,x)==0) {
                /* x is in the join-hull of the other base vectors */
                cp_vector_remove_element_at(data->base, i);
                N --;
                i --;
                vector_free(x);
            }
        }

        cp_vector_add_element(data->base, diff);

    } else
        vector_free(diff);

    vector_free(J);

    return 0;
}

/** @brief find an orded set of bit-vectors, that acts as base on vector_set
 *
 * takes cp_avltree* whose keys are of the type vector,
 * and calculates a base, such that each key is a the
 * join of a unique combination of these base elements.
 *
 * @returns a cp_vector* containing vectors,
 *       has to be freed by the caller using
 *        cp_vector_destroy_custom( find_join_base(...) ,vector_free);
 */
cp_vector *find_join_base(cp_avltree *vector_set) {
    s_fjb_data s_data;
    fjb_data data = &s_data;
    data->base = cp_vector_create(1);

    cp_avltree_callback(vector_set, fjb_callback, data);

    return data->base;
}

/** @brief get "solves"-matrix A and skill matrix B from response matrix C
 *
 * Input: C(.,.), where
 *   each row C(i,.) corresponds to the response pattern of
 *   one subject
 *
 * will generate appropriate matrices A and B, such that
 *
 *   A * B^T = C^T
 *
 * holds. The dimensions of A and B are determined and
 * then the matrices are created accordingly.
 *
 * The caller has to use matrix_free on both *A and *B
 * to free the memory allocated by this routine.
 *
 * In this case, the rows A(i,.) contain information
 * which skills will lead to a solution, whereas 
 * the rows B(i,.) contain information about the available
 * skills.
 */
void generate_skill_matrix(matrix *pA, matrix *pB, matrix C) {
    cp_avltree *responses = cp_avltree_create(vector_cmp);
    cp_vector *base;

    vector V;
    int i,N,j;


    N = C->rows;

    for (i = 0; i < N; i++) {
        V = vector_dup_matrix_row(C,i);
        if (cp_avltree_contains(responses,V)) {
            vector_free(V);
        } else {
            cp_avltree_insert(responses,V,1);
        }
    }


    base = find_join_base(responses);
    N = cp_vector_size(base);

    /*
    printf("Response Base: %d\n",N);
    for (i = 0; i < N; i++) {
        printf("x%2d = ",i);
        vector_fput(cp_vector_element_at(base,i),stdout);
    }
    */

    *pA = matrix_alloc(C->cols, N);
    *pB = matrix_alloc(C->rows, N);

    matrix A = *pA;
    matrix B = *pB;

    /* the response base gives
     * the information, which items 
     * can be solved by a corresponding
     * compound-skill
     */

    for (i = 0; i < A->rows; i++) {
        for (j = 0; j < A->cols; j++) {
            V = cp_vector_element_at(base,j);
            if (INCIDESV(V->incidence,i))
                CROSSV(ROW(i,A),j);
        }
    }

    /*
     * every candidate gets assigned all those
     * skills, that are coherent with his/her
     * response. I.e. if the candidate fails
     * an item that is solved by that skill,
     * s/he will not receive a cross for it.
     */

    for (i = 0; i < B->rows; i++) {
        for (j = 0; j < B->cols; j++) {
            V = vector_dup_matrix_row(C,i);

            if (vector_AsubsetB(cp_vector_element_at(base,j),V))
                CROSSV(ROW(i,B),j);

            vector_free(V);
        }
    }


    cp_vector_destroy_custom(base,vector_free);
    cp_avltree_destroy_custom(responses, vector_free,0);
}

typedef struct t_grmc_cb {
    cp_avltree *recent;
    cp_avltree *set;
    matrix C;

} s_grmc_cb;
typedef s_grmc_cb *grmc_cb;

static int grmc_callback(void* entry, void* user) {
    cp_avlnode *node = entry;
    grmc_cb data = user;
    vector V = node->key;
    vector W;
    int i;

    for (i = 0; i < data->C->rows; i++) {
        W = vector_dup_matrix_row(data->C,i);

        vector_andB(W,V);

        if (! cp_avltree_contains(data->set,W)) {
            cp_avltree_insert(data->set, W, 1);
            cp_avltree_insert(data->recent, W, 1);
            //puts("W=");
            //vector_fput(W,stdout);
        } else
            vector_free(W);
    }

    return 0;
}

/** @brief return the set of all possible, non-trivial meets of rows of C. 
 * I.e. (excluding the empty meet)
 *
 *
 * @returns a new avltree that contains all meets
 * of rows of the matrix C.
 * The result has to be destroyed by the caller, using
 *   cp_avltree_destroy_custom( .., vector_free, 0);
 *
 */
cp_avltree* get_row_meet_closure(matrix C) {
    cp_avltree *set = cp_avltree_create(vector_cmp);
    cp_avltree *recent = cp_avltree_create(vector_cmp);
    cp_avltree *last_recent=0;
    s_grmc_cb sdata;
    grmc_cb data=&sdata;
    vector V;
    int i;

    data->C = C;
    data->set = set;

    for (i = 0; i < C->rows; i++) {
        V = vector_dup_matrix_row(C,i);

        cp_avltree_insert(set,V,1);
        cp_avltree_insert(recent,V,1);
    }

    while (cp_avltree_count(recent)) {
        if (last_recent)
            cp_avltree_destroy(last_recent);

        last_recent = recent;
        recent = cp_avltree_create(vector_cmp);

        data->recent = recent;

        cp_avltree_callback(last_recent, grmc_callback, data);
        
    }

    if (last_recent)
        cp_avltree_destroy(last_recent);

    cp_avltree_destroy(recent);

    return set;
}

int solve_DNFxBT_T_eqC_row(positive_DNF_matrix DNF, matrix B, matrix C, int row) {
    vector Bi = vector_alloc(B->cols);
    vector newBi;
    vector Ci = vector_dup_matrix_row(C,row);
    int i,j,M,k;
    int found, still_good;
    int_vector DNFj;
    M = cp_vector_size(DNF->entries);
    cp_vector *stack = cp_vector_create(1);
    cp_vector *stacki = cp_vector_create(1);
    cp_vector *stackj = cp_vector_create(1);
    int count;

    for (i = 0; i < Ci->size; i++) {
        if (INCIDESV(Ci->incidence, i)) {
            
            /** check whether we need to adjust Bi */

            found = 0;
            for (j = 0; j < M; j++) {
                DNFj = cp_vector_element_at(DNF->entries, j);
                if (DNFj->i == i) {
                    if (vector_AsubsetB(DNFj->V,Bi))
                    {
                        found = 1; break;
                    }
                }
            }
            if (!found) {
                /** we have to add some skills to Bi */
                for (j = 0; j < M; j++) {
continue_search:
                    DNFj = cp_vector_element_at(DNF->entries, j);
                    if (DNFj->i == i) {
                        /** we found a candidate */
                        newBi = vector_dup(Bi);
                        vector_orB(newBi, DNFj->V);

                        /** check whether it still works **/
                        still_good = 1;
                        for (k = 0; k < M; k++) {
                            DNFj = cp_vector_element_at(DNF->entries, k);
                            if (INCIDESV(Ci->incidence, DNFj->i))
                                continue;

                            if (vector_AsubsetB(DNFj->V,newBi)) {
                                still_good = 0;
                                break;
                            }
                        }

                        if (still_good) {
                            cp_vector_add_element(stack,Bi);
                            cp_vector_add_element(stacki,i);
                            cp_vector_add_element(stackj,j);
                            Bi = newBi;
                            found = 1;
                            break;
                        } else {
                            vector_free(newBi);
                        }
                    }
                }
            }
            if (!found) {
                /** still no candidate --> unroll the stack */
                count = cp_vector_size(stack);
                if (count == 0) {
                    /** nope, we are through. */

                    cp_vector_destroy(stacki);
                    cp_vector_destroy(stackj);
                    cp_vector_destroy_custom(stack, vector_free);
                    vector_free(Ci);
                    vector_free(Bi);

                    return 0;
                } else {
                    vector_free(Bi);
                    
                    Bi = cp_vector_element_at(stack, count-1);
                    i = cp_vector_element_at(stacki, count-1);
                    j = cp_vector_element_at(stackj, count-1);

                    cp_vector_remove_element_at(stack, count-1);
                    cp_vector_remove_element_at(stacki, count-1);
                    cp_vector_remove_element_at(stackj, count-1);

                    /** reenter the next try of the previous recursion */
                    found = 0;
                    ++j;
                    goto continue_search;

                }
            }
        }
    }

    /** copy the solution */
    memcpy(ROW(row,B), Bi->incidence, sizeof(chunk_t)*Bi->width);

    cp_vector_destroy(stacki);
    cp_vector_destroy(stackj);
    cp_vector_destroy_custom(stack, vector_free);
    vector_free(Ci);
    vector_free(Bi);
    return 1;
}

/**
 * tries to fill in values into matrix B, 
 * such that the equation
 *   DNFxBTT = C
 * is satisfied.
 *
 * returns 0 on failure, otherwise B is filled with such a matrix
 */
int solve_DNFxBT_T_eqC(positive_DNF_matrix DNF, matrix B, matrix C) {
    assert(DNF);
    assert(B);
    assert(C);
    assert(DNF->cols == B->cols);
    assert(DNF->rows == C->cols);
    assert(B->rows == C->rows);

    int i;

    for (i = 0; i < B->rows; i++) {
        if (solve_DNFxBT_T_eqC_row(DNF,B,C,i)==0) {
            printf("Couldn't solve row %d\n",i);
            return 0;
        }
    }

    return 1;
}

static int solve_DNF_from_BandC_step(positive_DNF_matrix DNF, matrix B, matrix C, int row) {
    int i,j,M;

    vector Ci = vector_dup_matrix_row(C,row);
    vector Bi = vector_dup_matrix_row(B,row);
    int_vector Ai;
    vector Xi,Yi;
    int has_it;

    M = cp_vector_size(DNF->entries);

    for (i = 0; i < Ci->size; i++) {
        if (INCIDESV(Ci->incidence, i)) {
            has_it = 0;
            for (j = 0; j < M; j++) {
                Ai = cp_vector_element_at(DNF->entries,j);
                if (Ai->i == i) {
                    if (vector_AsubsetB(Ai->V,Bi))
                    {
                        has_it = 1; break;
                    }
                }
            }
            if (!has_it) {
                for (j=0; j<B->rows; ++j) {
                    if (j==row)
                        continue;

                    Xi = vector_dup_matrix_row(C,j);
                    if (!INCIDESV(Xi->incidence, i))
                    {
                        Yi = vector_dup_matrix_row(B,j);
                        if (vector_AsubsetB(Bi,Yi))
                        {
                            /** we're doomed ... */
                            vector_free(Xi);
                            vector_free(Yi);
                            vector_free(Bi);
                            vector_free(Ci);
                            /** cannot satisfy DNF*B=C */
                            return 0;
                        }
                        vector_free(Yi);
                    }
                    vector_free(Xi);
                    
                }

                Ai = int_vector_alloc2(i,vector_dup(Bi));
                cp_vector_add_element(DNF->entries,Ai);
            }
        }       
    }

    vector_free(Ci);
    vector_free(Bi);

    return 1;
}

/**
 * this routines searches for a DNF, such that
 *  DNFxBTT = C
 * has a solution.
 * Input are both B and C
 *
 * returns 0 on failure, otherwise DNF contains the desired solution.
 */
int solve_DNF_from_BandC(positive_DNF_matrix DNF, matrix B, matrix C) {
    assert(DNF);
    assert(B);
    assert(C);
    assert(DNF->cols == B->cols);
    assert(DNF->rows == C->cols);
    assert(B->rows == C->rows);

    DNF_clear(DNF);

    int i;
    for (i = 0; i < B->rows; i++) {
        if (! solve_DNF_from_BandC_step(DNF,B,C,i)) {
            return 0;
        }
    }

    DNF_clarify_join(DNF);

    return 1;
}

static int partially_solve_DNF_from_BandC_step(positive_DNF_matrix DNF, matrix B, matrix C, int row, int rows) {
    int i,j,M;

    vector Ci = vector_dup_matrix_row(C,row);
    vector Bi = vector_dup_matrix_row(B,row);
    int_vector Ai;
    vector Xi,Yi;
    int has_it;

    M = cp_vector_size(DNF->entries);

    for (i = 0; i < Ci->size; i++) {
        if (INCIDESV(Ci->incidence, i)) {
            has_it = 0;
            for (j = 0; j < M; j++) {
                Ai = cp_vector_element_at(DNF->entries,j);
                if (Ai->i == i) {
                    if (vector_AsubsetB(Ai->V,Bi))
                    {
                        has_it = 1; break;
                    }
                }
            }
            if (!has_it) {
                for (j=0; j< rows; ++j) {
                    if (j==row)
                        continue;

                    Xi = vector_dup_matrix_row(C,j);
                    if (!INCIDESV(Xi->incidence, i))
                    {
                        Yi = vector_dup_matrix_row(B,j);
                        if (vector_AsubsetB(Bi,Yi))
                        {
                            /** we're doomed ... */
                            vector_free(Xi);
                            vector_free(Yi);
                            vector_free(Bi);
                            vector_free(Ci);
                            /** cannot satisfy DNF*B=C */
                            return 0;
                        }
                        vector_free(Yi);
                    }
                    vector_free(Xi);
                    
                }

                Ai = int_vector_alloc2(i,vector_dup(Bi));
                cp_vector_add_element(DNF->entries,Ai);
            }
        }       
    }

    vector_free(Ci);
    vector_free(Bi);

    return 1;
}


/**
 * same as   solve_DNF_from_BandC,
 * but only considering the first rows of a matrix
 *
 */

int partially_solve_DNF_from_BandC(positive_DNF_matrix DNF, matrix B, matrix C, int rows) {
    assert(DNF);
    assert(B);
    assert(C);
    assert(DNF->cols == B->cols);
    assert(DNF->rows == C->cols);
    assert(B->rows == C->rows);

    DNF_clear(DNF);

    int i;
    for (i = 0; i < rows; i++) {
        if (! partially_solve_DNF_from_BandC_step(DNF,B,C,i,rows)) {
            return 0;
        }
    }

    return 1;
}

/**
 *
 * takes a matrix C and the schematics of DNF and B,
 * and tries to find a skill assignment as well as
 * a DNF matrix, such that those two components
 * produce the response given by C.
 *
 * 
 */
int find_DNF_B_from_DNFxBTT(positive_DNF_matrix DNF, matrix B, matrix C) {
    assert(DNF);
    assert(B);
    assert(C);
    assert(DNF->cols == B->cols);
    assert(DNF->rows == C->cols);
    assert(B->rows == C->rows);

    int i;
    int go_deeper;

    matrix_clear(B);

    for (i = 0; i < B->rows; i++) {
try_again:
        go_deeper = 0;
        do {
            if (partially_solve_DNF_from_BandC(DNF,B,C,i+1))
            {
                go_deeper = 1;
                break;
            }
        } while (matrix_increment_row(B,i));
        if (go_deeper == 0) {
            if (i == 0) {
                return 0;
            } else {
                i --;
                matrix_increment_row(B,i);
                goto try_again;
            }
        }
    }

    DNF_clarify_join(DNF);

    return 1;
}

/**
 * counts the number of possible skill matrices that allow solving C
 */

int count_Bs_from_DNFxBTT(positive_DNF_matrix DNF, matrix B, matrix C) {
    assert(DNF);
    assert(B);
    assert(C);
    assert(DNF->cols == B->cols);
    assert(DNF->rows == C->cols);
    assert(B->rows == C->rows);

    int i;
    int go_deeper;
    int count;

   // matrix check;

    matrix_clear(B);
    count = 0;

    for (i = 0; i < B->rows; i++) {
try_again:
        go_deeper = 0;
        do {
            if (partially_solve_DNF_from_BandC(DNF,B,C,i+1))
            {
                go_deeper = 1;
                break;
            }
        } while (matrix_increment_row(B,i));
        if ((go_deeper)&&(i== B->rows-1)) {

            /*check = matrix_dup_DNFxBT_T(DNF,B);
            assert(matrix_cmp(check,C)==0);
            matrix_free(check);*/

            ++count;
            if (count % 15000 == 0) {
                printf(".");
                fflush(stdout);
            }
            go_deeper = 0;
        }
        if (go_deeper == 0) {
            if (i == 0) {
                return 0;
            } else {
                i --;
                matrix_increment_row(B,i);
                goto try_again;
            }
        }
    }


    return count;
}
