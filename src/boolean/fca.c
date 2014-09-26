/**
 * fca.c, (c) 2014, Immanuel Albrecht; Dresden University of
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

#include "fca.h"


/**
 * take some input intent B\subseteqM, calculate B''.
 */
void vector_close_intent(vector input, matrix context, vector output) {
    assert(output->size == context->cols);
    assert(input->size == context->cols);
    int i;

    vector_ones(output);
    for (i = 0; i < context->rows; i++) {
        if (vector_AsubsetB_row(input,context,i)) {
            vector_andB_row(output,context,i);
        }
    }
//    vector_fput(input,stdout);
//    printf("-->");
//    vector_fput(output,stdout);
}


/**
 *
 * run the next closure algorithm on the given context,
 * call the callback for every attribute-vector that belongs
 * to a concept, and hand over user parameter
 *
 * returns 1 if aborted by the callback, 0 else
 */
int matrix_next_closure(matrix ctx, cb_next_closure callback, void *user)
{
    int i, j, good;

    assert(ctx);

    vector M = vector_alloc(ctx->cols);
    vector Y = vector_alloc(ctx->cols);
    vector DELTA;

	/**
	 * calculate the bottom intent of the concept lattice, i.e. {}''
	 */
	vector_close_intent(Y,ctx, M);

    if (callback(M,user,0))
    {
        vector_free(M);
        vector_free(Y);
        return 1;
    }

	/**
	 * begin of nextClosure function iteration
	 */
	nextClosure:

	for (i = ctx->cols; i > 0;)
	{
		--i;

		if (!INCIDESV(M->incidence,i))
		{
			CROSSV(M->incidence, i);
			vector_close_intent(M, ctx, Y);

			good = 1;

			for (j = 0; j < OFFSET(i); ++j)
			{
				if (Y->incidence[j] & (~(M->incidence[j])))
				{
					good = 0;
					break;
				}
			}
			if (good)
			{
				if (Y->incidence[OFFSET(i)] & (~M->incidence[OFFSET(i)]) & CRIMPVALUE(i))
				{
					good = 0;
				}
			}

			if (good)
			{
				/**
				 * we found the next intent
				 */

                if (callback(M,user,0)) {
                    vector_free(M);
                    vector_free(Y);
                    return 1;
                }

				/**
				 * continue with Y for M
				 */

				DELTA = M;
				M = Y;
				Y = DELTA;
				/**
				 * do the nextClosure
				 */
				goto nextClosure;
			}
		}

		CLEARV(M->incidence, i);
	}

	/**
	 * free up memory
	 */

	vector_free(M);
	vector_free(Y);

	return 0;
}

static int mcc_cb(vector V, void *user, int thread) {
    int *pI = user;
    (*pI) ++;

    return 0;
}

int matrix_count_concepts(matrix ctx) {
    int count;
    count = 0;
    matrix_next_closure(ctx,mcc_cb,&count);
    return count;
}

static int mccm_cb(vector V, void *user, int thread) {
    int *pI = user;
    int sum = vector_checksum(V);
    
    if (sum > *pI)
        *pI = sum;

    return 0;
}

int matrix_concept_checksum_max(matrix ctx) {
    int max;
    max = 0;
    matrix_next_closure(ctx,mccm_cb,&max);
    return max;
}

static int mcf_cb(vector V, void *user, int thread) {
    vector_fput(V,user);
    return 0;
}
void matrix_concept_fput(matrix ctx, FILE* f) {
    matrix_next_closure(ctx,mcf_cb,f);
}
