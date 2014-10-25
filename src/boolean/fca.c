/*
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


/** @brief calculate .'' for an intent-bit-vector
 *
 * @param output  <-- input''
 * @param input  attributes to be closed
 * @param context formal context matrix I
 *
 *
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


/** @brief run next_closure algorithm on ctx, uses callback
 *
 * @param ctx the context represented as bit-matrix
 * @param callback routine called on each intent-vector that
 *                 belongs to a concept of ctx.
 * @param user user data pointer
 *
 * run the next closure algorithm on the given context,
 * call the callback for every attribute-vector that belongs
 * to a concept, and hand over user parameter
 *
 * @returns 1 if aborted by the callback, 0 else
 */
int matrix_next_closure(matrix ctx, cb_next_closure callback, void *user)
{
    int i, j, good;

    assert(ctx);

    vector M = vector_alloc(ctx->cols);
    vector Y = vector_alloc(ctx->cols);
    vector DELTA;

	/*
	 * calculate the bottom intent of the concept lattice, i.e. {}''
	 */
	vector_close_intent(Y,ctx, M);

    if (callback(M,user,0))
    {
        vector_free(M);
        vector_free(Y);
        return 1;
    }

	/*
	 * begin of nextClosure function iteration
	 */
	nextClosure:

	for (i = ctx->cols; i > 0;)
	{
        //puts("M=");
        //vector_fput(M,stdout);
		--i;

		if (!INCIDESV(M->incidence,i))
		{
			CROSSV(M->incidence, i);
            //puts("M+i=");
            //vector_fput(M,stdout);

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
                //puts("Y=");
                //vector_fput(Y,stdout);
				if (Y->incidence[OFFSET(i)] 
                        & (~ M->incidence[OFFSET(i)])
                        & CRIMPVALUE(i))
				{
					good = 0;
                    //printf("%"PRIx64" CRIMP\n",CRIMPVALUE(i));
                    //printf("%"PRIx64" Y\n",Y->incidence[OFFSET(i)]);
                    //printf("%"PRIx64" M\n",M->incidence[OFFSET(i)]);
				}
			}

			if (good)
			{
				/*
				 * we found the next intent
				 */

                if (callback(Y,user,0)) {
                    vector_free(M);
                    vector_free(Y);
                    return 1;
                }

				/*
				 * continue with Y for M
				 */

				DELTA = M;
				M = Y;
				Y = DELTA;
				/*
				 * do the nextClosure
				 */
				goto nextClosure;
			}
		}

		CLEARV(M->incidence, i);
	}

	/*
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
/** @brief uses next_closure to count the number of concepts
 *
 *
 *
 */
int matrix_count_concepts(matrix ctx) {
    int count;
    count = 0;
    matrix_next_closure(ctx,mcc_cb,&count);
    return count;
}

static int mccm_cb(vector V, void *user, int thread) {
    int *pI = user;

    if (vector_isOnes(V)) /** ignore the bottom vector */
        return 0;

    int sum = vector_checksum(V);
    
    if (sum > *pI)
        *pI = sum;

    return 0;
}

/**
 * @brief max. number of set attribute bits for atoms of the concept lattice
 *
 * @returns the maximal number of set attributes of the 
 * atoms of the concept lattice.
 *
 * probably pretty dumb to go through all concepts for this.
 */
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

/** @brief print all concepts of ctx using next_closure
 *
 *
 *
 */
void matrix_concept_fput(matrix ctx, FILE* f) {
    matrix_next_closure(ctx,mcf_cb,f);
}
/** @brief like matrix_next_closure for an ideal
 *
 * @param attributes set that generates the ideal 
 *                   (i.e. subset of allowed concept intents)
 *
 * @param ctx the context represented as bit-matrix
 * @param callback routine called on each intent-vector that
 *                 belongs to a concept of ctx.
 * @param user user data pointer
 *
 * 
 * next closure that only returns intents that
 * have at least the given attributes
 *
 * @returns 1 if aborted by the callback, 0 else
 */
int matrix_next_closure_ideal(matrix ctx, vector attributes,
        cb_next_closure callback, void *user)
{
    int i, j, good;

    assert(attributes->size == ctx->cols);

    assert(ctx);

    vector M0 = vector_alloc(ctx->cols);
    vector M,Y;
    vector DELTA;

	/*
	 * calculate the bottom intent of the concept sub-lattice, i.e. {}''
	 */
	vector_close_intent(attributes,ctx, M0);

    if (callback(M0,user,0))
    {
        vector_free(M0);
        return 1;
    }

    /* initialize M and Y */

    Y = vector_alloc(ctx->cols);
    M = vector_dup(M0);

	/*
	 * begin of nextClosure function iteration
	 */
	nextClosure:

	for (i = ctx->cols; i > 0;)
	{
        //puts("M=");
        //vector_fput(M,stdout);
		--i;

        if (INCIDESV(M0->incidence, i))
            /* keep the attributes of M0 */
            continue; 

		if (!INCIDESV(M->incidence,i))
		{
			CROSSV(M->incidence, i);
            //puts("M+i=");
            //vector_fput(M,stdout);

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
                //puts("Y=");
                //vector_fput(Y,stdout);
				if (Y->incidence[OFFSET(i)] 
                        & (~ M->incidence[OFFSET(i)])
                        & CRIMPVALUE(i))
				{
					good = 0;
                    //printf("%"PRIx64" CRIMP\n",CRIMPVALUE(i));
                    //printf("%"PRIx64" Y\n",Y->incidence[OFFSET(i)]);
                    //printf("%"PRIx64" M\n",M->incidence[OFFSET(i)]);
				}
			}

			if (good)
			{
				/*
				 * we found the next intent
				 */

                if (callback(Y,user,0)) {
                    vector_free(M0);
                    vector_free(M);
                    vector_free(Y);
                    return 1;
                }

				/*
				 * continue with Y for M
				 */

				DELTA = M;
				M = Y;
				Y = DELTA;
				/*
				 * do the nextClosure
				 */
				goto nextClosure;
			}
		}

		CLEARV(M->incidence, i);
	}

	/*
	 * free up memory
	 */

    vector_free(M0);
	vector_free(M);
	vector_free(Y);

	return 0;
}

/** @brief like matrix_next_closure, but for some filter
 *
 * @param attributes  attributes that defines the filter 
 *                    (i.e. superset of allowed concept intents)
 *
 * @param ctx the context represented as bit-matrix
 * @param callback routine called on each intent-vector that
 *                 belongs to a concept of ctx.
 * @param user user data pointer
 *
 * run the next closure algorithm on the given context,
 * call the callback for every attribute-vector that belongs
 * to a concept, and hand over user parameter
 *
 * @returns 1 if aborted by the callback, 0 else
 */
int matrix_next_closure_filter(matrix ctx, vector attributes, cb_next_closure callback, void *user)
{
    int i, j, good;

    assert(ctx);
    assert(ctx->cols == attributes->size);

    vector M = vector_alloc(ctx->cols);
    vector Y = vector_alloc(ctx->cols);
    vector DELTA;

    /* store the forbidden attributes set */
    vector invM1 = vector_dup(attributes);
    vector_Inv(invM1);

	/*
	 * calculate the bottom intent of the concept lattice, i.e. {}''
	 */
	vector_close_intent(Y,ctx, M);

    /* check whether we are still in the filter */
    if (! vector_AdisjointB(invM1, M)) {
        vector_free(M);
        vector_free(Y);
        vector_free(invM1);

        return 0;
    }
    

    if (callback(M,user,0))
    {
        vector_free(M);
        vector_free(Y);
        vector_free(invM1);
        return 1;
    }

	/*
	 * begin of nextClosure function iteration
	 */
	nextClosure:

	for (i = ctx->cols; i > 0;)
	{
        //puts("M=");
        //vector_fput(M,stdout);
		--i;

        if (INCIDESV(invM1->incidence, i))
            continue;

		if (!INCIDESV(M->incidence,i))
		{
			CROSSV(M->incidence, i);
            //puts("M+i=");
            //vector_fput(M,stdout);

			vector_close_intent(M, ctx, Y);

            /* check whether there are attributes in Y that aren't in the
             * attributes-parameter
             */
			good = vector_AdisjointB(invM1,Y);

            if (good) {
                for (j = 0; j < OFFSET(i); ++j)
                {
                    if (Y->incidence[j] & (~(M->incidence[j])))
                    {
                        good = 0;
                        break;
                    }
                }
            }

			if (good)
			{
                //puts("Y=");
                //vector_fput(Y,stdout);
				if (Y->incidence[OFFSET(i)] 
                        & (~ M->incidence[OFFSET(i)])
                        & CRIMPVALUE(i))
				{
					good = 0;
                    //printf("%"PRIx64" CRIMP\n",CRIMPVALUE(i));
                    //printf("%"PRIx64" Y\n",Y->incidence[OFFSET(i)]);
                    //printf("%"PRIx64" M\n",M->incidence[OFFSET(i)]);
				}
			}

			if (good)
			{
				/*
				 * we found the next intent
				 */

                if (callback(Y,user,0)) {
                    vector_free(M);
                    vector_free(Y);
                    vector_free(invM1);
                    return 1;
                }

				/*
				 * continue with Y for M
				 */

				DELTA = M;
				M = Y;
				Y = DELTA;
				/*
				 * do the nextClosure
				 */
				goto nextClosure;
			}
		}

		CLEARV(M->incidence, i);
	}

	/*
	 * free up memory
	 */

	vector_free(M);
	vector_free(Y);
    vector_free(invM1);

	return 0;
}
/** @brief like matrix_next_closure, for a given interval
 *
 *
 * 
 *
 */
int matrix_next_closure_interval(matrix ctx, vector subset_of_attrs, vector superset_of_attrs, 
        cb_next_closure callback, void *user) 
{

    int i, j, good;

    assert(ctx);
    assert(ctx->cols == superset_of_attrs->size);
    assert(ctx->cols == subset_of_attrs->size);

    vector M0 = vector_alloc(ctx->cols);
	vector_close_intent(subset_of_attrs,ctx,M0);

	/*
	 * calculate the bottom intent of the concept sub-lattice, i.e. {}''
	 */
    vector M,Y;
    vector DELTA;

    /* store the forbidden superset_of_attrs set */
    vector invM1 = vector_dup(superset_of_attrs);
    vector_Inv(invM1);

	/*
	 * calculate the bottom intent of the concept lattice, i.e. {}''
	 */
	vector_close_intent(subset_of_attrs,ctx, M0);

    /* check whether we are still in the filter */
    if (! vector_AdisjointB(invM1, M0)) {
        vector_free(M0);
        vector_free(invM1);

        return 0;
    }

    if (callback(M0,user,0))
    {
        vector_free(M0);
        vector_free(invM1);
        return 1;
    }

    /* now, allocate the loop variables */
    Y = vector_alloc(ctx->cols);
    M = vector_dup(M0);

	/*
	 * begin of nextClosure function iteration
	 */
	nextClosure:

	for (i = ctx->cols; i > 0;)
	{
        //puts("M=");
        //vector_fput(M,stdout);
		--i;

        if ((INCIDESV(invM1->incidence, i))||(INCIDESV(M0->incidence,i)))
            continue;

		if (!INCIDESV(M->incidence,i))
		{
			CROSSV(M->incidence, i);
            //puts("M+i=");
            //vector_fput(M,stdout);

			vector_close_intent(M, ctx, Y);

            /* check whether there are superset_of_attrs in Y that aren't in the
             * superset_of_attrs-parameter
             */
			good = vector_AdisjointB(invM1,Y);

            if (good) {
                for (j = 0; j < OFFSET(i); ++j)
                {
                    if (Y->incidence[j] & (~(M->incidence[j])))
                    {
                        good = 0;
                        break;
                    }
                }
            }

			if (good)
			{
                //puts("Y=");
                //vector_fput(Y,stdout);
				if (Y->incidence[OFFSET(i)] 
                        & (~ M->incidence[OFFSET(i)])
                        & CRIMPVALUE(i))
				{
					good = 0;
                    //printf("%"PRIx64" CRIMP\n",CRIMPVALUE(i));
                    //printf("%"PRIx64" Y\n",Y->incidence[OFFSET(i)]);
                    //printf("%"PRIx64" M\n",M->incidence[OFFSET(i)]);
				}
			}

			if (good)
			{
				/*
				 * we found the next intent
				 */

                if (callback(Y,user,0)) {
                    vector_free(M);
                    vector_free(Y);
                    vector_free(M0);
                    vector_free(invM1);
                    return 1;
                }

				/*
				 * continue with Y for M
				 */

				DELTA = M;
				M = Y;
				Y = DELTA;
				/*
				 * do the nextClosure
				 */
				goto nextClosure;
			}
		}

		CLEARV(M->incidence, i);
	}

	/*
	 * free up memory
	 */

	vector_free(M);
	vector_free(M0);
	vector_free(Y);
    vector_free(invM1);

	return 0;
}
