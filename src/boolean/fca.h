/*
 * fca.h, (c) 2014, Immanuel Albrecht; Dresden University of
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

#ifndef FCA_8UNGMSLQ

#define FCA_8UNGMSLQ

#include "advanced.h"

/**
 * Callback routine for next closure algorithms
 *
 * attributes __ attribute vector of the last found concept,
 *               DO NOT ALTER THE VALUE! (you could skip some concepts)
 *
 * user       __ pointer to user data structure
 * thread     __ thread id (=0 in single threaded version)
 *
 * return 0 to continue next_closure algorithm, 1 to abort.
 */
typedef int fn_cb_next_closure(vector attributes, void* user, int thread);
typedef fn_cb_next_closure *cb_next_closure;

void vector_close_intent(vector input, matrix context, vector output);

int matrix_next_closure(matrix ctx, cb_next_closure callback, void *user);
int matrix_next_closure_ideal(matrix ctx, vector attributes,
        cb_next_closure callback, void *user);
int matrix_next_closure_filter(matrix ctx, vector attributes, 
        cb_next_closure callback, void *user);
int matrix_next_closure_interval(matrix ctx, vector subset_of_attrs, vector superset_of_attrs, 
        cb_next_closure callback, void *user) ;
int matrix_count_concepts(matrix ctx);

int matrix_concept_checksum_max(matrix ctx);
void matrix_concept_fput(matrix ctx, FILE* f);

#endif /* end of include guard: FCA_8UNGMSLQ */
