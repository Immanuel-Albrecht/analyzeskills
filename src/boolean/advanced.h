/*
 * advanced.h, (c) 2014, Immanuel Albrecht; Dresden University of
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

#ifndef ADVANCED_LHKYG0F1

#define ADVANCED_LHKYG0F1

#include "matrices.h"

int matrix_solveB_AimpliesBeqC(matrix A, matrix B, matrix C);

cp_vector *find_join_base(cp_avltree *vector_set);

void generate_skill_matrix(matrix *pA, matrix *pB, matrix C);

cp_avltree* get_row_meet_closure(matrix C);


int solve_DNFxBT_T_eqC(positive_DNF_matrix DNF, matrix B, matrix C);
int solve_DNFxBT_T_eqC_row(positive_DNF_matrix DNF, matrix B, matrix C, int row);
int solve_DNF_from_BandC(positive_DNF_matrix DNF, matrix B, matrix C);

int partially_solve_DNF_from_BandC(positive_DNF_matrix DNF, matrix B, matrix C, int rows);
int find_DNF_B_from_DNFxBTT(positive_DNF_matrix DNF, matrix B, matrix C);

int count_Bs_from_DNFxBTT(positive_DNF_matrix DNF, matrix B, matrix C);
#endif /* end of include guard: ADVANCED_LHKYG0F1 */



