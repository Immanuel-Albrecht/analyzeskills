/**
 * analyzeskills.c, (c) 2014, Immanuel Albrecht; Dresden University of
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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <err.h>

#include <cprops/heap.h>
#include <cprops/util.h>
#include <cprops/hashtable.h>
#include <cprops/vector.h>
#include <cprops/multimap.h>
#include <cprops/trie.h>
#include <cprops/avl.h>
#include <cprops/rb.h>

#include "boolean/advanced.h"

uint64_t largest_antichain(int N) {
    int Nhalf = N / 2;
    int i;
    uint64_t product = 1;
    uint64_t quotient = 1;
    for (i = 1; i <= Nhalf; i++) {
        product *= (N+i-1);
        quotient *= i;
    }
    return product/quotient;
}

/** thanks to http://rosettacode.org/wiki/Read_a_file_line_by_line#Using_mmap.28.29 */

int read_lines(const char * fname, int (*call_back)(const char*, const char*))
{
        int fd = open(fname, O_RDONLY);
        struct stat fs;
        char *buf, *buf_end;
        char *begin, *end, c;
 
        if (fd == -1) {
                err(1, "open: %s", fname);
                return 0;
        }
 
        if (fstat(fd, &fs) == -1) {
                err(1, "stat: %s", fname);
                return 0;
        }
 
        /* fs.st_size could have been 0 actually */
        buf = mmap(0, fs.st_size, PROT_READ, MAP_SHARED, fd, 0);
        if (buf == (void*) -1) {
                err(1, "mmap: %s", fname);
                close(fd);
                return 0;
        }
 
        buf_end = buf + fs.st_size;
 
        begin = end = buf;
        while (1) {
                if (! (*end == '\r' || *end == '\n')) {
                        if (++end < buf_end) continue;
                } else if (1 + end < buf_end) {
                        /* see if we got "\r\n" or "\n\r" here */
                        c = *(1 + end);
                        if ( (c == '\r' || c == '\n') && c != *end)
                                ++end;
                }
 
                /* call the call back and check error indication. Announce
                   error here, because we didn't tell call_back the file name */
                if (! call_back(begin, end)) {
                        err(1, "[callback] %s", fname);
                        break;
                }
 
                if ((begin = ++end) >= buf_end)
                        break;
        }
 
        munmap(buf, fs.st_size);
        close(fd);
        return 1;
}

int rows,cols;

matrix M;

static int get_size(const char* begin, const char* end){
    int i = end - begin;
    rows ++;
    if (i>cols)
        cols = i;
}

int gd_row;

static int get_data(const char* begin, const char* end){
    int i,N;
    if (gd_row >= rows)
        return;

    N = cols;
    if (end-begin < cols)
        N = end-begin;

    for (i = 0; i < N; i++) {
        switch(begin[i]) {
            case 'X':
            case '1':
            case 'x':
            case '+':
                CROSSV(ROW(gd_row,M),i);
        }
    }
    gd_row ++;
}


int main(int argc, const char *argv[])
{
    int dnf_skills;
    int run;
    char buf[1025];
    int i,j;
    int incomparables;
    uint64_t antichain;

    if (argc != 2) {
        puts("Usage: analyzeskills [response_patterns.txt]");
        exit(0);
    }

    rows = 0;
    cols = 0;
    gd_row = 0;

    read_lines(argv[1],get_size);
    printf("Response Dimensions: %d Subjects x %d Items\n",rows,cols);

    M = matrix_alloc(rows,cols);

    read_lines(argv[1],get_data);

    puts("M=");
    matrix_fput(M,stdout);

    matrix INCMP=matrix_alloc(M->rows,M->rows);
    vector V,W;
    for (i = 0; i < M->rows; i++) {
        CROSSV(ROW(i,INCMP),i);
        V = vector_dup_matrix_row(M,i);
        for (j = i+1; j < M->rows; j++) {
            W = vector_dup_matrix_row(M,j);
            
            if ((! vector_AsubsetB(V,W))&&(! vector_AsubsetB(W,V)))
            {
                CROSSV(ROW(i,INCMP),j);
                CROSSV(ROW(j,INCMP),i);
            }


            vector_free(W);
        }
        vector_free(V);
    }

    incomparables = matrix_concept_checksum_max(INCMP);

    printf("Need %d incomparable skill vectors.\n",incomparables);



    matrix_free(INCMP);





    matrix A,B;
    positive_DNF_matrix L;
    matrix S,C;

    generate_skill_matrix(&A,&B,M);

    C = matrix_dupAxBT_T(A,B);
    assert(matrix_cmp(M,C)==0);
    matrix_free(C);

    printf("Relation Product Factorization A*B=M: %d Factors\n",B->cols);

    puts("A=");
    matrix_fput(A,stdout);
    puts("B=");
    matrix_fput(B,stdout);


    dnf_skills = B->cols -1;
    run = 0;
    while(dnf_skills) {
        ++run;
        printf("Searching for pDNF L%d and skill matrix S%d with %d elements...\n",run,run,dnf_skills);
        antichain = largest_antichain(dnf_skills);
        printf("Largest antichain in %d-elementary powerset: %"PRIu64"\n",dnf_skills, antichain);
        if (antichain < incomparables)
        {
            printf("Too many incomparable responses for %d skills!\n",incomparables);
            break;
        }
        
        L = positive_DNF_matrix_alloc(A->rows,dnf_skills);
        S = matrix_alloc(B->rows, dnf_skills);
        if (!find_DNF_B_from_DNFxBTT(L,S,M)) {
            puts("None found.");

            matrix_free(S);
            DNF_free(L);
            break;
        }

        C = matrix_dup_DNFxBT_T(L,S);
        assert(matrix_cmp(M,C)==0);
        matrix_free(C);

        snprintf(buf,1024,"L%d",run);
        DNF_fput(buf,L,stdout);
        printf("S%d=\n",run);
        matrix_fput(S,stdout);

        matrix_free(S);
        DNF_free(L);
        dnf_skills--;
        fflush(stdout);
    }



    matrix_free(A);
    matrix_free(B);
    matrix_free(M);
    return 0;
}


