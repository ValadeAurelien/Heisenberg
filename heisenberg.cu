#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#define PI 3.141592653


/* 
 * DISTMAT
 */

typedef struct {
    int N;
    int a;
    float *d_distmat,
          *h_distmat;
} distmat;

__host__ 
float get_h_dist(distmat *d, int i, int j, int k, int l) {
    return d->h_distmat[ abs(i-k)*d->N+abs(j-l) ];
}

__device__ 
float get_d_dist(distmat *d, int i, int j, int k, int l) {
    return d->d_distmat[ abs(i-k)*d->N+abs(j-l) ];
}

__host__
int init_distmat(distmat *d, int N, int a) {
    if ( N == 0 ) {
        printf("init distmat : N is null");
        return 1;
    }
    d->N = N;
    d->a = a;
    d->h_distmat = (float*) malloc(sizeof(float)*d->N*d->N);
    cudaMalloc(&(d->d_distmat), sizeof(float)*d->N*d->N);
    for(int i=0; i<d->N; ++i) 
        for(int j=0; j<d->N; ++j) 
            if (i==0 && j==0)
                d->h_distmat[i*d->N+j] = 0;
            else
                d->h_distmat[i*d->N+j] = pow((float) i*i+j*j, -d->a/2.);
    cudaMempy(d->d_distmat, d->h_distmat, sizeof(float)*d->N*d->N, cudaMemcpyHostToDevice);
    /*for(int i=0; i<d->N; ++i) {*/
        /*printf("\n");*/
        /*for(int j=0; j<d->N; ++j) */
            /*printf("%f\t\t", get_dist(d, 0, 0, i, j));*/
    /*}*/

    return 0;
}

/*
 * LATTICE
 */

typedef enum {NONE=0, RANDOM=1, UP=2, DOWN=3, CHECKERBOARD=4} init_type;

typedef enum {LATTICE=1, SCALARS=2} print_type;

typedef struct {
    int N, neighb;
    float k, h;
    float *h_lattice, *d_lattice, 
          *h_latcos, *d_latcos,
          *h_latsin, *d_latsin,
          *h_latinter, *d_latinter;
    float Q, M, E;
} lattice;


__host__
void set_h_lat(lattice *l, int i, int j, float val) {
    l->h_lattice[i*l->N+j] = val;
}

__device__
void set_d_lat(lattice *l, int i, int j, float val) {
    l->d_lattice[i*l->N+j] = val;
}

__host__
float get_h_lat(lattice *l, int i, int j) {
    return l->h_lattice[i*l->N+j];
}

__device__
float get_d_lat(lattice *l, int i, int j) {
    return l->d_lattice[i*l->N+j];
}

__host__
void set_h_cos(lattice *l, int i, int j, float val) {
    l->h_latcos[i*l->N+j] = val;
}

__device__
void set_d_cos(lattice *l, int i, int j, float val) {
    l->d_latcos[i*l->N+j] = val;
}

__host__
float get_h_cos(lattice *l, int i, int j) {
    return l->h_latcos[i*l->N+j];
}

__device__
float get_d_cos(lattice *l, int i, int j) {
    return l->d_latcos[i*l->N+j];
}

__host__
void set_h_sin(lattice *l, int i, int j, float val) {
    l->h_latsin[i*l->N+j] = val;
}

__device__
void set_d_sin(lattice *l, int i, int j, float val) {
    l->d_latsin[i*l->N+j] = val;
}

__host__
float get_h_sin(lattice *l, int i, int j) {
    return l->h_latsin[i*l->N+j];
}

__device__
float get_d_sin(lattice *l, int i, int j) {
    return l->d_latsin[i*l->N+j];
}

__host__
void set_h_inter(lattice *l, int i, int j, float val) {
    l->h_latinter[i*l->N+j] = val;
}

__device__
void set_d_inter(lattice *l, int i, int j, float val) {
    l->d_latinter[i*l->N+j] = val;
}

__host__
void add_h_inter(lattice *l, int i, int j, float val) {
    l->h_latinter[i*l->N+j] += val;
} 

__device__
void add_d_inter(lattice *l, int i, int j, float val) {
    l->d_latinter[i*l->N+j] += val;
} 

__host__
float get_h_inter(lattice *l, int i, int j) {
    return l->h_latinter[i*l->N+j];
}

__device__
float get_d_inter(lattice *l, int i, int j) {
    return l->d_latinter[i*l->N+j];
}

/*void calc_cossin(lattice *l) {*/
    /*for(int i=0; i<l->N; ++i) */
        /*for(int j=0; j<l->N; ++j) {*/
            /*set_cos(l, i, j, cos(get_lat(l, i, j)));*/
            /*set_sin(l, i, j, sin(get_lat(l, i, j)));*/
        /*}*/
/*}*/


__global__
void calc_cossin(lattice *l) {
    int i = blockIdx.x*blockDim.x+threadIdx.x,
        j = blockIdx.y*blockDim.y+threadIdx.y;
    if (i>l->N || j>l->N) return; 
    set_d_cos(l, i, j, cosf(get_lat(l, i, j)));
    set_d_sin(l, i, j, sinf(get_lat(l, i, j)));
}

__host__
int init_lattice(lattice *l, int N, float h, float k, int neighb, init_type it) {
    if ( N == 0 ) {
        printf("init lattice : N is null");
        return 1;
    }
    l->N = N;
    l->h = h;
    l->k = k;
    l->neighb = neighb;
    l->h_lattice  = (float*) malloc(sizeof(float)*l->N*l->N);
    l->h_latcos   = (float*) malloc(sizeof(float)*l->N*l->N);
    l->h_latsin   = (float*) malloc(sizeof(float)*l->N*l->N);
    l->h_latinter = (float*) malloc(sizeof(float)*l->N*l->N);
    cudaMalloc(&(l->d_lattice),  sizeof(float)*l->N*l->N);
    cudaMalloc(&(l->d_latcos),   sizeof(float)*l->N*l->N);
    cudaMalloc(&(l->d_latsin),   sizeof(float)*l->N*l->N);
    cudaMalloc(&(l->d_latinter), sizeof(float)*l->N*l->N);
    l->Q = 0; l->M = 0; l->E = 0;
    if (l->h_lattice == NULL || 
            l->h_latcos == NULL ||
            l->h_latsin == NULL ||
            l->h_latinter == NULL ) {
        printf("init lattice : cannot allocate memory");
        return 1;
    }
    if (l->d_lattice == NULL || 
            l->d_latcos == NULL ||
            l->d_latsin == NULL ||
            l->d_latinter == NULL ) {
        printf("init lattice : cannot allocate memory (GPU)");
        return 1;
    }
    switch (it) {
        case NONE:
            break;
        case RANDOM:
            for(int i=0; i<l->N; ++i) 
                for(int j=0; j<l->N; ++j) 
                    set_h_lat(l, i, j, rand()*PI/RAND_MAX);
            break;
        case UP:
            for(int i=0; i<l->N; ++i) 
                for(int j=0; j<l->N; ++j) 
                    set_h_lat(l, i, j, 0);
            break;
        case DOWN:
            for(int i=0; i<l->N; ++i) 
                for(int j=0; j<l->N; ++j) 
                    set_h_lat(l, i, j, PI);
            break;
        case CHECKERBOARD:
            for(int i=0; i<l->N; i=i+2) 
                for(int j=0; j<l->N; j=j+2) 
                    set_h_lat(l, i, j, PI);
            for(int i=1; i<l->N; i=i+2) 
                for(int j=1; j<l->N; j=j+2) 
                    set_h_lat(l, i, j, PI);
            for(int i=1; i<l->N; i=i+2) 
                for(int j=0; j<l->N; j=j+2) 
                    set_h_lat(l, i, j, 0);
            for(int i=0; i<l->N; i=i+2) 
                for(int j=1; j<l->N; j=j+2) 
                    set_h_lat(l, i, j, 0);
            break;
    }
    cudaMempy(l->d_lattice, l->h_lattice, sizeof(float)*d->N*d->N, cudaMemcpyHostToDevice);

    return 0;
}

/*void calc_inter_one_spin(lattice *l, distmat *d, int i, int j) {*/
    /*set_inter(l, i, j, 0);*/
    /*int sii, sjj;*/
    /*for(int ii=i-l->neighb; ii<i+l->neighb+1; ii++) */
        /*for(int jj=j-l->neighb; jj<j+l->neighb+1; jj++) */
            /*if (!(ii == i && jj == j)) {*/
                /*if (ii<0) sii = ii + l->N;*/
                /*else if (ii>l->N) sii = ii - l->N;*/
                /*else sii = ii;*/
                /*if (jj<0) sjj = jj + l->N;*/
                /*else if (jj>l->N) sjj = jj - l->N;*/
                /*else sjj = jj;*/
                /*add_inter(l, i, j, get_dist(d, i, j, ii, jj)*/
                        /** ( get_cos(l, i, j)*get_cos(l, sii, sjj) +*/
                            /*get_sin(l, i, j)*get_sin(l, sii, sjj) ));*/
                /*[>printf("#%i %i %i %i %f %f %f %f %f %f\n", i, j, ii, jj, <]*/
                        /*[>get_dist(d, i, j, ii, jj), get_cos(l, i, j), get_cos(l, ii, jj),<]*/
                        /*[>get_sin(l, i, j), get_sin(l, ii, jj),<]*/
                        /*[>get_inter(l, i, j));<]*/
                /*}*/
/*}*/

/*void calc_inters(lattice *l, distmat *d) {*/
    /*for(int i=0; i<l->N; ++i) */
        /*for(int j=0; j<l->N; ++j) */
            /*calc_inter_one_spin(l, d, i, j);*/
/*}*/

__global__
void calc_inters(lattice *l, distmat *d) {
    int i = blockIdx.x*blockDim.x+threadIdx.x,
        j = blockIdx.y*blockDim.y+threadIdx.y;
    if (i>l->N || j>l->N) return; 
    set_d_inter(l, i, j, 0);
    int sii, sjj;
    for(int ii=i-l->neighb; ii<i+l->neighb+1; ii++) 
        for(int jj=j-l->neighb; jj<j+l->neighb+1; jj++) 
            if (!(ii == i && jj == j)) {
                if (ii<0) sii = ii + l->N;
                else if (ii>l->N) sii = ii - l->N;
                else sii = ii;
                if (jj<0) sjj = jj + l->N;
                else if (jj>l->N) sjj = jj - l->N;
                else sjj = jj;
                add_d_inter(l, i, j, get_d_dist(d, i, j, ii, jj)
                        * ( get_d_cos(l, i, j)*get_d_cos(l, sii, sjj) +
                            get_d_sin(l, i, j)*get_d_sin(l, sii, sjj) ));
                /*printf("#%i %i %i %i %f %f %f %f %f %f\n", i, j, ii, jj, */
                        /*get_dist(d, i, j, ii, jj), get_cos(l, i, j), get_cos(l, ii, jj),*/
                        /*get_sin(l, i, j), get_sin(l, ii, jj),*/
                        /*get_inter(l, i, j));*/
                }
}

__host__
void calc_Q(lattice *l) {
    l->Q = 0;
    for(int i=0; i<l->N; ++i) 
        for(int j=0; j<l->N; ++j) {
            //printf("%u %u %f\n", i, j, get_inter(l, i, j));
            l->Q += get_h_inter(l, i, j);
        }
}

__host__
void calc_M(lattice *l) {
    l->M = 0;
    for(int i=0; i<l->N; ++i) 
        for(int j=0; j<l->N; ++j) 
            l->M += get_h_cos(l, i, j);
}

__host__
void calc_E(lattice *l) {
    l->E = -( l->k*l->Q + l->h*l->M );
}

__host__
void calc_all(lattice *l, distmat *d, dim3 blockSize, dim3 gridSize) {
    calc_cossin<<<blockSize, gridSize>>>(l);
    calc_inters<<<blockSize, gridSize>>>(l, d);
    cudaMempy(l->h_lattices, l->d_lattices, sizeof(float)*d->N*d->N, cudaMemcpyDeviceToHost);
    cudaMempy(l->h_latcos,   l->d_latcos,   sizeof(float)*d->N*d->N, cudaMemcpyDeviceToHost);
    cudaMempy(l->h_latinter, l->d_latinter, sizeof(float)*d->N*d->N, cudaMemcpyDeviceToHost);
    calc_Q(l);
    calc_M(l);
    calc_E(l);
}

__host__
void print_lat(lattice *l, int n, int print_type,
               char *latfname, char*intfname, 
               char *ext, FILE *scaf) {
    char _latfname[32], 
         _intfname[32],
         number[16], header[128];
    memset(_latfname, 0, 32);
    memset(_intfname, 0, 32);
    memset(number, 0, 16);
    switch(print_type) {
        case LATTICE:
            sprintf(number, "%04u", n);
            strcat(_latfname, latfname);
            strcat(_latfname, number);
            strcat(_latfname, ext);
            FILE *latf = fopen(_latfname, "w");
            fprintf(latf, "# %u %u %u %f %f %f %f %f\n", n, l->N, l->neighb, 
                    l->h, l->k, l->M, l->Q, l->E);
            for(int i=0; i<l->N; ++i) {
                fprintf(latf, "\n");
                for(int j=0; j<l->N; ++j) 
                    fprintf(latf, "\t\t%.5f", get_lat(l, i, j));
            }
            fclose(latf);
            strcat(_intfname, intfname);
            strcat(_intfname, number);
            strcat(_intfname, ext);
            FILE *intf = fopen(_intfname, "w");
            fprintf(intf, "# %u %u %u %f %f %f %f %f\n", n, l->N, l->neighb, 
                    l->h, l->k, l->M, l->Q, l->E);
            for(int i=0; i<l->N; ++i) {
                fprintf(intf, "\n");
                for(int j=0; j<l->N; ++j) 
                    fprintf(intf, "\t\t%.5f", get_inter(l, i, j));
            }
            fclose(intf);
            break;
        case SCALARS:
            if (n==0) fprintf(scaf, "# %u %u %f %f \n", l->N, l->neighb, l->h, l->k);
            fprintf(scaf, "%03u %.10f %.10f %.10f\n", n, l->M, l->Q, l->E);
            break;
    }
}

/*
 * HEISENBERG
 */

typedef enum {METROPOLIS=1, GLAUBER=2} MC_type;

typedef struct {
    lattice *new_l,
            *old_l;
    distmat *d;
    int mct;
} heisenberg;

__host__
int init_heisenberg(heisenberg *he, 
        int N, float h, float k, int neighb, int a, 
        init_type it, int mct) {
    he->mct = mct;
    he->new_l = (lattice*) malloc(sizeof(lattice));
    he->old_l = (lattice*) malloc(sizeof(lattice));
    he->d     = (distmat*) malloc(sizeof(distmat));
    if (he->old_l == NULL || he->new_l == NULL || he->d == NULL) {
        printf("init heisenberg, could not alloc memory");
        return 1;
    }
    if (init_lattice(he->old_l, N, h, k, neighb, it)) return 1;
    if (init_lattice(he->new_l, N, h, k, neighb, NONE)) return 1;
    if (init_distmat(he->d, N, a)) return 1; 
    return 0;
}

__host__ 
void swap_d_lat(lattice *l, int i, int j, float *oldval, float newval) {
    cudaMemcpy(oldval, &(l->lattice[i*l->N + j]), sizeof(float), cudaMemcpyDeviceToHost); 
    cudaMemcpy(&(l->lattice[i*l->N + j]), &newval, sizeof(float), cudaMemcpyHostToDevice); 
}

__host__
void accept_step(heisenberg *h, int i, int j, float angle) {
    lattice *tmp;
    tmp = h->old_l; h->old_l = h->new_l; h->new_l = tmp;
    /*set_lat(h->new_l, i, j, angle);*/
    swap_d_lat(h->new_l, i, j, &(tmp->M), angle); //random float available
}

__host__
void update(heisenberg *h, dim3 blockSize, dim3 gridSize) {
    float angle = rand()*PI/RAND_MAX,
          old_angle;
    int i = floor((float) rand()/RAND_MAX*h->new_l->N),
        j = floor((float) rand()/RAND_MAX*h->new_l->N);
    old_angle = get_lat(h->new_l, i, j);
    swap_d_lat(h->new_l, i, j, old_angle, angle)
    float acc = (float) rand()/RAND_MAX,
          dE;

    calc_all(h->new_l, h->d, blockSize, gridSize);
    dE = h->new_l->E - h->old_l->E;
    switch(h->mct) {
        case METROPOLIS:
            if (dE < 0 || acc < exp(-dE))
                accept_step(h, i, j, angle);
            else 
                swap_d_lat(h->new_l, i, j, angle, old_angle);
            break;
        case GLAUBER:
            if (acc < 1/(1+exp(dE))) 
                accept_step(h, i, j, angle);
            else 
                swap_d_lat(h->new_l, i, j, angle, old_angle);
            break;
    }
}


int main(int argc, char **argv) {
    if (argc < 7) {
        printf("missing args");
        return 1;
    }

    int N          = atoi(argv[1]);
    float h        = atof(argv[2]),
          k        = atof(argv[3]);
    int neighb     = atoi(argv[4]),
        a          = atoi(argv[5]),
        tmax       = atoi(argv[6]),
        t_lat      = atoi(argv[7]),
        t_sca      = atoi(argv[8]);
    char *latfname = argv[9],
         *intfname = argv[10],
         *scafname = argv[11],
         *ext      = argv[12],
         _scafname[32];
 
    dim3 blockSize (32, 32),
         gridSize  (32, 32);

    srand(time(NULL));
    
    memset(_scafname, 0, 32);
    strcat(_scafname, scafname);
    strcat(_scafname, ext);
    FILE *scaf = fopen(_scafname, "w"); 

    int t_percent =  floor(tmax/100.);
    heisenberg *he = (heisenberg*) malloc(sizeof(heisenberg));
    printf("#initialisation...\n");
    init_heisenberg(he, N, h, k, neighb, a, UP, METROPOLIS);

    cudaMemcpy(&(he->old_l->d_lattice), &(he->old_l->h_lattice), 
               sizeof(float)*he->old_l->N*he->old_l->N, cudaMemcpyHostToDevice);
    memcpy(he->new_l->h_lattice, he->old_l->h_lattice,
           sizeof(float)*he->old_l->h_lattice->N*he->old_l->h_latticet->N);
    cudaMemcpy(&(he->old_l->d_lattice), &(he->new_l->d_lattice), 
               sizeof(float)*he->old_l->N*he->old_l->N, cudaMemcpyDeviceToDevice);

    calc_all(he->old_l, he->d, blockSize, gridSize);

    print_lat(he->old_l, 0, LATTICE, latfname, intfname, ext, scaf);
    print_lat(he->old_l, 0, SCALARS, latfname, intfname, ext, scaf);
    printf("#running...\n");
    for (int t=1; t<tmax; ++t) {
        update(he);
        if (t%t_percent==0) { 
            printf("\r#%03u", (int) floor((float) t/tmax*100));
            fflush(stdout);
        }
        if (t%t_lat==0)     print_lat(he->old_l, t, LATTICE, latfname, intfname, ext, scaf);
        if (t%t_sca==0)     print_lat(he->old_l, t, SCALARS, latfname, intfname,  ext, scaf);
    }
    fclose(scaf);
    printf("\n");

    return 0;
}

// Qstom plot
