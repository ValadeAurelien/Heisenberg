#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#define PI 3.141592653
//#define INT
//#define SAFE


int safe_check(const char* func, int i, int j, int N) {
    if (i<0 || i>=N || j<0 || j>N) {
        printf("Error while reading table in function %s, with i=%i, j=%i, N=%i \n",
                func, i, j, N);
        exit(1);
    }
    return 0;
}

/* 
 * DISTMAT
 */

typedef struct {
    int N;
    int a;
    float *distmat;
    float integ;
} distmat;


float get_dist(const distmat *d, int i, int j, int k, int l) {
#ifdef SAFE
    safe_check("get_dist", abs(i-k), abs(j-l), d->N);
#endif
    return d->distmat[ (int) abs(i-k)*d->N+ (int) abs(j-l) ];
}

int init_distmat(distmat *d, int n, int a) {
    if ( n == 0 ) {
        printf("init distmat : n is null");
        return 1;
    }
    d->N = n+1;
    d->a = a;
    d->distmat = (float*) malloc(sizeof(float)*d->N*d->N);
    d->integ = 0;
    for(int i=0; i<d->N; ++i) 
        for(int j=0; j<d->N; ++j) {
            if (i==0 && j==0)
                d->distmat[i*d->N+j] = 0;
            else
                d->distmat[i*d->N+j] = pow((float) i*i+j*j, -d->a/2.);
        }
    for(int i=0; i<d->N; ++i) 
        for(int j=1; j<d->N; ++j) 
            d->integ += d->distmat[i*d->N+j];
    d->integ *= 4;
    for(int i=0; i<d->N; ++i) 
        for(int j=0; j<d->N; ++j) 
            d->distmat[i*d->N+j] /= d->integ;
    /*for(int i=0; i<d->N; ++i) {*/
        /*printf("\n");*/
        /*for(int j=0; j<d->N; ++j) */
            /*printf("%f\t\t", get_dist(d, 0, 0, i, j));*/
    /*}*/

    return 0;
}

/*
 * BACKUP LATTICE
 */

typedef struct {
    int N, os_i, os_j;
    float * inters;
    float angle,  
          cos, sin,  
          M, Q, E;
} backup_lattice;


/*
 * LATTICE
 */

typedef enum {NONE=0, RANDOM=1, UP=2, DOWN=3, CHECKERBOARD=4} init_type;

typedef enum {LATTICE=1, SCALARS=2} print_type;

typedef struct {
    int N, neighb;
    float k, h;
    float *lattice, 
          *latcos,
          *latsin,
          *latinter;
    float Q, M, E;
} lattice;


void set_lat(lattice *l, int i, int j, float val) {
#ifdef SAFE
    safe_check("set_lat", i, j, l->N);
#endif
    l->lattice[i*l->N+j] = val;
}

float get_lat(const lattice *l, int i, int j) {
#ifdef SAFE
    safe_check("get_lat", i, j, l->N);
#endif
    return l->lattice[i*l->N+j];
}

void set_cos(lattice *l, int i, int j, float val) {
#ifdef SAFE
    safe_check("set_cos", i, j, l->N);
#endif
    l->latcos[i*l->N+j] = val;
}

float get_cos(const lattice *l, int i, int j) {
#ifdef SAFE
    safe_check("get_cos", i, j, l->N);
#endif
    return l->latcos[i*l->N+j];
}

void set_sin(lattice *l, int i, int j, float val) {
#ifdef SAFE
    safe_check("set_sin", i, j, l->N);
#endif
    l->latsin[i*l->N+j] = val;
}

float get_sin(const lattice *l, int i, int j) {
#ifdef SAFE
    safe_check("get_sin", i, j, l->N);
#endif
    return l->latsin[i*l->N+j];
}

void set_inter(lattice *l, int i, int j, float val) {
#ifdef SAFE
    safe_check("set_inter", i, j, l->N);
#endif
    l->latinter[i*l->N+j] = val;
}

void add_inter(lattice *l, int i, int j, float val) {
#ifdef SAFE
    safe_check("add_inter", i, j, l->N);
#endif
    l->latinter[i*l->N+j] += val;
} 

float get_inter(const lattice *l, int i, int j) {
#ifdef SAFE
    safe_check("get_inter", i, j, l->N);
#endif
    return l->latinter[i*l->N+j];
}

void set_backup_inter(backup_lattice *bl, int i, int j, float val) {
#ifdef SAFE
    safe_check("set_backup_inter", i, j, bl->N);
#endif
    bl->inters[i*bl->N+j] = val;
}

float get_backup_inter(const backup_lattice *bl, int i, int j) {
#ifdef SAFE
    safe_check("get_backup_inter", i, j, bl->N);
#endif
    return bl->inters[i*bl->N+j];
}

void periodic_conditions(int i, int j, int* si, int* sj, int N) {
    if (i<0) *si = i + N;
    else if (i>=N) *si = i - N;
    else *si = i;
    if (j<0) *sj = j + N;
    else if (j>=N) *sj = j - N;
    else *sj = j;
}


void save_to_backup(backup_lattice *bl, const lattice *l, int i, int j) {
    bl->os_i = i;
    bl->os_j = j;
    bl->angle = get_lat(l, i, j);
    bl->cos = get_cos(l, i, j);
    bl->sin = get_sin(l, i, j);
    bl->M = l->M;
    bl->Q = l->Q;
    bl->E = l->E;
    int ii, jj, sii, sjj;
    for ( int si=0; si<bl->N; si++)
        for ( int sj=0; sj<bl->N; sj++) {
            ii = i - l->neighb + si;
            jj = j - l->neighb + sj;
            periodic_conditions(ii, jj, &sii, &sjj, l->N);
            set_backup_inter(bl, si, sj, get_inter(l, sii, sjj));
        }
}

void write_from_backup(lattice *l, const backup_lattice *bl) {
    set_lat(l, bl->os_i, bl->os_j, bl->angle);
    set_cos(l, bl->os_i, bl->os_j, bl->cos);
    set_sin(l, bl->os_i, bl->os_j, bl->sin);
    l->M = bl->M;
    l->Q = bl->Q;
    l->E = bl->E;
    int ii, jj, sii, sjj;
    for ( int si=0; si<bl->N; si++)
        for ( int sj=0; sj<bl->N; sj++) {
            ii = bl->os_i - l->neighb + si;
            jj = bl->os_j - l->neighb + sj;
            periodic_conditions(ii, jj, &sii, &sjj, l->N);
            set_inter(l, sii, sjj, get_backup_inter(bl, si, sj));
        }
}

void calc_cossin(lattice *l) {
    for(int i=0; i<l->N; ++i) 
        for(int j=0; j<l->N; ++j) {
            set_cos(l, i, j, cos(get_lat(l, i, j)));
            set_sin(l, i, j, sin(get_lat(l, i, j)));
        }
/*    for(int i=0; i<l->N; ++i) */
        /*for(int j=0; j<l->N; ++j) */
            /*printf("%i %i %f %f\n", i, j, get_cos(l, i, j), get_sin(l, i, j));*/
}

int init_backup_lattice(backup_lattice * bl, int n) {
    if ( n == 0 ) {
        printf("init backup_lattice : N is null");
        return 1;
    }
    bl->N = 2*n+1;
    bl->inters = (float*) malloc(bl->N*bl->N*sizeof(float));
    if ( bl->inters == NULL ) {
        printf("Init backup_lattice : cannot allocate memory !");
        return 1;
    }
    return 0;
}
    
int init_lattice(lattice *l, int N, float h, float k, int neighb, init_type it) {
    if ( N == 0 ) {
        printf("init lattice : N is null");
        return 1;
    }
    l->N = N;
    l->h = h;
    l->k = k;
    l->neighb = neighb;
    l->lattice  = (float*) malloc(sizeof(float)*l->N*l->N);
    l->latcos   = (float*) malloc(sizeof(float)*l->N*l->N);
    l->latsin   = (float*) malloc(sizeof(float)*l->N*l->N);
    l->latinter = (float*) malloc(sizeof(float)*l->N*l->N);
    l->Q = 0; l->M = 0; l->E = 0;
    if (l->lattice == NULL || 
            l->latcos == NULL ||
            l->latsin == NULL ||
            l->latinter == NULL) {
        printf("init lattice : cannot allocate memory");
        return 1;
    }
    switch (it) {
        case NONE:
            break;
        case RANDOM:
            for(int i=0; i<l->N; ++i) 
                for(int j=0; j<l->N; ++j) 
                    set_lat(l, i, j, rand()*PI/RAND_MAX);
            break;
        case UP:
            for(int i=0; i<l->N; ++i) 
                for(int j=0; j<l->N; ++j) 
                    set_lat(l, i, j, 0);
            break;
        case DOWN:
            for(int i=0; i<l->N; ++i) 
                for(int j=0; j<l->N; ++j) 
                    set_lat(l, i, j, PI);
            break;
        case CHECKERBOARD:
            for(int i=0; i<l->N; i=i+2) 
                for(int j=0; j<l->N; j=j+2) 
                    set_lat(l, i, j, PI);
            for(int i=1; i<l->N; i=i+2) 
                for(int j=1; j<l->N; j=j+2) 
                    set_lat(l, i, j, PI);
            for(int i=1; i<l->N; i=i+2) 
                for(int j=0; j<l->N; j=j+2) 
                    set_lat(l, i, j, 0);
            for(int i=0; i<l->N; i=i+2) 
                for(int j=1; j<l->N; j=j+2) 
                    set_lat(l, i, j, 0);
            break;
    }

    return 0;
}

void calc_inter_one_spin(lattice *l, const distmat *d, int i, int j) {
    set_inter(l, i, j, 0);
    int sii, sjj;
    for(int ii=i-l->neighb; ii<i+l->neighb+1; ii++) 
        for(int jj=j-l->neighb; jj<j+l->neighb+1; jj++) 
            if (!(ii == i && jj == j)) {
                periodic_conditions(ii, jj, &sii, &sjj, l->N);
                add_inter(l, i, j, get_dist(d, i, j, ii, jj)
                        * ( get_cos(l, i, j)*get_cos(l, sii, sjj) +
                            get_sin(l, i, j)*get_sin(l, sii, sjj) ));

                /*printf("#%i %i %i %i %f %f %f %f %f %f\n", i, j, ii, jj, */
                        /*get_dist(d, i, j, ii, jj), get_cos(l, i, j), get_cos(l, ii, jj),*/
                        /*get_sin(l, i, j), get_sin(l, ii, jj),*/
                        /*get_inter(l, i, j));*/
                }
}

void calc_inters(lattice *l, const distmat *d) {
    for(int i=0; i<l->N; ++i) 
        for(int j=0; j<l->N; ++j) 
            calc_inter_one_spin(l, d, i, j);
}

void calc_Q(lattice *l) {
    l->Q = 0;
    for(int i=0; i<l->N; ++i) 
        for(int j=0; j<l->N; ++j) {
            //printf("%u %u %f\n", i, j, get_inter(l, i, j));
            l->Q += get_inter(l, i, j);
        }
}

void calc_M(lattice *l) {
    l->M = 0;
    for(int i=0; i<l->N; ++i) 
        for(int j=0; j<l->N; ++j) 
            l->M += get_cos(l, i, j);
}

void calc_E(lattice *l) {
    l->E = -( l->k*l->Q + l->h*l->M );
}

void calc_all(lattice *l, const distmat *d) {
    calc_cossin(l);
    calc_inters(l, d);
    calc_Q(l);
    calc_M(l);
    calc_E(l);
}

void calc_all_one_change(lattice *l, const distmat *d, backup_lattice *bl,
                         int i, int j, float angle) {
    save_to_backup(bl, l, i, j);
    set_lat(l, i, j, angle);
    set_cos(l, i, j, cos(angle));
    set_sin(l, i, j, sin(angle));
    float Q_loc_old = 0,
          Q_loc_new = 0;
    int sii, sjj;
    for(int ii=i-l->neighb; ii<i+l->neighb+1; ii++) 
        for(int jj=j-l->neighb; jj<j+l->neighb+1; jj++) {
            periodic_conditions(ii, jj, &sii, &sjj, l->N);
            Q_loc_old += get_inter(l, sii, sjj);
            calc_inter_one_spin(l, d, sii, sjj);
            Q_loc_new += get_inter(l, sii, sjj);
        }
    l->M += get_cos(l, i, j) - bl->cos;
    l->Q += Q_loc_new - Q_loc_old;
    calc_E(l);
}

void print_lat(lattice *l, int n, int print_type,
               char *latfname, char*intfname, 
               char *ext, FILE *scaf) {
    char _latfname[32], 
         _intfname[32],
         number[32], header[128];
    memset(_latfname, 0, 32);
    memset(_intfname, 0, 32);
    memset(number, 0, 32);
    switch(print_type) {
        case LATTICE:
            sprintf(number, "%016u", n);
            sprintf(_latfname, latfname);
            /*strcat(_latfname, latfname);*/
            strcat(_latfname, number);
            strcat(_latfname, ext);
            FILE *latf = fopen(_latfname, "w");
            if (latf==NULL) {
                printf("print lat : cannot open lat file !");
                return;
            }
            fprintf(latf, "# %u %u %u %f %f %f %f %f\n", n, l->N, l->neighb, 
                    l->h, l->k, l->M, l->Q, l->E);
            for(int i=0; i<l->N; ++i) {
                fprintf(latf, "\n");
                for(int j=0; j<l->N; ++j) 
                    fprintf(latf, "\t\t%.5f", get_lat(l, i, j));
            }
            fclose(latf);
#ifdef INT
            sprintf(_intfname, intfname);
            /*strcat(_intfname, intfname);*/
            strcat(_intfname, number);
            strcat(_intfname, ext);
            FILE *intf = fopen(_intfname, "w");
            if (intf==NULL) {
                printf("print int : cannot open int file !");
                return;
            }
            fprintf(intf, "# %u %u %u %f %f %f %f %f\n", n, l->N, l->neighb, 
                    l->h, l->k, l->M, l->Q, l->E);
            for(int i=0; i<l->N; ++i) {
                fprintf(intf, "\n");
                for(int j=0; j<l->N; ++j) 
                    fprintf(intf, "\t\t%.5f", get_inter(l, i, j));
            }
            fclose(intf);
#endif
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
    lattice *l;
    distmat *d;
    backup_lattice *bl;
    int mct;
} heisenberg;

int init_heisenberg(heisenberg *he, 
        int N, float h, float k, int neighb, int a, 
        init_type it, int mct) {
    he->mct = mct;
    he->l = (lattice*) malloc(sizeof(lattice));
    he->d = (distmat*) malloc(sizeof(distmat));
    he->bl = (backup_lattice*) malloc(sizeof(backup_lattice));
    if (he->l == NULL || he->d == NULL) {
        printf("init heisenberg, could not alloc memory");
        return 1;
    }
    if (init_lattice(he->l, N, h, k, neighb, it)) return 1;
    if (init_distmat(he->d, neighb, a)) return 1; 
    if (init_backup_lattice(he->bl, neighb)) return 1; 
    calc_all(he->l, he->d);
    return 0;
}

void update(heisenberg *h) {
    float angle = rand()*PI/RAND_MAX;
    int i = floor((float) rand()/RAND_MAX*h->l->N),
        j = floor((float) rand()/RAND_MAX*h->l->N);
    float acc = (float) rand()/RAND_MAX,
          dE;

    calc_all_one_change( h->l, h->d, h->bl, i, j, angle);
    dE = h->l->E - h->bl->E;
    switch(h->mct) {
        case METROPOLIS:
            if (!(dE < 0 || acc < exp(-dE)))
                write_from_backup(h->l, h->bl);
            break;
        case GLAUBER:
            if (!(acc < 1/(1+exp(dE))))
                write_from_backup(h->l, h->bl);
            break;
    }
}


int main(int argc, char **argv) {
    if (argc < 14) {
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
         *argfile  = argv[13],
         _scafname[32];
 
    printf("#got these args : ");
    for (int i=0; i<argc; ++i) printf("%s ", argv[i]);
    printf("\n");
    printf("#N        = %s \n",  argv[1]);
    printf("#h        = %s \n",  argv[2]);
    printf("#k        = %s \n",  argv[3]);
    printf("#neighb   = %s \n",  argv[4]);
    printf("#a        = %s \n",  argv[5]);
    printf("#tmax     = %s \n",  argv[6]);
    printf("#t_lat    = %s \n",  argv[7]);
    printf("#t_sca    = %s \n",  argv[8]);
    printf("#latfname = %s \n",  argv[9]);
    printf("#intfname = %s \n",  argv[10]);
    printf("#scafname = %s \n",  argv[11]);
    printf("#ext      = %s \n",  argv[12]);
    printf("#argfile  = %s \n",  argv[13]);
    FILE * fargs = fopen(argfile, "w");
    fprintf(fargs, "#N        = %s \n",  argv[1]);
    fprintf(fargs, "#h        = %s \n",  argv[2]);
    fprintf(fargs, "#k        = %s \n",  argv[3]);
    fprintf(fargs, "#neighb   = %s \n",  argv[4]);
    fprintf(fargs, "#a        = %s \n",  argv[5]);
    fprintf(fargs, "#tmax     = %s \n",  argv[6]);
    fprintf(fargs, "#t_lat    = %s \n",  argv[7]);
    fprintf(fargs, "#t_sca    = %s \n",  argv[8]);
    fprintf(fargs, "#latfname = %s \n",  argv[9]);
    fprintf(fargs, "#intfname = %s \n",  argv[10]);
    fprintf(fargs, "#scafname = %s \n",  argv[11]);
    fprintf(fargs, "#ext      = %s \n",  argv[12]);
    fclose(fargs);
    
    int seed = 0;
    for (int i=0; i<strlen(argfile); ++i)
        seed += pow((char) argfile[i], 2);
    printf("seed = %u\n", seed);
    srand(time(NULL)*seed);
    
    memset(_scafname, 0, 32);
    strcat(_scafname, scafname);
    strcat(_scafname, ext);
    FILE *scaf = fopen(_scafname, "w"); 
    if ( scaf == NULL ) {
        printf("Cannot open scalars file ! Filename : %s >>> Error nb : %d\n", _scafname, errno);
        return 1;
    }

    int t_percent =  floor(tmax/100.);
    t_percent = (t_percent < 1 ? 1 : t_percent);
    heisenberg *he = (heisenberg*) malloc(sizeof(heisenberg));
    printf("#initialisation.");
    init_heisenberg(he, N, h, k, neighb, a, CHECKERBOARD, METROPOLIS);
    printf(".");
    print_lat(he->l, 0, LATTICE, latfname, intfname, ext, scaf);
    printf(".");
    print_lat(he->l, 0, SCALARS, latfname, intfname, ext, scaf);
    printf("done\n");
    printf("#running...\n");
    for (int t=1; t<tmax; ++t) {
        update(he);
        if (t%t_percent==0) { 
            printf("\r#%03u", (int) floor((float) t/tmax*100));
            fflush(stdout);
        }
        if (t_lat != 0 && t%t_lat==0)     print_lat(he->l, t, LATTICE, latfname, intfname, ext, scaf);
        if (t_sca != 0 && t%t_sca==0)     print_lat(he->l, t, SCALARS, latfname, intfname,  ext, scaf);
    }
    fprintf(scaf, "\n");
    fclose(scaf);
    printf("\n");

    return 0;
}

// Qstom plot
