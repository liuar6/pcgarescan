#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <stdarg.h>

typedef struct position_s{
    double b1;
    double b2;
    double b3;
    double b4;
} position_t;

#define element_t int
#define VECTOR_EXTENSION_FACTOR 0.5
#define VECTOR_MIN_CAPACITY 1

#define VEC_INIT(name, element_t)   \
    \
typedef struct vec_##name##_s{   \
    size_t size;    \
    size_t capacity;    \
    element_t *data;    \
}vec_##name##_t; \
    \
vec_##name##_t *new_##name##_vec(){   \
    return(calloc(sizeof(vec_##name##_t), 1));   \
}   \
    \
int name##_vec_add(vec_##name##_t* v, element_t element){   \
    if (v->size+1 > v->capacity){   \
        if (v->data == NULL) {  \
            v->data = malloc(VECTOR_MIN_CAPACITY * sizeof(element_t));    \
            v->capacity = VECTOR_MIN_CAPACITY;\
        }  \
        else {  \
            int extension=(int)(v->capacity) * VECTOR_EXTENSION_FACTOR; \
            if (extension < 1) extension = 1;   \
            element_t* new_data = realloc(v->data, (v->capacity + extension) * sizeof(element_t));    \
            if (new_data == NULL) return -1;    \
            v->data = new_data; \
            v->capacity = v->capacity + extension;  \
        }   \
    }   \
    v->data[v->size] = element; \
    v->size++;  \
    return 0;   \
}   \
    \
int del_##name##_vec(vec_##name##_t* v){  \
    free(v->data);  \
    free(v);    \
    return 0; \
}   \

VEC_INIT(position, position_t)
VEC_INIT(motif, vec_position_t*)
#define motif_t vec_position_t
#define new_motif new_position_vec
#define motif_add position_vec_add
#define del_motif del_position_vec

vec_motif_t* read_motifs(const char* file){
    FILE* f=fopen(file, "r");
    char buffer[4096];
    position_t position;
    motif_t *motif=NULL;
    vec_motif_t *motifs=new_motif_vec();
    while(fgets(buffer, 4096, f)!=NULL) {
        if (buffer[0] == '>') {
            motif = new_motif();
            motif_vec_add(motifs, motif);
        } else if (strlen(buffer) == 0) {
            continue;
        } else {
            char* s = buffer;
            position.b1=strtod(s, &s);
            position.b2=strtod(s, &s);
            position.b3=strtod(s, &s);
            position.b4=strtod(s, &s);
            motif_add(motif, position);
        }
    }
    fclose(f);
    return motifs;
}
void print_motifs(vec_motif_t* motifs){
    motif_t * motif;
    for (int i=0; i < motifs->size; ++i) {
        motif=motifs->data[i];
        printf(">%d\n", i+1);
        for (int j=0; j<motif->size; j++) {
            position_t *pos=&(motif->data[j]);
            printf("%f\t%f\t%f\t%f\n", pos->b1, pos->b2, pos->b3, pos->b4);
        }
    }
}
void free_motifs(vec_motif_t* motifs){
    for (int i=0; i < motifs->size; ++i) del_motif(motifs->data[i]);
    del_motif_vec(motifs);
}
#define ent(a) (a>0?a*log(a)/0.69314718055995:0)
double pos_ent(position_t *pos){
    double total=pos->b1+pos->b2+pos->b3+pos->b4;
    double ret=2+ent(pos->b1/total)+ent(pos->b2/total)+ent(pos->b3/total)+ent(pos->b4/total);
    return ret;

}
double pos_score(void *_a, void *_b){
    position_t *a =  _a;
    position_t *b =  _b;
    double scale_factor=(pos_ent(a)+pos_ent(b))/4;
    //printf("%f\n", scale_factor);
    double sum_a=a->b1+a->b2+a->b3+a->b4;
    double sum_b=b->b1+b->b2+b->b3+b->b4;
    double raw_score=1-fabs(a->b1/sum_a-b->b1/sum_b)-fabs(a->b2/sum_a-b->b2/sum_b)-fabs(a->b3/sum_a-b->b3/sum_b)-fabs(a->b4/sum_a-b->b4/sum_b);
    return raw_score*scale_factor;
}
double pos_score1(void *_a, void *_b){
    position_t *a =  _a;
    position_t *b =  _b;
    double enta=pos_ent(a);
    double entb=pos_ent(b);
    double sum_a=a->b1+a->b2+a->b3+a->b4;
    double sum_b=b->b1+b->b2+b->b3+b->b4;
    double scale_factor=((enta>entb)?entb:enta)/2;
    double raw_score=1-fabs(a->b1/sum_a-b->b1/sum_b)-fabs(a->b2/sum_a-b->b2/sum_b)-fabs(a->b3/sum_a-b->b3/sum_b)-fabs(a->b4/sum_a-b->b4/sum_b);
    return raw_score*scale_factor;
}

typedef struct cmp_motif_args_s{
    double (*score_func)(void *a, void *b);
    double edge_penalty;
    double gap_penalty;
    int is_local;
} cmp_motif_args_t;

double cmp_motif(motif_t *a, motif_t *b,void *args){
    double score[101][101]; /* max motif length is 100! */
    double sum_score[101][101];
    double (*score_func)(void *a, void *b) = ((cmp_motif_args_t*) args)->score_func;
    double edge_penalty = ((cmp_motif_args_t*) args)->edge_penalty;
    double gap_penalty = ((cmp_motif_args_t*) args)->gap_penalty;
    int is_local = ((cmp_motif_args_t*) args)->is_local;
    /* initiate score matrix */
    for (int i=0; i < a->size; ++i) for (int j=0; j < b->size; ++j) score[i][j]=score_func(&(a->data[i]), &(b->data[j]));
    /* initiate main matrix */
    sum_score[0][0] = 0;
    for (int i=1; i < a->size+1; ++i) sum_score[i][0]=sum_score[i-1][0]-(is_local)?0:edge_penalty;
    for (int j=1; j < b->size+1; ++j) sum_score[0][j]=sum_score[0][j-1]-(is_local)?0:edge_penalty;
    /* start filling the matrix */
    for (int i=1; i < a->size+1; ++i) {
        for (int j=1; j < b->size+1; ++j) {
            double score1 = sum_score[i-1][j-1] + score[i-1][j-1];
            /* if i == a->size, all positions in motif a are consumed, thus use edge penalty */
            double score2 = sum_score[i][j-1]-((i == a->size)?edge_penalty:gap_penalty);
            /* if j == b->size, all positions in motif a are consumed, thus use edge penalty */
            double score3 = sum_score[i-1][j]-((j == b->size)?edge_penalty:gap_penalty);
            sum_score[i][j] = (score1>score2)?((score1>score3)?score1:score3):((score2>score3)?score2:score3);
            if (is_local && sum_score[i][j] < 0) sum_score[i][j] = 0;
        }
    }
    /* print the matrice for debugging.
    for (int j = 0; j < b->size + 1; ++j) {
        for (int i = 0; i < a->size + 1; ++i) printf("%f ", sum_score[i][j]);
        printf("\n");
    }
    printf("\n");
    for (int j = 0; j < b->size; ++j) {
        for (int i = 0; i < a->size; ++i) printf("%f ", score[i][j]);
        printf("\n");
    }
     */
    if (is_local) {
        double max_score = 0;
        for (int i = 0; i < a->size + 1; ++i)
            for (int j = 1; j < b->size + 1; ++j)
                if (sum_score[i][j] > max_score) max_score = sum_score[i][j];
        return max_score;
    } else {
        return sum_score[a->size][b->size];
    }
}

void cmp_motif_usage(const char *msg, ...){
    va_list argptr;
    va_start(argptr, msg);
    vfprintf(stderr, msg, argptr);
    va_end(argptr);
    fprintf(stderr, "%s", "\n\
Usage:  cmpmotif [options] -a <motif file> -b <motif file>\n\
 [options]\n\
-s/--score          : score function for calculating pairwise base similarity \n\
-g/--gap            : gap penalty, default: 2.\n\
-e/--edge           : gap penalty for two ends, default: 0.2.\n\
-l/--local          : perform local alignments instead of global alignments\n\
-h/--help           : show help informations\n\
");
}

int main(int argc, char *argv[]) {
    if (argc == 1) {cmp_motif_usage("Invalid arguments!\n"); exit(0);}
    const char *shortOptions = "hvla:b:s:g:e:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "score" , required_argument , NULL, 's' },
                    { "gap" , required_argument , NULL, 'g' },
                    { "edge" , required_argument , NULL, 'e' },
                    { "local" , no_argument , NULL, 'l' },
                    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
            };
    char c;
    cmp_motif_args_t cmp_args;
    const char *input_a=NULL;
    const char *input_b=NULL;
    cmp_args.score_func=&pos_score;
    cmp_args.edge_penalty=0.2;
    cmp_args.gap_penalty=2;
    cmp_args.is_local=0;

    while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
        switch(c){
            case 'h':
                cmp_motif_usage("");
                break;
            case 'v':
                break;
            case 'a':
                input_a=optarg;
                break;
            case 'b':
                input_b=optarg;
                break;
            case 's':
                switch(optarg[0]){
                    case '0':
                        cmp_args.score_func=pos_score;
                        break;
                    case '1':
                        cmp_args.score_func=pos_score1;
                        break;
                    default:
                        cmp_motif_usage("Unrecognizable score function: '%c'.\n", optarg[0]);
                }
                break;

            case 'g':
                cmp_args.gap_penalty=strtod(optarg, NULL);
                break;
            case 'e':
                cmp_args.edge_penalty=strtod(optarg, NULL);
                break;
            case 'l':
                cmp_args.is_local=1;
                break;
            default:
                cmp_motif_usage("Invalid arguments!\n");
                exit(0);
        }
    if (argc != optind) {
        cmp_motif_usage("Invalid arguments!\n");
        exit(0);
    }
    vec_motif_t *motifs_a, *motifs_b;
    if (input_a == NULL) cmp_motif_usage("The input file is required!");
    motifs_a=read_motifs(input_a);
    if (input_b == NULL)  motifs_b=motifs_a;
    else motifs_b=read_motifs(input_b);

    for (int i=0; i < motifs_a->size; ++i){
        for (int j=0; j < motifs_b->size; ++j){
            double score = cmp_motif(motifs_a->data[i], motifs_b->data[j], &cmp_args);
            printf("%.2f ", score);
        }
        printf("\n");
    }
    free_motifs(motifs_a);
    if (input_b != NULL) free_motifs(motifs_b);
}

