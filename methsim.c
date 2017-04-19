#include "utils.h"
#include "number.h"
#include "stats.h"
#include "sequence.h"
#include "htslib/hts.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

struct meth_node {
    char strand;
    uint32_t start; // 0 based
    uint32_t m_rate;
};

struct meth_chrom {
    int m, n;
    struct meth_node *a;
};

KHASH_MAP_INIT_STR(chr, struct meth_chrom *)

typedef kh_chr_t reghash_t;

typedef unsigned int mut_t;

#define MUT_SHIFT 2
static mut_t mutmsk = 3;
static mut_t mut_type_none = 0;
static mut_t mut_type_subs = 1;
static mut_t mut_type_del  = 2;
static mut_t mut_type_ins  = 3;

struct mutseq {
    int l, m;
    mut_t *s;
};

struct args {
    // const char *bed_fname;

    // reference sequence in fasta format
    const char *fasta_fname;

    // hash struct for caching methylated data, key is chromosome, value point to an dynamic
    // array of methylated bed record
    reghash_t *hash;

    // file handler point to simulated file
    gzFile read1_fp;
    gzFile read2_fp;
    gzFile meth1_fp;
    gzFile meth2_fp;

    // variants reprot, stdout for default
    FILE *report_fp;

    // parameters for simulation
    double err_rate;
    double mut_rate;
    double indel_frac;
    double indel_extend;

    // read length for paired reads
    int read1_length;    
    int read2_length;

    // predefined fragment size and standard deviation, assume fragment size normal distributed
    int insert_size;    
    int std_dev;

    // read pairs
    uint64_t n_pairs;
    
    // max ratio of N bases, parameter of wgsim, here I delete it, only generate nonN bases
    // double max_ratio;

    // assume haplotype mode in default, so I delete this parameter now
    int haplotype_mode;
} args = {
    // .bed_fname = NULL,
    .fasta_fname = NULL,
    .hash = NULL,
    .read1_fp = NULL,
    .read2_fp = NULL,
    .meth1_fp = NULL,
    .meth2_fp = NULL,
    .report_fp = NULL,
    .err_rate = 0.02,
    .mut_rate = 0.001,
    .indel_frac = 0.15,
    .indel_extend = 0.3,
    // .max_n_ratio = 0,
    .read1_length = 100,
    .read2_length = 100,
    .insert_size = 500,
    .std_dev = 50,
    .n_pairs = 1000000,
    .haplotype_mode = 0,
};

int usage()
{
    fprintf(stderr, "methsim : Short methylated reads simulator.\n");
    fprintf(stderr, "The program adapt from Liheng's whole genome short simulator wgsim.\n");
    fprintf(stderr, "Visit wgsim -> https://github.com/lh3/wgsim\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "methsim [options] <in.ref.fasta> <out.reads1.fq> <out.reads2.fq> <out.meth_reads1.fq> <out.meth_reads2.fq>\n");
    fprintf(stderr, "Options:  -bed         mandatory methylated bed file\n");
    fprintf(stderr, "          -report      variation report [stdout]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options same with liheng's wgsim:\n");
    fprintf(stderr, "         -e FLOAT      base error rate [%.3f]\n", args.err_rate);
    fprintf(stderr, "         -d INT        outer distance between the two ends [500]\n");
    fprintf(stderr, "         -s INT        standard deviation [50]\n");
    fprintf(stderr, "         -N INT        number of read pairs [1000000]\n");
    fprintf(stderr, "         -1 INT        length of the first read [70]\n");
    fprintf(stderr, "         -2 INT        length of the second read [70]\n");
    fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", args.mut_rate);
    fprintf(stderr, "         -R FLOAT      fraction of indels [%.2f]\n", args.indel_frac);
    fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", args.indel_extend);
    fprintf(stderr, "         -S INT        seed for random generator [-1]\n");
//    fprintf(stderr, "         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", args.max_n_ratio);
    fprintf(stderr, "         -h            haplotype mode\n");
    fprintf(stderr, "         -help         see this message\n");
    fprintf(stderr, "Homepage:\n");
    fprintf(stderr, "   https://github.com/shiquan/meth-sim\n");
    return 1;
}

// Accept format:
//  chr1    3111719 3111721 '1/9'   111     -
reghash_t *load_meth_bed(const char *fname) {
    reghash_t *hash = kh_init(chr);
    kstring_t str = {0, 0, 0};
    htsFile *fp = hts_open(fname, "r");
    for ( ;; ) {
        if ( hts_getline(fp, KS_SEP_LINE, &str) < 0 )
            break;
        if ( str.l == 0 )
            continue;
        if ( str.s[0] == '#' || str.s[0] == '/')
            continue;

        int *splits;
        int nfields;
        khint_t k;
        splits = ksplit(&str, 0, &nfields);

        char *name = str.s + splits[0];
        int start = str2int(str.s + splits[1]); 
        
        if ( nfields != 6 ) {
            warnings("malformed format. %s\t%d", name, start);
            free(splits);
            continue;
        }

        k = kh_get(chr, hash, name);
        if ( k == kh_end(hash) ) {
            int ret;
            char *str = strdup(name);
            k = kh_put(chr, hash, str, &ret);
            struct meth_chrom *chr = malloc(sizeof(struct meth_chrom));
            chr->m = chr->n = 0;
            chr->a = NULL;
            kh_val(hash, k) = chr;
        }

        struct meth_chrom *chr = kh_val(hash, k);
        
        if ( chr->n == chr->m ) {
            chr->m = chr->m == 0 ? 10 : chr->m << 1;
            chr->a = realloc(chr->a, chr->m *sizeof(struct meth_node));
        }

        struct meth_node *node = &chr->a[chr->n];
        node->start = start;
        node->m_rate = str2int(str.s + splits[4]);
        node->strand = *(str.s + splits[5]);
        if ( node->strand == '-' )
            node->start ++;
        
        chr->n++;
        str.l = 0;
    }
    free(str.s);
    hts_close(fp);
    return hash;
}

int cmp_func(const void *a, const void *b)
{
    struct meth_node *a1 = (struct meth_node*)a;
    struct meth_node *b1 = (struct meth_node*)b;
    return a1->start - b1->start;
}

void sort_meth_hash(reghash_t *hash)
{
    khint_t k;
    for ( k = kh_begin(hash); k != kh_end(hash); ++k ) {
        if ( !kh_exist(hash, k) )
            continue;

        struct meth_chrom *chr = kh_val(hash, k);
        qsort(chr->a, chr->n, sizeof(struct meth_node), cmp_func);
        // debug_print("%s\t%d", kh_key(hash, k), chr->n);
    }
}
int generate_vars(const kseq_t *ks, struct mutseq *hap1, struct mutseq *hap2)
{
    khint_t k;
    struct meth_chrom *mchr = NULL;
    k = kh_get(chr, args.hash, ks->name.s);
    
    if ( kh_exist(args.hash, k)) {
        mchr = kh_val(args.hash, k);
    } else {
        warnings("chromosome %s not found in methylated data.", ks->name.s);
    }

    struct mutseq *ret[2];
    int i, deleting = 0, iter = 0;
    ret[0] = hap1;
    ret[1] = hap2;
    ret[0]->l = ks->seq.l;
    ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m;
    ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t*)calloc(ret[0]->m, sizeof(mut_t));
    ret[1]->s = (mut_t*)calloc(ret[0]->m, sizeof(mut_t));
    
    for ( i = 0; mchr != NULL && i < mchr->n; ++i ) {
        struct meth_node *node = &mchr->a[i];
        if ( node->strand == '+' ) {
            if (seq2code4(ks->seq.s[node->start]) != 1 ) {
                error("The base in the methylated position %s:%d is not a C vs %c",
                      ks->name.s, node->start+1, ks->seq.s[node->start]);
            }
            ret[0]->s[node->start] = node->m_rate << MUT_SHIFT;
        } else {
            if ( seq2code4(ks->seq.s[node->start]) != 2) {
                error("The base in the methylated position %s:%d is not a C vs %c",
                      ks->name.s, node->start+1, ks->seq.s[node->start]);
            }
            ret[1]->s[node->start] = node->m_rate << MUT_SHIFT;
        }           
    }

    for ( i = 0; i < ks->seq.l; ++i ) {
        // Skip if methylated.
        if ( (((ret[0]->s[i] & mutmsk) == 0 ) && (ret[0]->s[i] >> MUT_SHIFT)) ||
             (((ret[1]->s[i] & mutmsk) == 0 ) && (ret[1]->s[i] >> MUT_SHIFT)) )
            continue;
        
        int c = seq2code4(ks->seq.s[i]);

        // delete
        if ( deleting ) {
            if ( drand48() < args.indel_extend ) {
                if ( deleting & 1 )
                    ret[0]->s[i] = mut_type_del;
                if ( deleting & 2 )
                    ret[1]->s[i] = mut_type_del;
                continue;
            } else {
                deleting = 0;
            }
        }

        // variantion
        // for variations methylate bits will be filled by alts, so check the vartype before check methylation ratio
        if ( c < 4 && drand48() < args.mut_rate ) {
            // substitution
            if ( drand48() >= args.indel_frac ) {
                int c;
                c = ((c + (int)(drand48()*3.0 + 1)) & mutmsk) << MUT_SHIFT;
                if ( drand48() < 0.333333 ) {
                    ret[0]->s[i] = ret[1]->s[i] = c | mut_type_subs;
                } else {
                    ret[drand48() < 0.5 ? 0 : 1]->s[i] = c | mut_type_subs;
                }
            } else { // indel
                if ( drand48() < 0.5 ) { // deletion
                    if ( drand48() < 0.333333) { // hom-del
                        ret[0]->s[i] = ret[1]->s[i] = mut_type_del;
                        deleting = 3;
                    } else {
                        deleting = drand48() < 0.5 ? 1 : 2;
                        ret[deleting-1]->s[i] = mut_type_del;                        
                    }
                } else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        num_ins++;
                        ins = (ins<<2)|(int)(drand48()*4.0);
                    } while (num_ins < 4 && drand48() < args.indel_extend);
                                        
                    if ( drand48() < 0.333333 ) { // hom-ins
                        ret[0]->s[i] = ret[1]->s[i] = (num_ins<<10) | (ins<<2) | mut_type_ins;
                    } else { // het-ins
                        ret[drand48() < 0.5 ? 0 : 1]->s[i] = (num_ins<<10) | (ins<<2) | mut_type_ins;
                    }                    
                }
            } //end indel
        } // end variantion
    }
    
    // print mutref report
    for ( i = 0; i != ks->seq.l; ++i ) {
        int c[3];
        c[0] = seq2code4(ks->seq.s[i]);
        if ( c[0] >= 4 )
            continue;
        
        c[1] = hap1->s[i];
        c[2] = hap2->s[i];
        mut_t type1 = c[1] & mutmsk;
        mut_t type2 = c[2] & mutmsk;
        int bits1 = c[1] >> MUT_SHIFT;
        int bits2 = c[2] >> MUT_SHIFT;
        int d2 = 0; // d2 keeps the last end of deletion
        if ( type1 || type2 ) {
            if (type1 == type2) { // hom
                if ( type1 == mut_type_subs ) {
                    fprintf(args.report_fp, "%s\t%d\t%c\t%c\t-\n", ks->name.s, i+1, "ACGTN"[c[0]], "ACGTN"[bits2&mutmsk]);
                } else if ( type1 == mut_type_del ) {
                    if ( i >= d2 ) {
                        fprintf(args.report_fp, "%s\t%d\t", ks->name.s, i+1);
                        for ( d2=i; d2 < ks->seq.l && hap1->s[d2] == hap2->s[d2] && (hap1->s[d2]&mutmsk) == mut_type_del; ++d2)
                            putc("ACGTN"[seq2code4(ks->seq.s[d2])], args.report_fp);
                        fprintf(args.report_fp, "\t-\t-\n");
                    }
                } else if ( type1 == mut_type_ins ) {
                    fprintf(args.report_fp, "%s\t%d\t-\t", ks->name.s, i+1);
                    int n = bits1>>8;
                    while ( n > 0 ) {
                        fputc("ACGTN"[bits1&mutmsk], args.report_fp);
                        bits1>>=2;
                        n--;
                    }
                    fprintf(args.report_fp, "\t-\n");
                }
            } else { // het
                if ( type1 == mut_type_subs || type2 == mut_type_subs) {
                    fprintf(args.report_fp, "%s\t%d\t%c\t%c\t-\n", ks->name.s, i+1, "ACGTN"[c[0]],  "XACMGRSVTWYHKDBN"[1<<(bits1&mutmsk)|1<<(bits2&mutmsk)]);
                } else if ( type1 == mut_type_del ) {
                    if ( i >= d2 ) {
                        fprintf(args.report_fp, "%s\t%d\t", ks->name.s, i+1);
                        for ( d2=i; d2 < ks->seq.l && hap1->s[d2] != hap2->s[d2] && (hap1->s[d2]&mutmsk) == mut_type_del; ++d2)
                            putc("ACGTN"[seq2code4(ks->seq.s[d2])], args.report_fp);
                        fprintf(args.report_fp, "\t-\t-\n");
                    }
                } else if ( type2 == mut_type_del) {
                    if ( i >= d2 ) {
                        fprintf(args.report_fp, "%s\t%d\t", ks->name.s, i+1);
                        for ( d2=i; d2 < ks->seq.l && hap1->s[d2] != hap2->s[d2] && (hap2->s[d2]&mutmsk) == mut_type_del; ++d2)
                            putc("ACGTN"[seq2code4(ks->seq.s[d2])], args.report_fp);
                        fprintf(args.report_fp, "\t-\t-\n");
                    }                    
                } else if ( type1 == mut_type_ins ) {
                    fprintf(args.report_fp, "%s\t%d\t-\t", ks->name.s, i+1);
                    int n = bits1>>8;
                    while ( n > 0 ) {
                        fputc("ACGTN"[bits1&mutmsk], args.report_fp);
                        bits1>>=2;
                        n--;
                    }
                    fprintf(args.report_fp, "\t-\n");
                } else if ( type2 == mut_type_ins ) {
                    fprintf(args.report_fp, "%s\t%d\t-\t", ks->name.s, i+1);
                    int n = bits2>>8;
                    while ( n > 0 ) {
                        fputc("ACGTN"[bits2&mutmsk], args.report_fp);
                        bits2>>=2;
                        n--;
                    }
                    fprintf(args.report_fp, "\t-\n");
                }
            }
        }
    }
    return 0;
}

void print_seqs(kseq_t *ks, int length, int n_pairs, struct mutseq *hap1, struct mutseq *hap2)
{
    int max_size, Q, ii;
    max_size = args.read1_length > args.read2_length ? args.read1_length : args.read2_length;
    Q = ( args.err_rate == 0.0 ) ? 'I' : (int)(-10.0 *log(args.err_rate)/log(10.0) + 0.499) + 33;

    int s[2]; 
    gzFile fp[2], mp[2];
    uint64_t n_sub[2] = {0, 0}, n_indel[2] = {0, 0}, n_err[2] = {0, 0}, ext_coor[2] = {0,0};
    uint8_t *temp_seq[2], *temp_ms[2];
    temp_seq[0] = (uint8_t*)calloc(max_size+2, 1);
    temp_seq[1] = (uint8_t*)calloc(max_size+2, 1);
    temp_ms[0]  = (uint8_t*)calloc(max_size+2, 1);
    temp_ms[1]  = (uint8_t*)calloc(max_size+2, 1);
    kstring_t temp1 = { 0, 0, 0};
    kstring_t temp2 = { 0, 0, 0};
    for ( ii = 0; ii != n_pairs; ++ii ) {
        // double ran;
        int pos;
        int dist;
        int is_flip = 0;
        int k, i, j;

      regenerate:
        // generate random start position
        do {
            dist = ran_normal((double)args.insert_size, (double)args.std_dev);
            dist = dist < max_size ? max_size : dist;
            pos = (int)((length-dist+1) *drand48());
            // emit all N's
            for ( i = 0; i < max_size; ++i ) { 
                if ( seq2code4(ks->seq.s[pos+i]) >= 4 ) {
                    pos = -1;
                    break;                    
                }                    
            }                
        } while (pos < 0 || pos >= length || pos+dist-1 >= length);

        // debug_print("pos %u ii %d", pos, ii);
        // flip or not
        if ( drand48() < 0.5 ) {
            fp[0] = args.read1_fp;
            fp[1] = args.read2_fp;
            mp[0] = args.meth1_fp;
            mp[1] = args.meth2_fp;
            s[0] = args.read1_length;
            s[1] = args.read2_length;
        } else {
            fp[0] = args.read2_fp;
            fp[1] = args.read1_fp;
            mp[0] = args.meth2_fp;
            mp[1] = args.meth1_fp;
            s[0] = args.read2_length;
            s[1] = args.read1_length;            
            is_flip = 1;
        }        
/*
  \ 3'  --> read1               <--  read2
   \---------TNNNNNNNNNNNNNNNNNNNNA---------
    ---------ANNNNNNNNNNNNNNNNNNNNT---------\
   5'   read2-->                  <-- read1  3'
 */

        // forward strand or reverse
        int is_strand = drand48() < 0.5 ? 1 : 0;
        if ( is_strand ) {  // forward strand

#define BRANCH(_pos, x) do {                                            \
                for (ext_coor[x] = -10, i = (_pos), k = 0; i < length && k < max_size; ) { \
                    int c = seq2code4(ks->seq.s[i]);                    \
                    if ( c>= 4 )                                        \
                        goto regenerate;                                \
                    mut_t type = hap1->s[i] & mutmsk;                   \
                    int bits = hap1->s[i] >> MUT_SHIFT;                    \
                    if ( ext_coor[x] < 0 ) {\
                        if ( type != mut_type_none && type != mut_type_subs ) \
                            continue;\
                        ext_coor[x] = i;\
                    }\
                    if ( type == mut_type_none ) {                      \
                        temp_seq[x][k] = c;                             \
                        if ( c == 1 && drand48() < (double)bits/1000.0 ) {\
                            temp_ms[x][k] = 1;                          \
                        } else {                                        \
                            temp_ms[x][k] = c == 1 ? 3 : c;             \
                        }                                               \
                        i++;                                            \
                        k++;                                            \
                    } else if ( type == mut_type_subs ) {               \
                        int mut = bits & mutmsk;                             \
                        temp_seq[x][k] = mut;                           \
                        temp_ms[x][k] = mut == 1 ? 3 : mut;          \
                        i++;                                            \
                        k++;                                            \
                        n_sub[x]++;                                     \
                    } else if ( type == mut_type_del) {                 \
                        i++;                                            \
                        n_indel[x]++;                                   \
                    } else if ( type == mut_type_ins) {                 \
                        int ins_num = bits>>8;                          \
                        for ( j = 0; j < ins_num; ++j ) {               \
                            int mut = (bits>>(j*2))&mutmsk;             \
                            temp_seq[x][k] = mut;                       \
                            temp_ms[x][k] = mut == 1 ? 3 : mut;         \
                            k++;                                        \
                        }                                               \
                        n_indel[x]++;                                   \
                    }                                                   \
                }                                                       \
            } while(0)
            // read 1
            BRANCH(pos, 0);
            // read 2
            BRANCH(pos+dist-args.read2_length+1, 1);
            
#undef BRANCH
            
        } else { // reverse strand
#define BRANCH(_pos, x) do {                                            \
                for (ext_coor[x] = -10, i = (_pos), k = 0; i < length && k < max_size; ) { \
                    int c = 3-seq2code4(ks->seq.s[i]);                  \
                    if ( c>= 4 )                                        \
                        goto regenerate;                                \
                    mut_t type = hap2->s[i] & mutmsk;                   \
                    int bits = hap2->s[i] >> MUT_SHIFT;                    \
                                        if ( ext_coor[x] < 0 ) {\
                        if ( type != mut_type_none && type != mut_type_subs ) \
                            continue;\
                        ext_coor[x] = i;\
                    }\
                    if ( type == mut_type_none ) {                      \
                        temp_seq[x][k] = c;                             \
                        if ( c == 1 && drand48() < (double)bits/1000.0 ) {\
                            temp_ms[x][k] = 1;                          \
                        } else {                                        \
                            temp_ms[x][k] = c == 1 ? 3 : c;             \
                        }                                               \
                        i++;                                            \
                        k++;                                            \
                    } else if ( type == mut_type_subs ) {               \
                        int mut = bits & mutmsk;                             \
                        temp_seq[x][k] = mut;                           \
                        temp_ms[x][k] = mut == 1 ? 3 : mut;             \
                        i++;                                            \
                        k++;                                            \
                        n_sub[x]++;                                     \
                    } else if ( type == mut_type_del) {                 \
                        i++;                                            \
                        n_indel[x]++;                                   \
                    } else if ( type == mut_type_ins) {                 \
                        int ins_num = bits>>8;                          \
                        for ( j = 0; j < ins_num; ++j ) {               \
                            int mut = (bits>>(j*2))&mutmsk;             \
                            temp_seq[x][k] = mut;                       \
                            temp_ms[x][k] = mut == 1 ? 3 : mut;         \
                            k++;                                        \
                        }                                               \
                        n_indel[x]++;                                   \
                    }                                                   \
                }                                                       \
            } while(0)
            // read 1
            BRANCH(pos, 0);
            // read 2
            BRANCH(pos+dist-args.read2_length+1, 1);

#undef BRANCH
        }

        for ( k = 0; k < s[1]; ++k ) {
            temp_seq[1][k] = 3 - temp_seq[1][k];
            temp_ms[1][k] = 3 - temp_ms[1][k];
        }
        
        for ( j = 0; j < 2; ++j ) {
            
            // header
            ksprintf(&temp1, "@%s_%llu_%llu_%llu:%llu:%llu_%llu:%llu:%llu_%llx/%d\n", ks->name.s, ext_coor[0]+1, ext_coor[1]+1,
                     n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
                     (long long)ii, j==0? is_flip+1 : 2-is_flip);
            ksprintf(&temp2, "@%s_%llu_%llu_%llu:%llu:%llu_%llu:%llu:%llu_%llx/%d\n", ks->name.s, ext_coor[0]+1, ext_coor[1]+1,
                     n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
                     (long long)ii, j==0? is_flip+1 : 2-is_flip);
            
            for (i = 0; i < s[j]; ++i) {
                kputc("ACGTN"[(int)temp_seq[j][i]], &temp1);
                kputc("ACGTN"[(int)temp_ms[j][i]], &temp2);
            }
            kputs("\n+\n", &temp1);
            kputs("\n+\n", &temp2);
            
            for (i = 0; i < s[j]; ++i) {
                kputc(Q, &temp1);
                kputc(Q, &temp2);
            }
            kputc('\n', &temp1);
            kputc('\n', &temp2);
            gzprintf(fp[j], temp1.s);
            gzprintf(mp[j], temp2.s);
            temp1.l = 0;
            temp2.l = 0;
        }
    }
    free(temp1.s);
    free(temp2.s);        
}

void generate_simulate_data()
{
    uint64_t total_length = 0;
    int l = 0;
    int n_ref = 0;    
    gzFile fp;
    kseq_t *ks;
    
    fp = gzopen(args.fasta_fname, "r");
    if ( fp == NULL )
        error("%s : %s.", args.fasta_fname, strerror(errno));
    
    ks = kseq_init(fp);
    while ( (l = kseq_read(ks)) >= 0 ) {
        total_length += l;
        ++n_ref;
    }

    LOG_print("%d sequences, total length : %llu", n_ref, total_length);

    kseq_destroy(ks);
    gzclose(fp);

    fp = gzopen(args.fasta_fname, "r");
    ks = kseq_init(fp);
    
    while ( (l = kseq_read(ks)) >= 0 ) {
        uint64_t n_pairs = (uint64_t)((long double) l / total_length * args.n_pairs + 0.5);

        if ( l < args.insert_size + 3 * args.std_dev ) {
            warnings("skip sequence '%s' as it is shorter than %d!\n", ks->name.s, args.insert_size * 3 * args.std_dev);
            continue;
        }
        LOG_print("Generate reads for sequence %s, read pairs : %lld", ks->name.s, n_pairs);
        struct mutseq rseq[2] = { {0, 0, 0}, {0, 0, 0}};
        // generate variants and methylation positions and print them out
        generate_vars(ks, rseq, rseq+1);
        LOG_print("Generate variants success.");
        // generate output files
        print_seqs(ks, l, n_pairs, rseq, rseq+1);
        free(rseq[0].s);
        free(rseq[1].s);
    }
    
}
int release_memory()
{
    khint_t k;
    for ( k = kh_begin(args.hash); k != kh_end(args.hash); ++k ) {
        if ( kh_exist(args.hash, k ) ) {
            struct meth_chrom *chr = kh_val(args.hash, k);
            free(chr->a);
            free(chr);
            kh_del(chr, args.hash, k);
        }
    }
    kh_destroy(chr, args.hash);
    args.hash = NULL;

    gzclose(args.read1_fp);
    gzclose(args.read2_fp);
    gzclose(args.meth1_fp);
    gzclose(args.meth2_fp);
    fclose(args.report_fp);
    return 0;
}

int parse_args(int ac, char **av)
{
    if ( ac < 6 )
        return usage();

    // const char *fasta_fname = NULL;
    const char *bed_fname = NULL;
    const char *reads1_fname = NULL;
    const char *reads2_fname = NULL;
    const char *meth1_fname = NULL;
    const char *meth2_fname = NULL;

    const char *report_fname = NULL;
    
    const char *err_rate = NULL;
    const char *dist = NULL;
    const char *dev = NULL;
    const char *pairs = NULL;
    const char *length_read1 = NULL;
    const char *length_read2 = NULL;
    const char *mut_rate = NULL;
    const char *indel_frac = NULL;
    const char *indel_extend = NULL;
    const char *seed_str = NULL;
    // const char *max_ratio = NULL;
    
    int i;
    for ( i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-bed") == 0 && bed_fname == NULL )            
            var = &bed_fname;
        else if ( strcmp(a, "-report") == 0 )
            var = &report_fname;
        else if ( strcmp(a, "-e") == 0 )
            var = &err_rate;
        else if ( strcmp(a, "-d") == 0 )
            var = &dist;
        else if ( strcmp(a, "-s") == 0 )
            var = &dev;
        else if ( strcmp(a, "-N") == 0 )
            var = &pairs;
        else if ( strcmp(a, "-1") == 0 )
            var = &length_read1;
        else if ( strcmp(a, "-2") == 0 )
            var = &length_read2;
        else if ( strcmp(a, "-r") == 0 )
            var = &mut_rate;
        else if ( strcmp(a, "-R") == 0 )
            var = &indel_frac;
        else if ( strcmp(a, "-X") == 0 )
            var = &indel_extend;
        else if ( strcmp(a, "-S") == 0 )
            var = &seed_str;
        /* else if ( strcmp(a, "-A") == 0 ) */
        /*     var = &max_ratio; */
        
        if ( var != 0 ) {
            if ( i == ac ) 
                error("Miss an argument after %s", a);
            *var = av[i++];
            continue;
        }
            
        if ( strcmp(a, "-h") == 0 ) {
            args.haplotype_mode = 1;
            continue;
        }

        if ( strcmp(a, "-help") == 0 ) {
            return usage();
        }
        
        if ( args.fasta_fname == NULL ) {
            args.fasta_fname = a;
            continue;
        }

        if ( reads1_fname == NULL ) {
            reads1_fname = a;
            continue;
        }

        if ( reads2_fname == NULL ) {
            reads2_fname = a;
            continue;
        }

        if ( meth1_fname == NULL ) {
            meth1_fname = a;
            continue;
        }
        
        if ( meth2_fname == NULL ) {
            meth2_fname = a;
            continue;
        }
        
        error("Unknown argument : %s.", a);
    }

    if ( args.fasta_fname == NULL ) 
        error("No reference file specified.");

    if ( bed_fname == NULL )
        error("No methylated data specified. bed format file required");
    
    if ( reads1_fname == NULL || reads2_fname == NULL || meth1_fname == NULL || meth2_fname == NULL )
        return usage();

    if ( err_rate )
        args.err_rate = str2float((char*)err_rate);

    if ( dist )
        args.insert_size = str2int((char*)dist);

    if ( dev )
        args.std_dev = str2int((char*)dev);

    if ( pairs )
        args.n_pairs = str2int((char*)pairs);

    if ( length_read1 )
        args.read1_length = str2int((char*)length_read1);

    if ( length_read2 )
        args.read2_length = str2int((char*)length_read2);

    if ( mut_rate )
        args.mut_rate = str2float((char*)mut_rate);

    if ( indel_frac )
        args.indel_frac = str2float((char*)indel_frac);

    if ( indel_extend )
        args.indel_extend = str2float((char*)indel_extend);

    /* if ( max_ratio ) */
    /*     args.max_ratio = str2float((char*)max_ratio); */

    int seed = -1;
    if ( seed_str ) {
        seed = str2int((char*)seed_str);
        update_seed(seed);
    }

    seed = get_seed();
    srand48(seed);

    args.report_fp = report_fname == NULL ? stdout : fopen(report_fname, "w");
    if ( args.report_fp == NULL ) {
        warnings("%s : %s", report_fname, strerror(errno));
        args.report_fp = stdout;
    }

    args.read1_fp = gzopen(reads1_fname, "w");
    if ( args.read1_fp == NULL)
        error("%s : %s", reads1_fname, strerror(errno));

    args.read2_fp = gzopen(reads2_fname, "w");
    if ( args.read2_fp == NULL)
        error("%s : %s", reads2_fname, strerror(errno));

    args.meth1_fp = gzopen(meth1_fname, "w");
    if ( args.meth1_fp == NULL)
        error("%s : %s", meth1_fname, strerror(errno));

    args.meth2_fp = gzopen(meth2_fname, "w");
    if ( args.meth2_fp == NULL)
        error("%s : %s", meth2_fname, strerror(errno));

    LOG_print("Load methylated data %s", bed_fname);
    args.hash = load_meth_bed(bed_fname);
    LOG_print("Load methylated data success.");
    // sort_meth_hash(args.hash);
    
    return 0;
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    generate_simulate_data();

    release_memory();
    
    return 0;
}
