/*
Developed by Luis Alberto
email: alberto.bsd@gmail.com
Modified to include threading support and Bloom filter
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "util.h"
#include <errno.h>
#include <limits.h>

#include "gmpecc.h"
#include "base58/libbase58.h"
#include "rmd160.h"
#include "sha256.h"
#include "bloom.h"

#define MAX_BLOOM_SIZE (5ULL << 30) // 5 GB max size

const char *version = "0.3.20240928";
const char *EC_constant_N = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
const char *EC_constant_P = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f";
const char *EC_constant_Gx = "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798";
const char *EC_constant_Gy = "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8";

const char *formats[3] = {"publickey", "rmd160", "address"};
const char *looks[2] = {"compress", "uncompress"};

void showhelp();
void set_format(char *param);
void set_look(char *param);
void set_bit_range(char *param);
void set_publickey(char *param);
void set_range(char *param);
void generate_straddress(struct Point *publickey, bool compress, char *dst);
void generate_strrmd160(struct Point *publickey, bool compress, char *dst);
void generate_strpublickey(struct Point *publickey, bool compress, char *dst);
void *thread_process(void *arg);
int load_bloom_filter(const char *filename, BloomFilter **bloom);

char *str_output = NULL;
char *bloom_file = NULL;

char str_publickey[131];
char str_rmd160[41];
char str_address[41];

struct Point target_publickey, base_publickey, sum_publickey, negated_publickey, dst_publickey;

int FLAG_RANGE = 0;
int FLAG_BIT = 0;
int FLAG_RANDOM = 0;
int FLAG_PUBLIC = 0;
int FLAG_FORMART = 0;
int FLAG_HIDECOMMENT = 0;
int FLAG_LOOK = 0;
int FLAG_MODE = 0;
int FLAG_N;
uint64_t N = 0, M;

mpz_t min_range, max_range, diff, TWO, base_key, sum_key, dst_key;
gmp_randstate_t state;

// New global variables for threading and Bloom filter
int num_threads = 1;
pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t found_mutex = PTHREAD_MUTEX_INITIALIZER;
FILE *OUTPUT;
BloomFilter *bloom_filter = NULL;
int found_match = 0;
char found_publickey[131];
mpz_t found_subtraction;

// Thread argument structure
struct ThreadArg
{
    int thread_id;
    uint64_t start;
    uint64_t end;
};

int main(int argc, char **argv)
{
    char c;
    uint64_t i = 0;
    mpz_init_set_str(EC.p, EC_constant_P, 16);
    mpz_init_set_str(EC.n, EC_constant_N, 16);
    mpz_init_set_str(G.x, EC_constant_Gx, 16);
    mpz_init_set_str(G.y, EC_constant_Gy, 16);
    init_doublingG(&G);

    mpz_init(min_range);
    mpz_init(max_range);
    mpz_init(diff);
    mpz_init_set_ui(TWO, 2);
    mpz_init(target_publickey.x);
    mpz_init_set_ui(target_publickey.y, 0);
    mpz_init(found_subtraction);

    while ((c = getopt(argc, argv, "hvxRb:n:o:p:r:f:l:t:B:")) != -1)
    {
        switch (c)
        {
        case 'x':
            FLAG_HIDECOMMENT = 1;
            break;
        case 'h':
            showhelp();
            exit(0);
            break;
        case 'b':
            set_bit_range((char *)optarg);
            FLAG_BIT = 1;
            break;
        case 'n':
            N = strtol((char *)optarg, NULL, 10);
            if (N <= 0)
            {
                fprintf(stderr, "[E] invalid bit N number %s\n", optarg);
                exit(0);
            }
            FLAG_N = 1;
            break;
        case 'o':
            str_output = (char *)optarg;
            break;
        case 'p':
            set_publickey((char *)optarg);
            FLAG_PUBLIC = 1;
            break;
        case 'r':
            set_range((char *)optarg);
            FLAG_RANGE = 1;
            break;
        case 'R':
            FLAG_RANDOM = 1;
            break;
        case 'v':
            printf("version %s\n", version);
            exit(0);
            break;
        case 'l':
            set_look((char *)optarg);
            break;
        case 'f':
            set_format((char *)optarg);
            break;
        case 't':
            num_threads = atoi(optarg);
            if (num_threads <= 0)
            {
                fprintf(stderr, "[E] invalid number of threads: %s\n", optarg);
                exit(0);
            }
            break;
        case 'B':
            bloom_file = (char *)optarg;
            break;
        }
    }

    if ((FLAG_BIT || FLAG_RANGE) && FLAG_PUBLIC && FLAG_N && bloom_file)
    {
        if (str_output)
        {
            OUTPUT = fopen(str_output, "a");
            if (OUTPUT == NULL)
            {
                fprintf(stderr, "can't open file %s\n", str_output);
                OUTPUT = stdout;
            }
        }
        else
        {
            OUTPUT = stdout;
        }

        // Initialize and load Bloom filter
        if (load_bloom_filter(bloom_file, &bloom_filter) != 0)
        {
            fprintf(stderr, "Error loading Bloom filter\n");
            exit(1);
        }

        if (N % 2 == 1)
        {
            N++;
        }
        M = N / 2;
        mpz_sub(diff, max_range, min_range);
        mpz_init(base_publickey.x);
        mpz_init(base_publickey.y);
        mpz_init(sum_publickey.x);
        mpz_init(sum_publickey.y);
        mpz_init(negated_publickey.x);
        mpz_init(negated_publickey.y);
        mpz_init(dst_publickey.x);
        mpz_init(dst_publickey.y);
        mpz_init(base_key);
        mpz_init(sum_key);

        pthread_t threads[num_threads];
        struct ThreadArg thread_args[num_threads];

        uint64_t chunk_size = M / num_threads;
        uint64_t remainder = M % num_threads;

        for (i = 0; i < num_threads; i++)
        {
            thread_args[i].thread_id = i;
            thread_args[i].start = i * chunk_size;
            thread_args[i].end = (i + 1) * chunk_size;

            if (i == num_threads - 1)
            {
                thread_args[i].end += remainder;
            }

            pthread_create(&threads[i], NULL, thread_process, (void *)&thread_args[i]);
        }

        for (i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], NULL);
        }

        if (found_match)
        {
            printf("\nMatch found!\n");
            printf("%s - %s\n", found_publickey, mpz_get_str(NULL, 10, found_subtraction));
        }
        else
        {
            printf("\nNo match found.\n");
        }

        mpz_clear(base_publickey.x);
        mpz_clear(base_publickey.y);
        mpz_clear(sum_publickey.x);
        mpz_clear(sum_publickey.y);
        mpz_clear(negated_publickey.x);
        mpz_clear(negated_publickey.y);
        mpz_clear(dst_publickey.x);
        mpz_clear(dst_publickey.y);
        mpz_clear(base_key);
        mpz_clear(sum_key);
        mpz_clear(found_subtraction);

        if (OUTPUT != stdout)
        {
            fclose(OUTPUT);
        }

        free_bloom_filter(bloom_filter);
    }
    else
    {
        fprintf(stderr, "Version: %s\n", version);
        fprintf(stderr, "[E] there are some missing parameters\n");
        showhelp();
    }
    return 0;
}

void *thread_process(void *arg)
{
    struct ThreadArg *thread_arg = (struct ThreadArg *)arg;
    uint64_t i;
    mpz_t thread_base_key, thread_sum_key;
    struct Point thread_base_publickey, thread_sum_publickey, thread_negated_publickey, thread_dst_publickey;

    mpz_init(thread_base_key);
    mpz_init(thread_sum_key);
    mpz_init(thread_base_publickey.x);
    mpz_init(thread_base_publickey.y);
    mpz_init(thread_sum_publickey.x);
    mpz_init(thread_sum_publickey.y);
    mpz_init(thread_negated_publickey.x);
    mpz_init(thread_negated_publickey.y);
    mpz_init(thread_dst_publickey.x);
    mpz_init(thread_dst_publickey.y);

    char thread_publickey[131];

    if (FLAG_RANDOM)
    {
        gmp_randstate_t thread_state;
        gmp_randinit_mt(thread_state);
        gmp_randseed_ui(thread_state, ((int)clock()) + ((int)time(NULL)) + thread_arg->thread_id);

        for (i = thread_arg->start; i < thread_arg->end && !found_match; i++)
        {
            mpz_urandomm(thread_base_key, thread_state, diff);
            Scalar_Multiplication(G, &thread_base_publickey, thread_base_key);
            Point_Negation(&thread_base_publickey, &thread_negated_publickey);
            Point_Addition(&thread_base_publickey, &target_publickey, &thread_dst_publickey);

            generate_strpublickey(&thread_dst_publickey, FLAG_LOOK == 0, thread_publickey);

            if (contains(bloom_filter, thread_publickey))
            {
                pthread_mutex_lock(&found_mutex);
                if (!found_match)
                {
                    found_match = 1;
                    strcpy(found_publickey, thread_publickey);
                    mpz_set(found_subtraction, thread_base_key);
                }
                pthread_mutex_unlock(&found_mutex);
                break;
            }

            pthread_mutex_lock(&output_mutex);
            printf("\r%s - %s", thread_publickey, mpz_get_str(NULL, 10, thread_base_key));
            fflush(stdout);
            pthread_mutex_unlock(&output_mutex);

            if (found_match)
                break;

            Point_Addition(&thread_negated_publickey, &target_publickey, &thread_dst_publickey);
            generate_strpublickey(&thread_dst_publickey, FLAG_LOOK == 0, thread_publickey);

            if (contains(bloom_filter, thread_publickey))
            {
                pthread_mutex_lock(&found_mutex);
                if (!found_match)
                {
                    found_match = 1;
                    strcpy(found_publickey, thread_publickey);
                    mpz_neg(found_subtraction, thread_base_key);
                }
                pthread_mutex_unlock(&found_mutex);
                break;
            }

            pthread_mutex_lock(&output_mutex);
            printf("\r%s - %s", thread_publickey, mpz_get_str(NULL, 10, thread_base_key));
            fflush(stdout);
            pthread_mutex_unlock(&output_mutex);
        }

        gmp_randclear(thread_state);
    }
    else
    {
        mpz_cdiv_q_ui(thread_base_key, diff, M);
        mpz_mul_ui(thread_sum_key, thread_base_key, thread_arg->start);
        Scalar_Multiplication(G, &thread_sum_publickey, thread_sum_key);

        for (i = thread_arg->start; i < thread_arg->end && !found_match; i++)
        {
            Point_Negation(&thread_sum_publickey, &thread_negated_publickey);
            Point_Addition(&thread_sum_publickey, &target_publickey, &thread_dst_publickey);

            generate_strpublickey(&thread_dst_publickey, FLAG_LOOK == 0, thread_publickey);

            if (contains(bloom_filter, thread_publickey))
            {
                pthread_mutex_lock(&found_mutex);
                if (!found_match)
                {
                    found_match = 1;
                    strcpy(found_publickey, thread_publickey);
                    mpz_set(found_subtraction, thread_sum_key);
                }
                pthread_mutex_unlock(&found_mutex);
                break;
            }

            pthread_mutex_lock(&output_mutex);
            printf("\r%s - %s", thread_publickey, mpz_get_str(NULL, 10, thread_sum_key));
            fflush(stdout);
            pthread_mutex_unlock(&output_mutex);

            if (found_match)
                break;

            Point_Addition(&thread_negated_publickey, &target_publickey, &thread_dst_publickey);
            generate_strpublickey(&thread_dst_publickey, FLAG_LOOK == 0, thread_publickey);

            if (contains(bloom_filter, thread_publickey))
            {
                pthread_mutex_lock(&found_mutex);
                if (!found_match)
                {
                    found_match = 1;
                    strcpy(found_publickey, thread_publickey);
                    mpz_neg(found_subtraction, thread_sum_key);
                }
                pthread_mutex_unlock(&found_mutex);
                break;
            }

            pthread_mutex_lock(&output_mutex);
            printf("\r%s - %s", thread_publickey, mpz_get_str(NULL, 10, thread_sum_key));
            fflush(stdout);
            pthread_mutex_unlock(&output_mutex);

            Point_Addition(&thread_sum_publickey, &thread_base_publickey, &thread_dst_publickey);
            mpz_set(thread_sum_publickey.x, thread_dst_publickey.x);
            mpz_set(thread_sum_publickey.y, thread_dst_publickey.y);
            mpz_add(thread_sum_key, thread_sum_key, thread_base_key);
        }
    }

    mpz_clear(thread_base_key);
    mpz_clear(thread_sum_key);
    mpz_clear(thread_base_publickey.x);
    mpz_clear(thread_base_publickey.y);
    mpz_clear(thread_sum_publickey.x);
    mpz_clear(thread_sum_publickey.y);
    mpz_clear(thread_negated_publickey.x);
    mpz_clear(thread_negated_publickey.y);
    mpz_clear(thread_dst_publickey.x);
    mpz_clear(thread_dst_publickey.y);

    return NULL;
}

void showhelp()
{
    printf("\nUsage:\n-h\t\tshow this help\n");
    printf("-b bits\t\tSet bit range for the search\n");
    printf("-f format\tOutput format <publickey, rmd160, address>. Default: publickey\n");
    printf("-l look\t\tOutput <compress, uncompress>. Default: compress\n");
    printf("-n number\tNumber of publickeys to be generated, this number will be even\n");
    printf("-o file\t\tOutput file, if you omit this option the output will go to the standard output\n");
    printf("-p key\t\tPublickey to be subtracted compress or uncompress\n");
    printf("-r A:B\t\trange A to B\n");
    printf("-R\t\tSet the publickey subtraction Random instead of sequential\n");
    printf("-t threads\tNumber of threads to use for processing\n");
    printf("-B file\t\tBloom filter file containing public keys\n");
    printf("-x\t\tExclude comment\n\n");
    printf("Developed by albertobsd\n");
    printf("Modified to include threading support and Bloom filter\n\n");
}

void set_bit_range(char *param)
{
    mpz_t MPZAUX;
    int bitrange = strtol(param, NULL, 10);
    if (bitrange > 0 && bitrange <= 256)
    {
        mpz_init(MPZAUX);
        mpz_pow_ui(MPZAUX, TWO, bitrange - 1);
        mpz_set(min_range, MPZAUX);
        mpz_pow_ui(MPZAUX, TWO, bitrange);
        mpz_sub_ui(MPZAUX, MPZAUX, 1);
        mpz_set(max_range, MPZAUX);
        gmp_fprintf(stderr, "[+] Min range: %Zx\n", min_range);
        gmp_fprintf(stderr, "[+] Max range: %Zx\n", max_range);
        mpz_clear(MPZAUX);
    }
    else
    {
        fprintf(stderr, "[E] invalid bit param: %s\n", param);
        exit(0);
    }
}

void set_publickey(char *param)
{
    char hexvalue[65];
    char *dest;
    int len;
    len = strlen(param);
    dest = (char *)calloc(len + 1, 1);
    if (dest == NULL)
    {
        fprintf(stderr, "[E] Error calloc\n");
        exit(0);
    }
    memset(hexvalue, 0, 65);
    memcpy(dest, param, len);
    trim(dest, " \t\n\r");
    len = strlen(dest);
    switch (len)
    {
    case 66:
        mpz_set_str(target_publickey.x, dest + 2, 16);
        break;
    case 130:
        memcpy(hexvalue, dest + 2, 64);
        mpz_set_str(target_publickey.x, hexvalue, 16);
        memcpy(hexvalue, dest + 66, 64);
        mpz_set_str(target_publickey.y, hexvalue, 16);
        break;
    }
    if (mpz_cmp_ui(target_publickey.y, 0) == 0)
    {
        mpz_t mpz_aux, mpz_aux2, Ysquared;
        mpz_init(mpz_aux);
        mpz_init(mpz_aux2);
        mpz_init(Ysquared);
        mpz_pow_ui(mpz_aux, target_publickey.x, 3);
        mpz_add_ui(mpz_aux2, mpz_aux, 7);
        mpz_mod(Ysquared, mpz_aux2, EC.p);
        mpz_add_ui(mpz_aux, EC.p, 1);
        mpz_fdiv_q_ui(mpz_aux2, mpz_aux, 4);
        mpz_powm(target_publickey.y, Ysquared, mpz_aux2, EC.p);
        mpz_sub(mpz_aux, EC.p, target_publickey.y);
        switch (dest[1])
        {
        case '2':
            if (mpz_tstbit(target_publickey.y, 0) == 1)
            {
                mpz_set(target_publickey.y, mpz_aux);
            }
            break;
        case '3':
            if (mpz_tstbit(target_publickey.y, 0) == 0)
            {
                mpz_set(target_publickey.y, mpz_aux);
            }
            break;
        default:
            fprintf(stderr, "[E] Some invalid bit in the publickey: %s\n", dest);
            exit(0);
            break;
        }
        mpz_clear(mpz_aux);
        mpz_clear(mpz_aux2);
        mpz_clear(Ysquared);
    }
    free(dest);
}

void set_range(char *param)
{
    Tokenizer tk;
    char *dest;
    int len;
    len = strlen(param);
    dest = (char *)calloc(len + 1, 1);
    if (dest == NULL)
    {
        fprintf(stderr, "[E] Error calloc\n");
        exit(0);
    }
    memcpy(dest, param, len);
    dest[len] = '\0';
    stringtokenizer(dest, &tk);
    if (tk.n == 2)
    {
        mpz_init_set_str(min_range, nextToken(&tk), 16);
        mpz_init_set_str(max_range, nextToken(&tk), 16);
    }
    else
    {
        fprintf(stderr, "%i\n", tk.n);
        fprintf(stderr, "[E] Invalid range expected format A:B\n");
        exit(0);
    }
    freetokenizer(&tk);
    free(dest);
}

void set_format(char *param)
{
    int index = indexOf(param, formats, 3);
    if (index == -1)
    {
        fprintf(stderr, "[E] Unknown format: %s\n", param);
    }
    else
    {
        FLAG_FORMART = index;
    }
}

void set_look(char *param)
{
    int index = indexOf(param, looks, 2);
    if (index == -1)
    {
        fprintf(stderr, "[E] Unknown look: %s\n", param);
    }
    else
    {
        FLAG_LOOK = index;
    }
}

void generate_strpublickey(struct Point *publickey, bool compress, char *dst)
{
    memset(dst, 0, 132);
    if (compress)
    {
        if (mpz_tstbit(publickey->y, 0) == 0)
        { // Even
            gmp_snprintf(dst, 67, "02%0.64Zx", publickey->x);
        }
        else
        {
            gmp_snprintf(dst, 67, "03%0.64Zx", publickey->x);
        }
    }
    else
    {
        gmp_snprintf(dst, 131, "04%0.64Zx%0.64Zx", publickey->x, publickey->y);
    }
}

void generate_strrmd160(struct Point *publickey, bool compress, char *dst)
{
    char str_publickey[131];
    char bin_publickey[65];
    char bin_sha256[32];
    char bin_rmd160[20];
    memset(dst, 0, 42);
    if (compress)
    {
        if (mpz_tstbit(publickey->y, 0) == 0)
        { // Even
            gmp_snprintf(str_publickey, 67, "02%0.64Zx", publickey->x);
        }
        else
        {
            gmp_snprintf(str_publickey, 67, "03%0.64Zx", publickey->x);
        }
        hexs2bin(str_publickey, bin_publickey);
        sha256(bin_publickey, 33, bin_sha256);
    }
    else
    {
        gmp_snprintf(str_publickey, 131, "04%0.64Zx%0.64Zx", publickey->x, publickey->y);
        hexs2bin(str_publickey, bin_publickey);
        sha256(bin_publickey, 65, bin_sha256);
    }
    RMD160Data((const unsigned char *)bin_sha256, 32, bin_rmd160);
    tohex_dst(bin_rmd160, 20, dst);
}

void generate_straddress(struct Point *publickey, bool compress, char *dst)
{
    char str_publickey[131];
    char bin_publickey[65];
    char bin_sha256[32];
    char bin_digest[60];
    size_t pubaddress_size = 42;
    memset(dst, 0, 42);
    if (compress)
    {
        if (mpz_tstbit(publickey->y, 0) == 0)
        { // Even
            gmp_snprintf(str_publickey, 67, "02%0.64Zx", publickey->x);
        }
        else
        {
            gmp_snprintf(str_publickey, 67, "03%0.64Zx", publickey->x);
        }
        hexs2bin(str_publickey, bin_publickey);
        sha256(bin_publickey, 33, bin_sha256);
    }
    else
    {
        gmp_snprintf(str_publickey, 131, "04%0.64Zx%0.64Zx", publickey->x, publickey->y);
        hexs2bin(str_publickey, bin_publickey);
        sha256(bin_publickey, 65, bin_sha256);
    }
    RMD160Data((const unsigned char *)bin_sha256, 32, bin_digest + 1);

    /* First byte 0, this is for the Address beginning with 1.... */
    bin_digest[0] = 0;

    /* Double sha256 checksum */
    sha256(bin_digest, 21, bin_digest + 21);
    sha256(bin_digest + 21, 32, bin_digest + 21);

    /* Get the address */
    if (!b58enc(dst, &pubaddress_size, bin_digest, 25))
    {
        fprintf(stderr, "error b58enc\n");
    }
}

int load_bloom_filter(const char *filename, BloomFilter **bloom)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file: %s (Error: %s)\n", filename, strerror(errno));
        return 1;
    }

    // Count the total number of keys
    size_t total_keys = 0;
    char line[256];
    while (fgets(line, sizeof(line), file))
    {
        total_keys++;
        if (total_keys % 10000000 == 0)
        {
            printf("Counted %zu million keys...\n", total_keys / 1000000);
        }
    }
    printf("Total number of keys in file: %zu\n", total_keys);

    // Reset file pointer to the beginning
    if (fseek(file, 0, SEEK_SET) != 0)
    {
        fprintf(stderr, "Error resetting file pointer: %s\n", strerror(errno));
        fclose(file);
        return 1;
    }

    // Ask user how many lines to load
    size_t lines_to_load;
    printf("How many lines do you want to load into the Bloom filter? (0 for all): ");
    scanf("%zu", &lines_to_load);

    if (lines_to_load == 0 || lines_to_load > total_keys)
    {
        lines_to_load = total_keys;
    }

    // Ask for desired memory
    size_t size_in_gb;
    printf("Enter the desired Bloom filter size in GB: ");
    scanf("%zu", &size_in_gb);

    // Print Bloom filter options
    print_bloom_filter_options(size_in_gb, lines_to_load);

    // Ask user to choose an option
    size_t option;
    printf("Choose an option (1-%d): ", MAX_HASH_FUNCTIONS);
    scanf("%zu", &option);

    if (option < 1 || option > MAX_HASH_FUNCTIONS)
    {
        fprintf(stderr, "Invalid option.\n");
        fclose(file);
        return 1;
    }

    // Create Bloom filter
    *bloom = create_bloom_filter(size_in_gb, lines_to_load);
    if (*bloom == NULL)
    {
        fprintf(stderr, "Failed to create Bloom filter.\n");
        fclose(file);
        return 1;
    }

    (*bloom)->num_hash_functions = option;
    calculate_optimal_parameters(*bloom);

    printf("Initializing Bloom filter with parameters:\n");
    printf("  - Number of keys to load: %zu\n", lines_to_load);
    printf("  - Bloom filter size: %zu bits (%.2f GB)\n", (*bloom)->size, (float)(*bloom)->size / (8.0 * 1024 * 1024 * 1024));
    printf("  - Hash functions: %zu\n", (*bloom)->num_hash_functions);
    printf("  - Estimated false positive rate: %.20f\n", (*bloom)->false_positive_rate);

    // Add keys to the Bloom filter
    size_t count = 0;
    while (fgets(line, sizeof(line), file) && count < lines_to_load)
    {
        size_t len = strlen(line);
        if (len > 0 && line[len - 1] == '\n')
        {
            line[len - 1] = '\0';
            len--;
        }

        insert(*bloom, line);
        count++;

        if (count % 1000000 == 0)
        {
            printf("Added %zu million keys to the Bloom filter...\n", count / 1000000);
        }
    }

    fclose(file);

    printf("Bloom filter initialization complete:\n");
    printf("  - Keys loaded: %zu\n", count);
    printf("  - Bloom filter size: %.2f GB\n", (float)(*bloom)->size / (8.0 * 1024 * 1024 * 1024));
    printf("  - Hash functions: %zu\n", (*bloom)->num_hash_functions);
    printf("  - Actual false positive rate: %.20f\n", (*bloom)->false_positive_rate);

    return 0;
}
