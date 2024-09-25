#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <gmp.h>
#include <string>
#include <stdbool.h>
#include <pthread.h>
#include <vector>
#include <atomic>

#include "util.h"
#include "bloom.h"
#include "gmpecc.h"
#include "base58/libbase58.h"
#include "hash/ripemd160.h"
#include "hash/sha256.h"

#define MAX_MEMORY_GB 6
#define GB_TO_BYTES (1024ULL * 1024 * 1024)
#define MAX_BLOOM_BITS (MAX_MEMORY_GB * 8ULL * GB_TO_BYTES)
#define MAX_PUBLIC_KEYS 10
#define INITIAL_SEGMENTS 100

const char *version = "0.4.20240901";
const char *EC_constant_N = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
const char *EC_constant_P = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f";
const char *EC_constant_Gx = "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798";
const char *EC_constant_Gy = "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8";

const char *formats[3] = {"publickey", "rmd160", "address"};
const char *looks[2] = {"compress", "uncompress"};

// Function prototypes
void showhelp();
void set_format(char *param);
void set_look(char *param);
void set_bit(char *param);
void set_publickey(char *param, int index);
void set_range(char *param);
void generate_straddress(struct Point *publickey, bool compress, char *dst);
void generate_strrmd160(struct Point *publickey, bool compress, char *dst);
void generate_strpublickey(struct Point *publickey, bool compress, char *dst);
int init_bloom_filter(const char *filename);
void segmented_subtraction(mpz_t start, mpz_t end, int segments);
bool check_and_subtract(mpz_t start, mpz_t end, uint64_t *checks_performed);
void *random_search_thread(void *arg);
void parallel_random_search();

struct bloom *bf;

char *str_output = NULL;

char str_publickey[131];
char str_rmd160[41];
char str_address[41];

int num_threads = 1;
std::atomic<bool> found_match(false);
pthread_mutex_t gmp_mutex = PTHREAD_MUTEX_INITIALIZER;

struct ThreadData
{
    mpz_t start;
    mpz_t end;
    uint64_t checks_performed;
    int thread_id;
};

struct Point target_publickeys[MAX_PUBLIC_KEYS],
    base_publickey, sum_publickey, dst_publickey;
int num_public_keys = 0;

int FLAG_RANGE = 0;
int FLAG_BIT = 0;
int FLAG_RANDOM = 0;
int FLAG_PUBLIC = 0;
int FLAG_FORMART = 0;
int FLAG_HIDECOMMENT = 0;
int FLAG_LOOK = 0;
int FLAG_MODE = 0;
int FLAG_N;
int FLAG_OVERWRITE = 0;
int FLAG_SUBTRACTION = 0;
uint64_t N = 0;

mpz_t min_range, max_range, diff, TWO, base_key, sum_key, dst_key;
gmp_randstate_t state;

int main(int argc, char **argv)
{
    FILE *OUTPUT;
    char c;
    uint64_t i = 0;

    // Initialize EC constants
    mpz_init_set_str(EC.p, EC_constant_P, 16);
    mpz_init_set_str(EC.n, EC_constant_N, 16);
    mpz_init_set_str(G.x, EC_constant_Gx, 16);
    mpz_init_set_str(G.y, EC_constant_Gy, 16);
    init_doublingG(&G);

    // Initialize variables
    mpz_init(min_range);
    mpz_init(max_range);
    mpz_init(diff);
    mpz_init_set_ui(TWO, 2);

    // Initialize bloom filter and load public keys
    if (init_bloom_filter("public_keys.txt") != 0)
    {
        fprintf(stderr, "Failed to initialize Bloom filter\n");
        return 1;
    }

    // Parse command-line arguments
    while ((c = getopt(argc, argv, "hvxeRsb:n:o:p:r:f:l:t:")) != -1)
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
            set_bit((char *)optarg);
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
            if (num_public_keys < MAX_PUBLIC_KEYS)
            {
                set_publickey((char *)optarg, num_public_keys);
                num_public_keys++;
                FLAG_PUBLIC = 1;
            }
            else
            {
                fprintf(stderr, "[E] Maximum number of public keys reached\n");
                exit(0);
            }
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
        case 'e':
            FLAG_OVERWRITE = 1;
            break;
        case 's':
            FLAG_SUBTRACTION = 1;
            break;
        case 't':
            num_threads = atoi(optarg);
            if (num_threads <= 0)
            {
                fprintf(stderr, "[E] Invalid number of threads: %s\n", optarg);
                exit(0);
            }
            break;
        }
    }

    // Ensure necessary flags are set
    if ((FLAG_BIT || FLAG_RANGE) && FLAG_PUBLIC)
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

        // Calculate the difference between max and min range
        mpz_sub(diff, max_range, min_range);

        // Initialize variables for public key operations
        mpz_init(base_publickey.x);
        mpz_init(base_publickey.y);
        mpz_init(sum_publickey.x);
        mpz_init(sum_publickey.y);
        mpz_init(dst_publickey.x);
        mpz_init(dst_publickey.y);
        mpz_init(base_key);
        mpz_init(sum_key);

        if (FLAG_SUBTRACTION)
        {
            segmented_subtraction(min_range, max_range, INITIAL_SEGMENTS);
        }
        else if (FLAG_RANDOM)
        {
            parallel_random_search();
        }
        else
        {
            mpz_cdiv_q_ui(base_key, diff, N);
            mpz_set(sum_key, min_range);

            for (i = 0; i < N; i++)
            {
                Scalar_Multiplication(G, &base_publickey, sum_key);

                for (int j = 0; j < num_public_keys; j++)
                {
                    Point_Addition(&base_publickey, &target_publickeys[j], &dst_publickey);

                    switch (FLAG_FORMART)
                    {
                    case 0:
                        generate_strpublickey(&dst_publickey, FLAG_LOOK == 0, str_publickey);
                        if (bloom_check(bf, str_publickey, strlen(str_publickey)))
                        {
                            printf("Match found for public key %d!\n", j);
                            printf("Public key: %s\n", str_publickey);
                            gmp_printf("Subtraction: %Zd\n", sum_key);
                            goto cleanup;
                        }
                        if (FLAG_OVERWRITE)
                        {
                            fprintf(OUTPUT, "\r");
                            if (FLAG_HIDECOMMENT)
                            {
                                fprintf(OUTPUT, "%s", str_publickey);
                            }
                            else
                            {
                                gmp_fprintf(OUTPUT, "%s # - %Zd (Key %d)", str_publickey, sum_key, j);
                            }
                            fflush(OUTPUT);
                        }
                        else
                        {
                            if (FLAG_HIDECOMMENT)
                            {
                                fprintf(OUTPUT, "%s\n", str_publickey);
                            }
                            else
                            {
                                gmp_fprintf(OUTPUT, "%s # - %Zd (Key %d)\n", str_publickey, sum_key, j);
                            }
                        }
                        break;
                    case 1:
                        generate_strrmd160(&dst_publickey, FLAG_LOOK == 0, str_rmd160);
                        if (FLAG_OVERWRITE)
                        {
                            fprintf(OUTPUT, "\r");
                            if (FLAG_HIDECOMMENT)
                            {
                                fprintf(OUTPUT, "%s", str_rmd160);
                            }
                            else
                            {
                                gmp_fprintf(OUTPUT, "%s # - %Zd (Key %d)", str_rmd160, sum_key, j);
                            }
                            fflush(OUTPUT);
                        }
                        else
                        {
                            if (FLAG_HIDECOMMENT)
                            {
                                fprintf(OUTPUT, "%s\n", str_rmd160);
                            }
                            else
                            {
                                gmp_fprintf(OUTPUT, "%s # - %Zd (Key %d)\n", str_rmd160, sum_key, j);
                            }
                        }
                        break;
                    case 2:
                        generate_straddress(&dst_publickey, FLAG_LOOK == 0, str_address);
                        if (FLAG_OVERWRITE)
                        {
                            fprintf(OUTPUT, "\r");
                            if (FLAG_HIDECOMMENT)
                            {
                                fprintf(OUTPUT, "%s", str_address);
                            }
                            else
                            {
                                gmp_fprintf(OUTPUT, "%s # - %Zd (Key %d)", str_address, sum_key, j);
                            }
                            fflush(OUTPUT);
                        }
                        else
                        {
                            if (FLAG_HIDECOMMENT)
                            {
                                fprintf(OUTPUT, "%s\n", str_address);
                            }
                            else
                            {
                                gmp_fprintf(OUTPUT, "%s # - %Zd (Key %d)\n", str_address, sum_key, j);
                            }
                        }
                        break;
                    }
                }

                mpz_add(sum_key, sum_key, base_key);
            }
        }

        // Handle the final public key output based on format
        if (FLAG_OVERWRITE)
        {
            fprintf(OUTPUT, "\n");
        }
        for (int j = 0; j < num_public_keys; j++)
        {
            switch (FLAG_FORMART)
            {
            case 0:
                generate_strpublickey(&target_publickeys[j], FLAG_LOOK == 0, str_publickey);
                if (FLAG_HIDECOMMENT)
                {
                    fprintf(OUTPUT, "%s\n", str_publickey);
                }
                else
                {
                    fprintf(OUTPUT, "%s # target (Key %d)\n", str_publickey, j);
                }
                break;
            case 1:
                generate_strrmd160(&target_publickeys[j], FLAG_LOOK == 0, str_rmd160);
                if (FLAG_HIDECOMMENT)
                {
                    fprintf(OUTPUT, "%s\n", str_rmd160);
                }
                else
                {
                    fprintf(OUTPUT, "%s # target (Key %d)\n", str_rmd160, j);
                }
                break;
            case 2:
                generate_straddress(&target_publickeys[j], FLAG_LOOK == 0, str_address);
                if (FLAG_HIDECOMMENT)
                {
                    fprintf(OUTPUT, "%s\n", str_address);
                }
                else
                {
                    fprintf(OUTPUT, "%s # target (Key %d)\n", str_address, j);
                }
                break;
            }
        }
    }
    else
    {
        printf("[E] Invalid mode, please check N, bit, range, or public key...\n");
        showhelp();
    }

cleanup:
    if (OUTPUT != stdout)
    {
        fclose(OUTPUT);
    }

    // Clear memory
    bloom_free(bf);
    free(bf);
    mpz_clear(TWO);
    for (int j = 0; j < num_public_keys; j++)
    {
        mpz_clear(target_publickeys[j].x);
        mpz_clear(target_publickeys[j].y);
    }
    mpz_clear(base_publickey.x);
    mpz_clear(base_publickey.y);
    mpz_clear(sum_publickey.x);
    mpz_clear(sum_publickey.y);
    mpz_clear(dst_publickey.x);
    mpz_clear(dst_publickey.y);
    mpz_clear(base_key);
    mpz_clear(sum_key);

    return 0;
}

void showhelp()
{
    printf("\nUsage:\n-h\t\tshow this help\n");
    printf("-b bits\t\tFor some puzzles you only need a bit range\n");
    printf("-f format\tOutput format <publickey, rmd160, address>. Default: publickey\n");
    printf("-l look\t\tOutput <compress, uncompress>. Default: compress\n");
    printf("-n number\tNumber of public keys to be generated, this number will be even\n");
    printf("-o file\t\tOutput file, if you omit this option the output will go to the standard output\n");
    printf("-p key\t\tPublic key to be subtracted (compress or uncompress). Can be used multiple times for multiple keys\n");
    printf("-r A:B\t\trange A to B\n");
    printf("-R\t\tSet the public key subtraction Random instead of sequential\n");
    printf("-s\t\tUse segmented subtraction mode\n");
    printf("-x\t\tExclude comment\n");
    printf("-e\t\tOverwrite output instead of printing on new lines\n");
    printf("-t threads\tNumber of threads to use for parallel processing\n\n");
    printf("Developed by albertobsd, modified for multiple keys, segmented subtraction, and multi-threading\n\n");
}

void set_bit(char *param)
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

void set_publickey(char *param, int index)
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
    mpz_init(target_publickeys[index].x);
    mpz_init(target_publickeys[index].y);
    switch (len)
    {
    case 66:
        mpz_set_str(target_publickeys[index].x, dest + 2, 16);
        break;
    case 130:
        memcpy(hexvalue, dest + 2, 64);
        mpz_set_str(target_publickeys[index].x, hexvalue, 16);
        memcpy(hexvalue, dest + 66, 64);
        mpz_set_str(target_publickeys[index].y, hexvalue, 16);
        break;
    }
    if (mpz_cmp_ui(target_publickeys[index].y, 0) == 0)
    {
        mpz_t mpz_aux, mpz_aux2, Ysquared;
        mpz_init(mpz_aux);
        mpz_init(mpz_aux2);
        mpz_init(Ysquared);
        mpz_pow_ui(mpz_aux, target_publickeys[index].x, 3);
        mpz_add_ui(mpz_aux2, mpz_aux, 7);
        mpz_mod(Ysquared, mpz_aux2, EC.p);
        mpz_add_ui(mpz_aux, EC.p, 1);
        mpz_fdiv_q_ui(mpz_aux2, mpz_aux, 4);
        mpz_powm(target_publickeys[index].y, Ysquared, mpz_aux2, EC.p);
        mpz_sub(mpz_aux, EC.p, target_publickeys[index].y);
        switch (dest[1])
        {
        case '2':
            if (mpz_tstbit(target_publickeys[index].y, 0) == 1)
            {
                mpz_set(target_publickeys[index].y, mpz_aux);
            }
            break;
        case '3':
            if (mpz_tstbit(target_publickeys[index].y, 0) == 0)
            {
                mpz_set(target_publickeys[index].y, mpz_aux);
            }
            break;
        default:
            fprintf(stderr, "[E] Some invalid bit in the public key: %s\n", dest);
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
    unsigned char bin_publickey[65];
    unsigned char bin_sha256[32];
    unsigned char bin_rmd160[20];

    // Clear the destination buffer
    memset(dst, 0, 41); // 40 hex characters + null terminator

    // Step 1: Create the public key binary representation
    if (compress)
    {
        bin_publickey[0] = mpz_tstbit(publickey->y, 0) == 0 ? 0x02 : 0x03;
        mpz_export(bin_publickey + 1, NULL, 1, 32, 1, 0, publickey->x);
        sha256(bin_publickey, 33, bin_sha256);
    }
    else
    {
        bin_publickey[0] = 0x04;
        mpz_export(bin_publickey + 1, NULL, 1, 32, 1, 0, publickey->x);
        mpz_export(bin_publickey + 33, NULL, 1, 32, 1, 0, publickey->y);
        sha256(bin_publickey, 65, bin_sha256);
    }

    // Step 2: Perform RIPEMD-160 hash on the SHA-256 result
    ripemd160(bin_sha256, 32, bin_rmd160);

    // Step 3: Convert the binary RIPEMD-160 hash to hexadecimal string
    for (int i = 0; i < 20; i++)
    {
        sprintf(dst + (i * 2), "%02x", bin_rmd160[i]);
    }
}

void generate_straddress(struct Point *publickey, bool compress, char *dst)
{
    unsigned char bin_publickey[65];
    unsigned char bin_sha256[32];
    unsigned char bin_ripemd160[20];
    unsigned char bin_address[25];
    size_t pubaddress_size = 42;

    // Clear the destination buffer
    memset(dst, 0, 42);

    // Step 1: Create the public key binary representation
    if (compress)
    {
        bin_publickey[0] = mpz_tstbit(publickey->y, 0) == 0 ? 0x02 : 0x03;
        mpz_export(bin_publickey + 1, NULL, 1, 32, 1, 0, publickey->x);
        sha256(bin_publickey, 33, bin_sha256);
    }
    else
    {
        bin_publickey[0] = 0x04;
        mpz_export(bin_publickey + 1, NULL, 1, 32, 1, 0, publickey->x);
        mpz_export(bin_publickey + 33, NULL, 1, 32, 1, 0, publickey->y);
        sha256(bin_publickey, 65, bin_sha256);
    }

    // Step 2: Perform RIPEMD-160 hash on the SHA-256 result
    ripemd160(bin_sha256, 32, bin_ripemd160);

    // Step 3: Add version byte in front (0x00 for Main Network)
    bin_address[0] = 0x00;
    memcpy(bin_address + 1, bin_ripemd160, 20);

    // Step 4: Perform double SHA-256 hash on the extended RIPEMD-160 result
    sha256(bin_address, 21, bin_sha256);
    sha256(bin_sha256, 32, bin_sha256);

    // Step 5: Take the first 4 bytes of the second SHA-256 hash for checksum
    memcpy(bin_address + 21, bin_sha256, 4);

    // Step 6: Base58 encode the result
    if (!b58enc(dst, &pubaddress_size, bin_address, 25))
    {
        fprintf(stderr, "Error: Base58 encoding failed\n");
        // Handle the error as appropriate for your application
    }
}

int init_bloom_filter(const char *filename)
{
    uint64_t max_elements = 1000000001;           // Initial estimate
    double false_positive_rate = 0.0000000000001; // Very low false positive rate

    bf = (struct bloom *)malloc(sizeof(struct bloom));
    if (bf == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for Bloom filter\n");
        return 1;
    }

    // Calculate optimal Bloom filter size (in bits) based on max_elements and false_positive_rate
    uint64_t optimal_bits = (uint64_t)(-((double)max_elements * log(false_positive_rate)) / (log(2) * log(2)));

    // Ensure we don't exceed MAX_BLOOM_BITS
    uint64_t actual_bits = (optimal_bits < MAX_BLOOM_BITS) ? optimal_bits : MAX_BLOOM_BITS;

    // Recalculate max_elements based on actual_bits to maintain desired false positive rate
    max_elements = (uint64_t)((actual_bits * (log(2) * log(2))) / -log(false_positive_rate));

    if (bloom_init2(bf, max_elements, false_positive_rate) != 0)
    {
        fprintf(stderr, "Error initializing Bloom filter\n");
        free(bf);
        return 1;
    }

    printf("Bloom filter initialized with %llu bits (%.2f GB)\n",
           (unsigned long long)actual_bits, (double)actual_bits / (8.0 * GB_TO_BYTES));
    printf("Maximum elements: %llu\n", (unsigned long long)max_elements);
    printf("False positive rate: %.8f%%\n", false_positive_rate * 100);

    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Error opening file");
        bloom_free(bf);
        free(bf);
        return 1;
    }

    char line[256];
    uint64_t added_elements = 0;
    uint64_t total_lines = 0;

    // Count total lines in the file
    while (fgets(line, sizeof(line), file) != NULL)
    {
        total_lines++;
    }
    rewind(file);

    while (fgets(line, sizeof(line), file) && added_elements < max_elements)
    {
        line[strcspn(line, "\n")] = 0; // Remove newline
        bloom_add(bf, line, strlen(line));
        added_elements++;

        if (added_elements % 1000 == 0)
        {
            double progress = (double)added_elements / total_lines * 100;
            printf("\rLoading Bloom filter: %.2f%% (%llu / %llu)",
                   progress, (unsigned long long)added_elements, (unsigned long long)total_lines);
            fflush(stdout);
        }
    }

    fclose(file);
    printf("\nFinished loading %llu public keys into the Bloom filter\n", (unsigned long long)added_elements);

    return 0;
}

void segmented_subtraction(mpz_t start, mpz_t end, int segments)
{
    mpz_t range, segment_size, current, step;
    mpz_init(range);
    mpz_init(segment_size);
    mpz_init(current);
    mpz_init(step);

    mpz_sub(range, end, start);
    mpz_cdiv_q_ui(segment_size, range, segments);
    mpz_set(current, start);

    uint64_t total_checks = 0;
    uint64_t last_report = 0;

    while (mpz_cmp(current, end) <= 0)
    {
        for (int i = 0; i < segments; i++)
        {
            mpz_add(step, current, segment_size);
            if (mpz_cmp(step, end) > 0)
            {
                mpz_set(step, end);
            }

            uint64_t checks_performed = 0;
            if (!check_and_subtract(current, step, &checks_performed))
            {
                // If no match found, recursively search between current and step
                segmented_subtraction(current, step, segments);
            }
            total_checks += checks_performed;

            // Report progress every 1,000,000 checks
            if (total_checks - last_report > 1000000)
            {
                gmp_printf("Progress: Checked up to %Zd\n", current);
                last_report = total_checks;
            }

            if (mpz_cmp(step, end) == 0)
            {
                break;
            }
            mpz_set(current, step);
        }

        if (mpz_cmp(current, end) == 0)
        {
            break;
        }
    }

    gmp_printf("Completed search. Total checks performed: %llu\n", total_checks);

    mpz_clear(range);
    mpz_clear(segment_size);
    mpz_clear(current);
    mpz_clear(step);
}

bool check_and_subtract(mpz_t start, mpz_t end, uint64_t *checks_performed)
{
    mpz_t current;
    mpz_init(current);
    mpz_set(current, start);

    *checks_performed = 0;

    while (mpz_cmp(current, end) <= 0)
    {
        Scalar_Multiplication(G, &base_publickey, current);

        for (int j = 0; j < num_public_keys; j++)
        {
            Point_Addition(&base_publickey, &target_publickeys[j], &dst_publickey);

            switch (FLAG_FORMART)
            {
            case 0:
                generate_strpublickey(&dst_publickey, FLAG_LOOK == 0, str_publickey);
                if (bloom_check(bf, str_publickey, strlen(str_publickey)))
                {
                    printf("Match found for public key %d!\n", j);
                    printf("Public key: %s\n", str_publickey);
                    gmp_printf("Subtraction: %Zd\n", current);
                    mpz_clear(current);
                    return true;
                }
                break;
            case 1:
                generate_strrmd160(&dst_publickey, FLAG_LOOK == 0, str_rmd160);
                if (bloom_check(bf, str_rmd160, strlen(str_rmd160)))
                {
                    printf("Match found for public key %d!\n", j);
                    printf("RMD160: %s\n", str_rmd160);
                    gmp_printf("Subtraction: %Zd\n", current);
                    mpz_clear(current);
                    return true;
                }
                break;
            case 2:
                generate_straddress(&dst_publickey, FLAG_LOOK == 0, str_address);
                if (bloom_check(bf, str_address, strlen(str_address)))
                {
                    printf("Match found for public key %d!\n", j);
                    printf("Address: %s\n", str_address);
                    gmp_printf("Subtraction: %Zd\n", current);
                    mpz_clear(current);
                    return true;
                }
                break;
            }
        }

        mpz_add_ui(current, current, 1);
        (*checks_performed)++;

        // Optionally, you can add more frequent progress updates here
        if (*checks_performed % 10000 == 0)
        {
            gmp_printf("Checked: %Zd\r", current);
            fflush(stdout);
        }
    }

    mpz_clear(current);
    return false;
}

void *random_search_thread(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    mpz_t random_key, thread_diff;
    mpz_init(random_key);
    mpz_init(thread_diff);

    gmp_randstate_t thread_state;
    pthread_mutex_lock(&gmp_mutex);
    gmp_randinit_mt(thread_state);
    gmp_randseed_ui(thread_state, time(NULL) + data->thread_id);
    pthread_mutex_unlock(&gmp_mutex);

    mpz_sub(thread_diff, data->end, data->start);

    struct Point base_publickey, dst_publickey;
    mpz_init(base_publickey.x);
    mpz_init(base_publickey.y);
    mpz_init(dst_publickey.x);
    mpz_init(dst_publickey.y);

    char thread_output[256];

    while (!found_match.load())
    {
        pthread_mutex_lock(&gmp_mutex);
        mpz_urandomm(random_key, thread_state, thread_diff);
        mpz_add(random_key, random_key, data->start);
        pthread_mutex_unlock(&gmp_mutex);

        Scalar_Multiplication(G, &base_publickey, random_key);

        for (int j = 0; j < num_public_keys; j++)
        {
            Point_Addition(&base_publickey, &target_publickeys[j], &dst_publickey);

            bool match_found = false;
            switch (FLAG_FORMART)
            {
            case 0:
                generate_strpublickey(&dst_publickey, FLAG_LOOK == 0, str_publickey);
                match_found = bloom_check(bf, str_publickey, strlen(str_publickey));
                break;
            case 1:
                generate_strrmd160(&dst_publickey, FLAG_LOOK == 0, str_rmd160);
                match_found = bloom_check(bf, str_rmd160, strlen(str_rmd160));
                break;
            case 2:
                generate_straddress(&dst_publickey, FLAG_LOOK == 0, str_address);
                match_found = bloom_check(bf, str_address, strlen(str_address));
                break;
            }

            if (match_found)
            {
                found_match.store(true);
                printf("\nMatch found by thread %d for public key %d!\n", data->thread_id, j);
                switch (FLAG_FORMART)
                {
                case 0:
                    printf("Public key: %s\n", str_publickey);
                    break;
                case 1:
                    printf("RMD160: %s\n", str_rmd160);
                    break;
                case 2:
                    printf("Address: %s\n", str_address);
                    break;
                }
                gmp_printf("Subtraction: %Zd\n", random_key);
                mpz_clear(random_key);
                mpz_clear(thread_diff);
                gmp_randclear(thread_state);
                return NULL;
            }
        }

        // Display current key being processed
        pthread_mutex_lock(&gmp_mutex);
        switch (FLAG_FORMART)
        {
        case 0:
            gmp_snprintf(thread_output, sizeof(thread_output), "\rThread %d: %s # - %Zd",
                         data->thread_id, str_publickey, random_key);
            break;
        case 1:
            gmp_snprintf(thread_output, sizeof(thread_output), "\rThread %d: %s # - %Zd",
                         data->thread_id, str_rmd160, random_key);
            break;
        case 2:
            gmp_snprintf(thread_output, sizeof(thread_output), "\rThread %d: %s # - %Zd",
                         data->thread_id, str_address, random_key);
            break;
        }
        printf("%s", thread_output);
        fflush(stdout);
        pthread_mutex_unlock(&gmp_mutex);
    }

    mpz_clear(random_key);
    mpz_clear(thread_diff);
    gmp_randclear(thread_state);
    return NULL;
}

void parallel_random_search()
{
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> thread_data(num_threads);

    mpz_t thread_range;
    mpz_init(thread_range);
    mpz_sub(thread_range, max_range, min_range);
    mpz_fdiv_q_ui(thread_range, thread_range, num_threads);

    for (int i = 0; i < num_threads; i++)
    {
        mpz_init(thread_data[i].start);
        mpz_init(thread_data[i].end);

        mpz_mul_ui(thread_data[i].start, thread_range, i);
        mpz_add(thread_data[i].start, thread_data[i].start, min_range);

        if (i == num_threads - 1)
        {
            mpz_set(thread_data[i].end, max_range);
        }
        else
        {
            mpz_mul_ui(thread_data[i].end, thread_range, i + 1);
            mpz_add(thread_data[i].end, thread_data[i].end, min_range);
        }

        thread_data[i].thread_id = i;

        pthread_create(&threads[i], NULL, random_search_thread, &thread_data[i]);
    }

    for (int i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    mpz_clear(thread_range);

    if (!found_match.load())
    {
        printf("\nNo match found. Search terminated.\n");
    }
}