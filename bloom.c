#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xxhash.h"
#include <float.h>

#define MAX_HASH_FUNCTIONS 20

typedef struct
{
    unsigned char *bits;
    size_t size;
    size_t num_elements;
    size_t num_hash_functions;
    double false_positive_rate;
} BloomFilter;

BloomFilter *create_bloom_filter(size_t size_in_gb, size_t num_elements)
{
    BloomFilter *filter = (BloomFilter *)malloc(sizeof(BloomFilter));
    if (!filter)
        return NULL;

    size_t size_in_bits = size_in_gb * 8ULL * 1024 * 1024 * 1024;
    filter->size = size_in_bits;
    filter->num_elements = num_elements;

    filter->bits = (unsigned char *)calloc((size_in_bits + 7) / 8, sizeof(unsigned char));
    if (!filter->bits)
    {
        free(filter);
        return NULL;
    }

    return filter;
}

void free_bloom_filter(BloomFilter *filter)
{
    if (filter)
    {
        free(filter->bits);
        free(filter);
    }
}

void set_bit(BloomFilter *filter, size_t index)
{
    filter->bits[index / 8] |= 1 << (index % 8);
}

int test_bit(BloomFilter *filter, size_t index)
{
    return (filter->bits[index / 8] & (1 << (index % 8))) != 0;
}

void calculate_optimal_parameters(BloomFilter *filter)
{
    double best_fpr = 1.0;
    size_t best_k = 1;

    double m = (double)filter->size;
    double n = (double)filter->num_elements;

    for (size_t k = 1; k <= MAX_HASH_FUNCTIONS; k++)
    {
        double exp_part = -((double)k * n) / m;
        double fpr;

        if (exp_part < -DBL_MAX_EXP)
        {
            fpr = 0.0; // Effectively zero
        }
        else
        {
            double inner = 1.0 - exp(exp_part);
            if (inner < 0.0)
                inner = 0.0; // Ensure non-negative
            fpr = pow(inner, (double)k);
        }

        if (fpr < best_fpr)
        {
            best_fpr = fpr;
            best_k = k;
        }
    }

    filter->num_hash_functions = best_k;
    filter->false_positive_rate = best_fpr;
}

void insert(BloomFilter *filter, const char *key)
{
    for (size_t i = 0; i < filter->num_hash_functions; i++)
    {
        XXH64_hash_t hash = XXH3_64bits_withSeed(key, strlen(key), i);
        size_t index = hash % filter->size;
        set_bit(filter, index);
    }
}

int contains(BloomFilter *filter, const char *key)
{
    for (size_t i = 0; i < filter->num_hash_functions; i++)
    {
        XXH64_hash_t hash = XXH3_64bits_withSeed(key, strlen(key), i);
        size_t index = hash % filter->size;
        if (!test_bit(filter, index))
        {
            return 0;
        }
    }
    return 1;
}

void print_bloom_filter_options(size_t size_in_gb, size_t num_elements)
{
    BloomFilter *filter = create_bloom_filter(size_in_gb, num_elements);
    if (!filter)
    {
        printf("Failed to create Bloom filter.\n");
        return;
    }

    printf("Bloom Filter Options:\n");
    printf("Size: %zu GB\n", size_in_gb);
    printf("Number of elements: %zu\n", num_elements);

    double size_in_bits = (double)filter->size;
    double m = size_in_bits;
    double n = (double)num_elements;

    for (size_t k = 1; k <= MAX_HASH_FUNCTIONS; k++)
    {
        filter->num_hash_functions = k;
        double exp_part = -((double)k * n) / m;
        double fpr;

        if (exp_part < -DBL_MAX_EXP)
        {
            fpr = 0.0; // Effectively zero
        }
        else
        {
            double inner = 1.0 - exp(exp_part);
            if (inner < 0.0)
                inner = 0.0; // Ensure non-negative
            fpr = pow(inner, (double)k);
        }

        printf("Option %zu: %zu hash functions, False Positive Rate: %.10f\n", k, k, fpr);
    }

    free_bloom_filter(filter);
}
