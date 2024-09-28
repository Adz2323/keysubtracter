#ifndef BLOOM_H
#define BLOOM_H

#include <stddef.h>

#define MAX_HASH_FUNCTIONS 20

typedef struct
{
    unsigned char *bits;
    size_t size;
    size_t num_elements;
    size_t num_hash_functions;
    double false_positive_rate;
} BloomFilter;

// Create a new Bloom filter
BloomFilter *create_bloom_filter(size_t size_in_gb, size_t num_elements);

// Free the memory allocated for the Bloom filter
void free_bloom_filter(BloomFilter *filter);

// Calculate optimal parameters for the Bloom filter
void calculate_optimal_parameters(BloomFilter *filter);

// Insert an element into the Bloom filter
void insert(BloomFilter *filter, const char *key);

// Check if an element might be in the Bloom filter
int contains(BloomFilter *filter, const char *key);

// Print Bloom filter options based on size and number of elements
void print_bloom_filter_options(size_t size_in_gb, size_t num_elements);

// Set a bit in the Bloom filter
void set_bit(BloomFilter *filter, size_t index);

// Test a bit in the Bloom filter
int test_bit(BloomFilter *filter, size_t index);

#endif // BLOOM_H
