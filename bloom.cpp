#include <assert.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <inttypes.h>
#include <sys/types.h>
#include <unistd.h>

#include "bloom.h"
#include "xxhash.h"

#define MAKESTRING(n) STRING(n)
#define STRING(n) #n
#define BLOOM_MAGIC "libbloom2"
#define BLOOM_VERSION_MAJOR 2
#define BLOOM_VERSION_MINOR 202

class OptimizedBloom
{
private:
    uint64_t *bits;
    uint64_t *secondaryBits;
    uint64_t *tertiaryBits;
    size_t numBits;
    size_t numSecondaryBits;
    size_t numTertiaryBits;
    size_t numHashes;
    static const size_t BITS_PER_ELEMENT = 64;

public:
    OptimizedBloom(size_t expectedElements, double falsePositiveRate)
    {
        numBits = static_cast<size_t>(-log(falsePositiveRate) * expectedElements / (log(2) * log(2)));
        numBits = ((numBits + BITS_PER_ELEMENT - 1) / BITS_PER_ELEMENT) * BITS_PER_ELEMENT;
        numSecondaryBits = numBits / 2;
        numTertiaryBits = numBits / 4;
        numHashes = static_cast<size_t>(log(2) * numBits / expectedElements);

        bits = new uint64_t[numBits / BITS_PER_ELEMENT]();
        secondaryBits = new uint64_t[numSecondaryBits / BITS_PER_ELEMENT]();
        tertiaryBits = new uint64_t[numTertiaryBits / BITS_PER_ELEMENT]();
    }

    ~OptimizedBloom()
    {
        delete[] bits;
        delete[] secondaryBits;
        delete[] tertiaryBits;
    }

    void add(const void *data, size_t len)
    {
        uint64_t hash1 = XXH64(data, len, 0x59f2815b16f81798);
        uint64_t hash2 = XXH64(data, len, hash1);
        uint64_t hash3 = XXH64(data, len, hash2);

        for (size_t i = 0; i < numHashes; ++i)
        {
            size_t bitIndex = (hash1 + i * hash2) % numBits;
            size_t elementIndex = bitIndex / BITS_PER_ELEMENT;
            uint64_t bitMask = 1ULL << (bitIndex % BITS_PER_ELEMENT);
            bits[elementIndex] |= bitMask;

            // Secondary filter
            bitIndex = (hash2 + i * hash3) % numSecondaryBits;
            elementIndex = bitIndex / BITS_PER_ELEMENT;
            bitMask = 1ULL << (bitIndex % BITS_PER_ELEMENT);
            secondaryBits[elementIndex] |= bitMask;

            // Tertiary filter
            bitIndex = (hash3 + i * hash1) % numTertiaryBits;
            elementIndex = bitIndex / BITS_PER_ELEMENT;
            bitMask = 1ULL << (bitIndex % BITS_PER_ELEMENT);
            tertiaryBits[elementIndex] |= bitMask;
        }
    }

    bool possiblyContains(const void *data, size_t len) const
    {
        uint64_t hash1 = XXH64(data, len, 0x59f2815b16f81798);
        uint64_t hash2 = XXH64(data, len, hash1);
        uint64_t hash3 = XXH64(data, len, hash2);

        for (size_t i = 0; i < numHashes; ++i)
        {
            // Primary check
            size_t bitIndex = (hash1 + i * hash2) % numBits;
            size_t elementIndex = bitIndex / BITS_PER_ELEMENT;
            uint64_t bitMask = 1ULL << (bitIndex % BITS_PER_ELEMENT);
            if (!(bits[elementIndex] & bitMask))
            {
                return false;
            }

            // Secondary check (more thorough)
            bitIndex = (hash2 + i * hash3) % numSecondaryBits;
            elementIndex = bitIndex / BITS_PER_ELEMENT;
            bitMask = 1ULL << (bitIndex % BITS_PER_ELEMENT);
            if (!(secondaryBits[elementIndex] & bitMask))
            {
                return false;
            }

            // Tertiary check (most thorough)
            bitIndex = (hash3 + i * hash1) % numTertiaryBits;
            elementIndex = bitIndex / BITS_PER_ELEMENT;
            bitMask = 1ULL << (bitIndex % BITS_PER_ELEMENT);
            if (!(tertiaryBits[elementIndex] & bitMask))
            {
                return false;
            }
        }
        return true;
    }

    void addBatch(const void **items, const int *lengths, int count)
    {
        for (int i = 0; i < count; ++i)
        {
            add(items[i], lengths[i]);
        }
    }

    size_t getTotalBitCount() const
    {
        return numBits + numSecondaryBits + numTertiaryBits;
    }

    size_t getHashCount() const
    {
        return numHashes;
    }
};

int bloom_init(struct bloom *bloom, uint64_t entries, long double error)
{
    return bloom_init2(bloom, entries, error);
}

int bloom_init2(struct bloom *bloom, uint64_t entries, long double error)
{
    memset(bloom, 0, sizeof(struct bloom));
    if (entries < 1000 || error <= 0 || error >= 1)
    {
        return 1;
    }

    bloom->entries = entries;
    bloom->error = error;

    long double num = -log(bloom->error);
    long double denom = 0.480453013918201; // ln(2)^2
    bloom->bpe = (num / denom);

    long double dentries = (long double)entries;
    long double allbits = dentries * bloom->bpe;
    bloom->bits = (uint64_t)allbits;

    bloom->bytes = (uint64_t)bloom->bits / 8;
    if (bloom->bits % 8)
    {
        bloom->bytes += 1;
    }

    bloom->hashes = (uint8_t)ceil(0.693147180559945 * bloom->bpe); // ln(2)

    bloom->bf = (uint8_t *)calloc(bloom->bytes, sizeof(uint8_t));
    if (bloom->bf == NULL)
    {
        return 1;
    }

    bloom->ready = 1;
    bloom->major = BLOOM_VERSION_MAJOR;
    bloom->minor = BLOOM_VERSION_MINOR;

    // Initialize OptimizedBloom
    bloom->optimizedBloom = new OptimizedBloom(entries, error);

    return 0;
}

int bloom_check(const struct bloom *bloom, const void *buffer, int len)
{
    if (bloom->ready == 0)
    {
        printf("bloom at %p not initialized!\n", (void *)bloom);
        return -1;
    }

    return static_cast<OptimizedBloom *>(bloom->optimizedBloom)->possiblyContains(buffer, len) ? 1 : 0;
}

int bloom_add(struct bloom *bloom, const void *buffer, int len)
{
    if (bloom->ready == 0)
    {
        printf("bloom at %p not initialized!\n", (void *)bloom);
        return -1;
    }

    if (!static_cast<OptimizedBloom *>(bloom->optimizedBloom)->possiblyContains(buffer, len))
    {
        static_cast<OptimizedBloom *>(bloom->optimizedBloom)->add(buffer, len);
        return 0; // Element was not present and was added
    }
    return 1; // Element (or a collision) was already present
}

void bloom_print(struct bloom *bloom)
{
    printf("bloom at %p\n", (void *)bloom);
    if (!bloom->ready)
    {
        printf(" *** NOT READY ***\n");
    }
    printf(" ->version = %d.%d\n", bloom->major, bloom->minor);
    printf(" ->entries = %" PRIu64 "\n", bloom->entries);
    printf(" ->error = %Lf\n", bloom->error);
    printf(" ->bits = %" PRIu64 "\n", bloom->bits);
    printf(" ->bits per elem = %f\n", bloom->bpe);
    printf(" ->bytes = %" PRIu64 "\n", bloom->bytes);
    unsigned int KB = bloom->bytes / 1024;
    unsigned int MB = KB / 1024;
    printf(" (%u KB, %u MB)\n", KB, MB);
    printf(" ->hash functions = %d\n", bloom->hashes);

    // Print OptimizedBloom stats
    OptimizedBloom *optimizedBloom = static_cast<OptimizedBloom *>(bloom->optimizedBloom);
    printf(" ->OptimizedBloom total bits = %zu\n", optimizedBloom->getTotalBitCount());
    printf(" ->OptimizedBloom hash functions = %zu\n", optimizedBloom->getHashCount());
}

void bloom_free(struct bloom *bloom)
{
    if (bloom->ready)
    {
        free(bloom->bf);
        delete static_cast<OptimizedBloom *>(bloom->optimizedBloom);
    }
    bloom->ready = 0;
}

int bloom_reset(struct bloom *bloom)
{
    if (!bloom->ready)
        return 1;
    memset(bloom->bf, 0, bloom->bytes);
    delete static_cast<OptimizedBloom *>(bloom->optimizedBloom);
    bloom->optimizedBloom = new OptimizedBloom(bloom->entries, bloom->error);
    return 0;
}

void bloom_add_batch(struct bloom *bloom, const void **items, const int *lengths, int count)
{
    if (!bloom->ready)
        return;
    static_cast<OptimizedBloom *>(bloom->optimizedBloom)->addBatch(items, lengths, count);
}

const char *bloom_version()
{
    return MAKESTRING(BLOOM_VERSION);
}