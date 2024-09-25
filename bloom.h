#ifndef _BLOOM_H
#define _BLOOM_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdint.h>

    struct bloom
    {
        uint64_t entries;
        uint64_t bits;
        uint64_t bytes;
        uint8_t hashes;
        long double error;
        uint8_t ready;
        uint8_t major;
        uint8_t minor;
        double bpe;
        uint8_t *bf;
        void *optimizedBloom;
    };

    int bloom_init2(struct bloom *bloom, uint64_t entries, long double error);
    int bloom_init(struct bloom *bloom, uint64_t entries, long double error);
    int bloom_check(const struct bloom *bloom, const void *buffer, int len);
    int bloom_add(struct bloom *bloom, const void *buffer, int len);
    void bloom_print(struct bloom *bloom);
    void bloom_free(struct bloom *bloom);
    int bloom_reset(struct bloom *bloom);
    const char *bloom_version();

    void bloom_add_batch(struct bloom *bloom, const void **items, const int *lengths, int count);

#ifdef __cplusplus
}
#endif

#endif