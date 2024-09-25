/*
 * This file is part of the VanitySearch distribution (https://github.com/JeanLucPons/VanitySearch).
 * Copyright (c) 2019 Jean Luc PONS.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string>
#include <immintrin.h>

#include "sha256.h"

#define BSWAP

/// Internal SHA-256 implementation.
namespace _sha256
{
    static const unsigned char pad[64] = { 0x80 };

    static const uint32_t K[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
        0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
        0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
        0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };

#ifndef WIN64
#define _byteswap_ulong __builtin_bswap32
#define _byteswap_uint64 __builtin_bswap64
#endif

#ifdef BSWAP
#define WRITEBE32(ptr,x) *((uint32_t *)(ptr)) = _byteswap_ulong(x)
#define WRITEBE64(ptr,x) *((uint64_t *)(ptr)) = _byteswap_uint64(x)
#define READBE32(ptr) (uint32_t)_byteswap_ulong(*(uint32_t *)(ptr))
#else
#define WRITEBE32(ptr,x) *(ptr) = x
#define WRITEBE64(ptr,x) *(ptr) = x
#define READBE32(ptr) *(uint32_t *)(ptr)
#endif

    inline void Initialize(uint32_t *s) {
        s[0] = 0x6a09e667ul;
        s[1] = 0xbb67ae85ul;
        s[2] = 0x3c6ef372ul;
        s[3] = 0xa54ff53aul;
        s[4] = 0x510e527ful;
        s[5] = 0x9b05688cul;
        s[6] = 0x1f83d9abul;
        s[7] = 0x5be0cd19ul;
    }

    void Transform(uint32_t* s, const unsigned char* chunk)
    {
        __m128i STATE0, STATE1;
        __m128i MSG, TMP;
        __m128i MSG0, MSG1, MSG2, MSG3;
        __m128i ABEF_SAVE, CDGH_SAVE;
        const __m128i MASK = _mm_set_epi64x(0x0c0d0e0f08090a0bULL, 0x0405060700010203ULL);

        // Load initial values
        TMP = _mm_loadu_si128((const __m128i*) &s[0]);
        STATE1 = _mm_loadu_si128((const __m128i*) &s[4]);
        
        TMP = _mm_shuffle_epi32(TMP, 0xB1);          // CDAB
        STATE1 = _mm_shuffle_epi32(STATE1, 0x1B);    // EFGH
        STATE0 = _mm_alignr_epi8(TMP, STATE1, 8);    // ABEF
        STATE1 = _mm_blend_epi16(STATE1, TMP, 0xF0); // CDGH

        // Save current state
        ABEF_SAVE = STATE0;
        CDGH_SAVE = STATE1;

        // Rounds 0-3
        MSG = _mm_loadu_si128((const __m128i*) &chunk[0]);
        MSG0 = _mm_shuffle_epi8(MSG, MASK);
        MSG = _mm_add_epi32(MSG0, _mm_set_epi64x(0xE9B5DBA5B5C0FBCFULL, 0x71374491428A2F98ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

        // Rounds 4-7
        MSG1 = _mm_loadu_si128((const __m128i*) &chunk[16]);
        MSG1 = _mm_shuffle_epi8(MSG1, MASK);
        MSG = _mm_add_epi32(MSG1, _mm_set_epi64x(0xAB1C5ED5923F82A4ULL, 0x59F111F13956C25BULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG0 = _mm_sha256msg1_epu32(MSG0, MSG1);

        // Rounds 8-11
        MSG2 = _mm_loadu_si128((const __m128i*) &chunk[32]);
        MSG2 = _mm_shuffle_epi8(MSG2, MASK);
        MSG = _mm_add_epi32(MSG2, _mm_set_epi64x(0x550C7DC3243185BEULL, 0x12835B01D807AA98ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG1 = _mm_sha256msg1_epu32(MSG1, MSG2);

        // Rounds 12-15
        MSG3 = _mm_loadu_si128((const __m128i*) &chunk[48]);
        MSG3 = _mm_shuffle_epi8(MSG3, MASK);
        MSG = _mm_add_epi32(MSG3, _mm_set_epi64x(0xC19BF1749BDC06A7ULL, 0x80DEB1FE72BE5D74ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG3, MSG2, 4);
        MSG0 = _mm_add_epi32(MSG0, TMP);
        MSG0 = _mm_sha256msg2_epu32(MSG0, MSG3);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG2 = _mm_sha256msg1_epu32(MSG2, MSG3);

        // Rounds 16-19
        MSG = _mm_add_epi32(MSG0, _mm_set_epi64x(0x240CA1CC0FC19DC6ULL, 0xEFBE4786E49B69C1ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG0, MSG3, 4);
        MSG1 = _mm_add_epi32(MSG1, TMP);
        MSG1 = _mm_sha256msg2_epu32(MSG1, MSG0);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG3 = _mm_sha256msg1_epu32(MSG3, MSG0);

        // Rounds 20-23
        MSG = _mm_add_epi32(MSG1, _mm_set_epi64x(0x76F988DA5CB0A9DCULL, 0x4A7484AA2DE92C6FULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG1, MSG0, 4);
        MSG2 = _mm_add_epi32(MSG2, TMP);
        MSG2 = _mm_sha256msg2_epu32(MSG2, MSG1);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG0 = _mm_sha256msg1_epu32(MSG0, MSG1);

        // Rounds 24-27
        MSG = _mm_add_epi32(MSG2, _mm_set_epi64x(0xBF597FC7B00327C8ULL, 0xA831C66D983E5152ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG2, MSG1, 4);
        MSG3 = _mm_add_epi32(MSG3, TMP);
        MSG3 = _mm_sha256msg2_epu32(MSG3, MSG2);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG1 = _mm_sha256msg1_epu32(MSG1, MSG2);

        // Rounds 28-31
        MSG = _mm_add_epi32(MSG3, _mm_set_epi64x(0x1429296706CA6351ULL, 0xD5A79147C6E00BF3ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG3, MSG2, 4);
        MSG0 = _mm_add_epi32(MSG0, TMP);
        MSG0 = _mm_sha256msg2_epu32(MSG0, MSG3);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG2 = _mm_sha256msg1_epu32(MSG2, MSG3);

        // Rounds 32-35
        MSG = _mm_add_epi32(MSG0, _mm_set_epi64x(0x53380D134D2C6DFCULL, 0x2E1B213827B70A85ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG0, MSG3, 4);
        MSG1 = _mm_add_epi32(MSG1, TMP);
        MSG1 = _mm_sha256msg2_epu32(MSG1, MSG0);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG3 = _mm_sha256msg1_epu32(MSG3, MSG0);

        // Rounds 36-39
        MSG = _mm_add_epi32(MSG1, _mm_set_epi64x(0x92722C8581C2C92EULL, 0x766A0ABB650A7354ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG1, MSG0, 4);
        MSG2 = _mm_add_epi32(MSG2, TMP);
        MSG2 = _mm_sha256msg2_epu32(MSG2, MSG1);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG0 = _mm_sha256msg1_epu32(MSG0, MSG1);

        // Rounds 40-43
        MSG = _mm_add_epi32(MSG2, _mm_set_epi64x(0xC76C51A3C24B8B70ULL, 0xA81A664BA2BFE8A1ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG2, MSG1, 4);
        MSG3 = _mm_add_epi32(MSG3, TMP);
        MSG3 = _mm_sha256msg2_epu32(MSG3, MSG2);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG1 = _mm_sha256msg1_epu32(MSG1, MSG2);

        // Rounds 44-47
        MSG = _mm_add_epi32(MSG3, _mm_set_epi64x(0x106AA070F40E3585ULL, 0xD6990624D192E819ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG3, MSG2, 4);
        MSG0 = _mm_add_epi32(MSG0, TMP);
        MSG0 = _mm_sha256msg2_epu32(MSG0, MSG3);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG2 = _mm_sha256msg1_epu32(MSG2, MSG3);

        // Rounds 48-51
        MSG = _mm_add_epi32(MSG0, _mm_set_epi64x(0x34B0BCB52748774CULL, 0x1E376C0819A4C116ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG0, MSG3, 4);
        MSG1 = _mm_add_epi32(MSG1, TMP);
        MSG1 = _mm_sha256msg2_epu32(MSG1, MSG0);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);
        MSG3 = _mm_sha256msg1_epu32(MSG3, MSG0);

        // Rounds 52-55
        MSG = _mm_add_epi32(MSG1, _mm_set_epi64x(0x682E6FF35B9CCA4FULL, 0x4ED8AA4A391C0CB3ULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG1, MSG0, 4);
        MSG2 = _mm_add_epi32(MSG2, TMP);
        MSG2 = _mm_sha256msg2_epu32(MSG2, MSG1);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

        // Rounds 56-59
        MSG = _mm_add_epi32(MSG2, _mm_set_epi64x(0x8CC7020884C87814ULL, 0x78A5636F748F82EEULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        TMP = _mm_alignr_epi8(MSG2, MSG1, 4);
        MSG3 = _mm_add_epi32(MSG3, TMP);
        MSG3 = _mm_sha256msg2_epu32(MSG3, MSG2);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

        // Rounds 60-63
        MSG = _mm_add_epi32(MSG3, _mm_set_epi64x(0xC67178F2BEF9A3F7ULL, 0xA4506CEB90BEFFFAULL));
        STATE1 = _mm_sha256rnds2_epu32(STATE1, STATE0, MSG);
        MSG = _mm_shuffle_epi32(MSG, 0x0E);
        STATE0 = _mm_sha256rnds2_epu32(STATE0, STATE1, MSG);

        // Combine state
        STATE0 = _mm_add_epi32(STATE0, ABEF_SAVE);
        STATE1 = _mm_add_epi32(STATE1, CDGH_SAVE);

        TMP = _mm_shuffle_epi32(STATE0, 0x1B);       // FEBA
        STATE1 = _mm_shuffle_epi32(STATE1, 0xB1);    // DCHG
        STATE0 = _mm_blend_epi16(TMP, STATE1, 0xF0); // DCBA
        STATE1 = _mm_alignr_epi8(STATE1, TMP, 8);    // ABEF

        // Save state
        _mm_storeu_si128((__m128i*) &s[0], STATE0);
        _mm_storeu_si128((__m128i*) &s[4], STATE1);
    }

    void Transform2(uint32_t* s, const unsigned char* chunk) {
        uint32_t t[8];
        memcpy(t, s, 32);
        Transform(t, chunk);
        s[0] = t[0];
    }
}

class CSHA256
{
private:
    alignas(16) uint32_t s[8];
    unsigned char buf[64];
    uint64_t bytes;

public:
    static const size_t OUTPUT_SIZE = 32;

    CSHA256();
    void Write(const unsigned char* data, size_t len);
    void Finalize(unsigned char hash[OUTPUT_SIZE]);
};

CSHA256::CSHA256() {
    bytes = 0;
    _sha256::Initialize(s);
}

void CSHA256::Write(const unsigned char* data, size_t len)
{
    const unsigned char* end = data + len;
    size_t bufsize = bytes % 64;
    if (bufsize && bufsize + len >= 64) {
        memcpy(buf + bufsize, data, 64 - bufsize);
        bytes += 64 - bufsize;
        data += 64 - bufsize;
        _sha256::Transform(s, buf);
        bufsize = 0;
    }
    while (end >= data + 64) {
        _sha256::Transform(s, data);
        bytes += 64;
        data += 64;
    }
    if (end > data) {
        memcpy(buf + bufsize, data, end - data);
        bytes += end - data;
    }
}

void CSHA256::Finalize(unsigned char hash[OUTPUT_SIZE])
{
    unsigned char sizedesc[8];
    WRITEBE64(sizedesc, bytes << 3);
    Write(_sha256::pad, 1 + ((119 - (bytes % 64)) % 64));
    Write(sizedesc, 8);
    WRITEBE32(hash, s[0]);
    WRITEBE32(hash + 4, s[1]);
    WRITEBE32(hash + 8, s[2]);
    WRITEBE32(hash + 12, s[3]);
    WRITEBE32(hash + 16, s[4]);
    WRITEBE32(hash + 20, s[5]);
    WRITEBE32(hash + 24, s[6]);
    WRITEBE32(hash + 28, s[7]);
}

void sha256(unsigned char *input, size_t length, unsigned char *digest) {
    CSHA256 sha;
    sha.Write(input, length);
    sha.Finalize(digest);
}

const uint8_t sizedesc_32[8] = { 0,0,0,0,0,0,1,0 };
const uint8_t sizedesc_33[8] = { 0,0,0,0,0,0,1,8 };
const uint8_t sizedesc_65[8] = { 0,0,0,0,0,0,2,8 };

void sha256_33(unsigned char *input, unsigned char *digest) {
    alignas(16) uint32_t s[8];
    _sha256::Initialize(s);
    memcpy(input + 33, _sha256::pad, 23);
    memcpy(input + 56, sizedesc_33, 8);
    _sha256::Transform(s, input);
    WRITEBE32(digest, s[0]);
    WRITEBE32(digest + 4, s[1]);
    WRITEBE32(digest + 8, s[2]);
    WRITEBE32(digest + 12, s[3]);
    WRITEBE32(digest + 16, s[4]);
    WRITEBE32(digest + 20, s[5]);
    WRITEBE32(digest + 24, s[6]);
    WRITEBE32(digest + 28, s[7]);
}

void sha256_65(unsigned char *input, unsigned char *digest) {
    alignas(16) uint32_t s[8];
    memcpy(input + 65, _sha256::pad, 55);
    memcpy(input + 120, sizedesc_65, 8);
    _sha256::Initialize(s);
    _sha256::Transform(s, input);
    _sha256::Transform(s, input+64);
    WRITEBE32(digest, s[0]);
    WRITEBE32(digest + 4, s[1]);
    WRITEBE32(digest + 8, s[2]);
    WRITEBE32(digest + 12, s[3]);
    WRITEBE32(digest + 16, s[4]);
    WRITEBE32(digest + 20, s[5]);
    WRITEBE32(digest + 24, s[6]);
    WRITEBE32(digest + 28, s[7]);
}

void sha256_checksum(uint8_t *input, int length, uint8_t *checksum) {
    alignas(16) uint32_t s[8];
    alignas(16) uint8_t b[64];
    memcpy(b, input, length);
    memcpy(b + length, _sha256::pad, 56 - length);
    WRITEBE64(b + 56, length << 3);
    _sha256::Transform2(s, b);
    WRITEBE32(checksum, s[0]);
}

std::string sha256_hex(unsigned char *digest) {
    char buf[2*32+1];
    buf[2*32] = 0;
    for (int i = 0; i < 32; i++)
        sprintf(buf+i*2,"%02x",digest[i]);
    return std::string(buf);
}

bool sha256_file(const char* file_name, uint8_t* checksum) {
    FILE* file = fopen(file_name, "rb");
    if (file == NULL) {
        printf("Failed to open file: %s\n", file_name);
        return false;
    }
    CSHA256 sha;
    uint8_t buffer[8192];
    size_t bytes_read;

    while ((bytes_read = fread(buffer, 1, sizeof(buffer), file)) > 0) {
        sha.Write(buffer, bytes_read);
    }

    sha.Finalize(checksum);
    fclose(file);
    return true;
}
