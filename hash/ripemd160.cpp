/*
 * This file is part of the VanitySearch distribution (https://github.com/JeanLucPons/VanitySearch).
 * Copyright (c) 2019 Jean Luc PONS.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "ripemd160.h"
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace _ripemd160 {

static const uint32_t IV[5] = {
    0x67452301ul, 0xEFCDAB89ul, 0x98BADCFEul, 0x10325476ul, 0xC3D2E1F0ul
};

#ifdef _MSC_VER
#define ROL(x, n) _rotl(x, n)
#else
#define ROL(x, n) (((x) << (n)) | ((x) >> (32 - (n))))
#endif

#define f1(x, y, z) ((x) ^ (y) ^ (z))
#define f2(x, y, z) (((x) & (y)) | (~(x) & (z)))
#define f3(x, y, z) (((x) | ~(y)) ^ (z))
#define f4(x, y, z) (((x) & (z)) | ((y) & ~(z)))
#define f5(x, y, z) ((x) ^ ((y) | ~(z)))

#define R(a, b, c, d, e, f, x, k, r) \
    a = ROL(a + f(b, c, d) + x + k, r) + e; \
    c = ROL(c, 10);

#define R1(a,b,c,d,e,x,r)  R(a, b, c, d, e, f1, x, 0x00000000ul, r)
#define R2(a,b,c,d,e,x,r)  R(a, b, c, d, e, f2, x, 0x5A827999ul, r)
#define R3(a,b,c,d,e,x,r)  R(a, b, c, d, e, f3, x, 0x6ED9EBA1ul, r)
#define R4(a,b,c,d,e,x,r)  R(a, b, c, d, e, f4, x, 0x8F1BBCDCul, r)
#define R5(a,b,c,d,e,x,r)  R(a, b, c, d, e, f5, x, 0xA953FD4Eul, r)

#define R1_(a,b,c,d,e,x,r) R(a, b, c, d, e, f5, x, 0x50A28BE6ul, r)
#define R2_(a,b,c,d,e,x,r) R(a, b, c, d, e, f4, x, 0x5C4DD124ul, r)
#define R3_(a,b,c,d,e,x,r) R(a, b, c, d, e, f3, x, 0x6D703EF3ul, r)
#define R4_(a,b,c,d,e,x,r) R(a, b, c, d, e, f2, x, 0x7A6D76E9ul, r)
#define R5_(a,b,c,d,e,x,r) R(a, b, c, d, e, f1, x, 0x00000000ul, r)

static inline void Transform(uint32_t* s, const uint32_t* chunk)
{
    uint32_t a1 = s[0], b1 = s[1], c1 = s[2], d1 = s[3], e1 = s[4];
    uint32_t a2 = a1, b2 = b1, c2 = c1, d2 = d1, e2 = e1;

    // Rounds 1-16
    R1(a1,b1,c1,d1,e1,chunk[ 0],11); R1_(a2,b2,c2,d2,e2,chunk[ 5], 8);
    R1(e1,a1,b1,c1,d1,chunk[ 1],14); R1_(e2,a2,b2,c2,d2,chunk[14], 9);
    R1(d1,e1,a1,b1,c1,chunk[ 2],15); R1_(d2,e2,a2,b2,c2,chunk[ 7], 9);
    R1(c1,d1,e1,a1,b1,chunk[ 3],12); R1_(c2,d2,e2,a2,b2,chunk[ 0],11);
    R1(b1,c1,d1,e1,a1,chunk[ 4], 5); R1_(b2,c2,d2,e2,a2,chunk[ 9],13);
    R1(a1,b1,c1,d1,e1,chunk[ 5], 8); R1_(a2,b2,c2,d2,e2,chunk[ 2],15);
    R1(e1,a1,b1,c1,d1,chunk[ 6], 7); R1_(e2,a2,b2,c2,d2,chunk[11],15);
    R1(d1,e1,a1,b1,c1,chunk[ 7], 9); R1_(d2,e2,a2,b2,c2,chunk[ 4], 5);
    R1(c1,d1,e1,a1,b1,chunk[ 8],11); R1_(c2,d2,e2,a2,b2,chunk[13], 7);
    R1(b1,c1,d1,e1,a1,chunk[ 9],13); R1_(b2,c2,d2,e2,a2,chunk[ 6], 7);
    R1(a1,b1,c1,d1,e1,chunk[10],14); R1_(a2,b2,c2,d2,e2,chunk[15], 8);
    R1(e1,a1,b1,c1,d1,chunk[11],15); R1_(e2,a2,b2,c2,d2,chunk[ 8],11);
    R1(d1,e1,a1,b1,c1,chunk[12], 6); R1_(d2,e2,a2,b2,c2,chunk[ 1],14);
    R1(c1,d1,e1,a1,b1,chunk[13], 7); R1_(c2,d2,e2,a2,b2,chunk[10],14);
    R1(b1,c1,d1,e1,a1,chunk[14], 9); R1_(b2,c2,d2,e2,a2,chunk[ 3],12);
    R1(a1,b1,c1,d1,e1,chunk[15], 8); R1_(a2,b2,c2,d2,e2,chunk[12], 6);

    // Rounds 17-32
    R2(e1,a1,b1,c1,d1,chunk[ 7], 7); R2_(e2,a2,b2,c2,d2,chunk[ 6], 9);
    R2(d1,e1,a1,b1,c1,chunk[ 4], 6); R2_(d2,e2,a2,b2,c2,chunk[11],13);
    R2(c1,d1,e1,a1,b1,chunk[13], 8); R2_(c2,d2,e2,a2,b2,chunk[ 3],15);
    R2(b1,c1,d1,e1,a1,chunk[ 1],13); R2_(b2,c2,d2,e2,a2,chunk[ 7], 7);
    R2(a1,b1,c1,d1,e1,chunk[10],11); R2_(a2,b2,c2,d2,e2,chunk[ 0],12);
    R2(e1,a1,b1,c1,d1,chunk[ 6], 9); R2_(e2,a2,b2,c2,d2,chunk[13], 8);
    R2(d1,e1,a1,b1,c1,chunk[15], 7); R2_(d2,e2,a2,b2,c2,chunk[ 5], 9);
    R2(c1,d1,e1,a1,b1,chunk[ 3],15); R2_(c2,d2,e2,a2,b2,chunk[10],11);
    R2(b1,c1,d1,e1,a1,chunk[12], 7); R2_(b2,c2,d2,e2,a2,chunk[14], 7);
    R2(a1,b1,c1,d1,e1,chunk[ 0],12); R2_(a2,b2,c2,d2,e2,chunk[15], 7);
    R2(e1,a1,b1,c1,d1,chunk[ 9],15); R2_(e2,a2,b2,c2,d2,chunk[ 8],12);
    R2(d1,e1,a1,b1,c1,chunk[ 5], 9); R2_(d2,e2,a2,b2,c2,chunk[12], 7);
    R2(c1,d1,e1,a1,b1,chunk[ 2],11); R2_(c2,d2,e2,a2,b2,chunk[ 4], 6);
    R2(b1,c1,d1,e1,a1,chunk[14], 7); R2_(b2,c2,d2,e2,a2,chunk[ 9],15);
    R2(a1,b1,c1,d1,e1,chunk[11],13); R2_(a2,b2,c2,d2,e2,chunk[ 1],13);
    R2(e1,a1,b1,c1,d1,chunk[ 8],12); R2_(e2,a2,b2,c2,d2,chunk[ 2],11);

    // Rounds 33-48
    R3(d1,e1,a1,b1,c1,chunk[ 3],11); R3_(d2,e2,a2,b2,c2,chunk[15], 9);
    R3(c1,d1,e1,a1,b1,chunk[10],13); R3_(c2,d2,e2,a2,b2,chunk[ 5], 7);
    R3(b1,c1,d1,e1,a1,chunk[14], 6); R3_(b2,c2,d2,e2,a2,chunk[ 1],15);
    R3(a1,b1,c1,d1,e1,chunk[ 4], 7); R3_(a2,b2,c2,d2,e2,chunk[ 3],11);
    R3(e1,a1,b1,c1,d1,chunk[ 9],14); R3_(e2,a2,b2,c2,d2,chunk[ 7], 8);
    R3(d1,e1,a1,b1,c1,chunk[15], 9); R3_(d2,e2,a2,b2,c2,chunk[14], 6);
    R3(c1,d1,e1,a1,b1,chunk[ 8],13); R3_(c2,d2,e2,a2,b2,chunk[ 6], 6);
    R3(b1,c1,d1,e1,a1,chunk[ 1],15); R3_(b2,c2,d2,e2,a2,chunk[ 9],14);
    R3(a1,b1,c1,d1,e1,chunk[ 2],14); R3_(a2,b2,c2,d2,e2,chunk[11],12);
    R3(e1,a1,b1,c1,d1,chunk[ 7], 8); R3_(e2,a2,b2,c2,d2,chunk[ 8],13);
    R3(d1,e1,a1,b1,c1,chunk[ 0],13); R3_(d2,e2,a2,b2,c2,chunk[12], 5);
    R3(c1,d1,e1,a1,b1,chunk[ 6], 6); R3_(c2,d2,e2,a2,b2,chunk[ 2],14);
    R3(b1,c1,d1,e1,a1,chunk[13], 5); R3_(b2,c2,d2,e2,a2,chunk[10],13);
    R3(a1,b1,c1,d1,e1,chunk[11],12); R3_(a2,b2,c2,d2,e2,chunk[ 0],13);
    R3(e1,a1,b1,c1,d1,chunk[ 5], 7); R3_(e2,a2,b2,c2,d2,chunk[ 4], 7);
    R3(d1,e1,a1,b1,c1,chunk[12], 5); R3_(d2,e2,a2,b2,c2,chunk[13], 5);

    // Rounds 49-64
    R4(c1,d1,e1,a1,b1,chunk[ 1],11); R4_(c2,d2,e2,a2,b2,chunk[ 8],15);
    R4(b1,c1,d1,e1,a1,chunk[ 9],12); R4_(b2,c2,d2,e2,a2,chunk[ 6], 5);
    R4(a1,b1,c1,d1,e1,chunk[11],14); R4_(a2,b2,c2,d2,e2,chunk[ 4], 8);
    R4(e1,a1,b1,c1,d1,chunk[10],15); R4_(e2,a2,b2,c2,d2,chunk[ 1],11);
    R4(d1,e1,a1,b1,c1,chunk[ 0],14); R4_(d2,e2,a2,b2,c2,chunk[ 3],14);
    R4(c1,d1,e1,a1,b1,chunk[ 8],15); R4_(c2,d2,e2,a2,b2,chunk[11],14);
    R4(c1,d1,e1,a1,b1,chunk[ 8],15); R4_(c2,d2,e2,a2,b2,chunk[11],14);
    R4(b1,c1,d1,e1,a1,chunk[12], 9); R4_(b2,c2,d2,e2,a2,chunk[15], 6);
    R4(a1,b1,c1,d1,e1,chunk[ 4], 8); R4_(a2,b2,c2,d2,e2,chunk[ 0],14);
    R4(e1,a1,b1,c1,d1,chunk[13], 9); R4_(e2,a2,b2,c2,d2,chunk[ 5], 6);
    R4(d1,e1,a1,b1,c1,chunk[ 3],14); R4_(d2,e2,a2,b2,c2,chunk[12], 9);
    R4(c1,d1,e1,a1,b1,chunk[ 7], 5); R4_(c2,d2,e2,a2,b2,chunk[ 2],12);
    R4(b1,c1,d1,e1,a1,chunk[15], 6); R4_(b2,c2,d2,e2,a2,chunk[13], 9);
    R4(a1,b1,c1,d1,e1,chunk[14], 8); R4_(a2,b2,c2,d2,e2,chunk[ 9],12);
    R4(e1,a1,b1,c1,d1,chunk[ 5], 6); R4_(e2,a2,b2,c2,d2,chunk[ 7], 5);
    R4(d1,e1,a1,b1,c1,chunk[ 6], 5); R4_(d2,e2,a2,b2,c2,chunk[10],15);
    R4(c1,d1,e1,a1,b1,chunk[ 2],12); R4_(c2,d2,e2,a2,b2,chunk[14], 8);

    // Rounds 65-80
    R5(b1,c1,d1,e1,a1,chunk[ 4], 9); R5_(b2,c2,d2,e2,a2,chunk[12], 8);
    R5(a1,b1,c1,d1,e1,chunk[ 0],15); R5_(a2,b2,c2,d2,e2,chunk[15], 5);
    R5(e1,a1,b1,c1,d1,chunk[ 5], 5); R5_(e2,a2,b2,c2,d2,chunk[10],12);
    R5(d1,e1,a1,b1,c1,chunk[ 9],11); R5_(d2,e2,a2,b2,c2,chunk[ 4], 9);
    R5(c1,d1,e1,a1,b1,chunk[ 7], 6); R5_(c2,d2,e2,a2,b2,chunk[ 1],12);
    R5(b1,c1,d1,e1,a1,chunk[12], 8); R5_(b2,c2,d2,e2,a2,chunk[ 5], 5);
    R5(a1,b1,c1,d1,e1,chunk[ 2],13); R5_(a2,b2,c2,d2,e2,chunk[ 8],14);
    R5(e1,a1,b1,c1,d1,chunk[10],12); R5_(e2,a2,b2,c2,d2,chunk[ 7], 6);
    R5(d1,e1,a1,b1,c1,chunk[14], 5); R5_(d2,e2,a2,b2,c2,chunk[ 6], 8);
    R5(c1,d1,e1,a1,b1,chunk[ 1],12); R5_(c2,d2,e2,a2,b2,chunk[ 2],13);
    R5(b1,c1,d1,e1,a1,chunk[ 3],13); R5_(b2,c2,d2,e2,a2,chunk[13], 6);
    R5(a1,b1,c1,d1,e1,chunk[ 8],14); R5_(a2,b2,c2,d2,e2,chunk[14], 5);
    R5(e1,a1,b1,c1,d1,chunk[11],11); R5_(e2,a2,b2,c2,d2,chunk[ 0],15);
    R5(d1,e1,a1,b1,c1,chunk[ 6], 8); R5_(d2,e2,a2,b2,c2,chunk[ 3],13);
    R5(c1,d1,e1,a1,b1,chunk[15], 5); R5_(c2,d2,e2,a2,b2,chunk[ 9],11);
    R5(b1,c1,d1,e1,a1,chunk[13], 6); R5_(b2,c2,d2,e2,a2,chunk[11],11);

    uint32_t t = s[1] + c1 + d2;
    s[1] = s[2] + d1 + e2;
    s[2] = s[3] + e1 + a2;
    s[3] = s[4] + a1 + b2;
    s[4] = s[0] + b1 + c2;
    s[0] = t;
}

} // namespace _ripemd160

CRIPEMD160::CRIPEMD160() : bytes(0) {
    memcpy(s, _ripemd160::IV, 20);
}

void CRIPEMD160::Write(const unsigned char* data, size_t len) {
    const unsigned char* end = data + len;
    size_t bufsize = bytes % 64;
    if (bufsize && bufsize + len >= 64) {
        memcpy(buf + bufsize, data, 64 - bufsize);
        _ripemd160::Transform(s, reinterpret_cast<const uint32_t*>(buf));
        data += 64 - bufsize;
        bytes += 64 - bufsize;
        bufsize = 0;
    }
    while (end >= data + 64) {
        _ripemd160::Transform(s, reinterpret_cast<const uint32_t*>(data));
        data += 64;
        bytes += 64;
    }
    if (end > data) {
        memcpy(buf + bufsize, data, end - data);
        bytes += end - data;
    }
}

void CRIPEMD160::Finalize(unsigned char hash[20]) {
    static const unsigned char pad[64] = {0x80};
    uint64_t total_bits = bytes << 3;
    uint32_t bufsize = bytes % 64;
    
    Write(pad, 1 + ((119 - bufsize) % 64));
    
    uint64_t size_le = total_bits;
    #ifdef __LITTLE_ENDIAN__
    Write(reinterpret_cast<const unsigned char*>(&size_le), 8);
    #else
    unsigned char size_bytes[8];
    size_bytes[0] = size_le >> 0;
    size_bytes[1] = size_le >> 8;
    size_bytes[2] = size_le >> 16;
    size_bytes[3] = size_le >> 24;
    size_bytes[4] = size_le >> 32;
    size_bytes[5] = size_le >> 40;
    size_bytes[6] = size_le >> 48;
    size_bytes[7] = size_le >> 56;
    Write(size_bytes, 8);
    #endif
    
    memcpy(hash, s, 20);
}

void ripemd160_32(unsigned char *input, unsigned char *digest) {
    uint32_t s[5];
    memcpy(s, _ripemd160::IV, 20);
    
    uint32_t chunk[16];
    memcpy(chunk, input, 32);
    chunk[8] = 0x80;
    memset(chunk + 9, 0, 6 * sizeof(uint32_t));
    chunk[14] = 32 << 3;
    chunk[15] = 0;
    
    _ripemd160::Transform(s, chunk);
    memcpy(digest, s, 20);
}

void ripemd160(unsigned char *input, int length, unsigned char *digest) {
    CRIPEMD160 cripe;
    cripe.Write(input, length);
    cripe.Finalize(digest);
}

std::string ripemd160_hex(unsigned char *digest) {
    static const char hex_digits[] = "0123456789abcdef";
    std::string result;
    result.reserve(40);
    for (int i = 0; i < 20; ++i) {
        result.push_back(hex_digits[digest[i] >> 4]);
        result.push_back(hex_digits[digest[i] & 0xf]);
    }
    return result;
}