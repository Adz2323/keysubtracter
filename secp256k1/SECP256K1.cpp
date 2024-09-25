#include <cstdio>
#include <cstring>
#include <immintrin.h>
#include "SECP256k1.h"
#include "Point.h"
#include "../util.h"
#include "../hash/sha256.h"
#include "../hash/ripemd160.h"
#include "Int.h"

#ifdef _MSC_VER
#define ALIGN16 __declspec(align(16))
#else
#define ALIGN16 __attribute__((aligned(16)))
#endif

Secp256K1::Secp256K1() {}

void Secp256K1::Init()
{
    P.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");
    Int::SetupField(&P);

    G.x.SetBase16("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798");
    G.y.SetBase16("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8");
    G.z.SetInt32(1);
    order.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141");

    Int::InitK1(&order);

    // Precompute Generator table
    Point N(G);
    for (int i = 0; i < 32; i++)
    {
        GTable[i * 256] = N;
        for (int j = 1; j < 255; j++)
        {
            GTable[i * 256 + j] = AddDirect(N, GTable[i * 256]);
        }
        GTable[i * 256 + 255] = N = DoubleDirect(N);
    }
}

Secp256K1::~Secp256K1() {}

Point Secp256K1::ComputePublicKey(Int *privKey)
{
    Point Q;
    Q.Clear();

    int i = 0;
    uint8_t b;

    // Find first significant byte
    while (i < 32 && !(b = privKey->GetByte(i)))
        i++;

    if (i < 32)
    {
        Q = GTable[256 * i + (b - 1)];
        i++;

        for (; i < 32; i++)
        {
            if ((b = privKey->GetByte(i)))
            {
                Q = Add2(Q, GTable[256 * i + (b - 1)]);
            }
        }
    }

    Q.Reduce();
    return Q;
}

inline Point Secp256K1::NextKey(Point &key)
{
    return AddDirect(key, G);
}

uint8_t Secp256K1::GetByte(const char *str, int idx)
{
    int val;
    if (sscanf(str + 2 * idx, "%2X", &val) != 1)
    {
        fprintf(stderr, "ParsePublicKeyHex: Error invalid public key specified (unexpected hexadecimal digit)\n");
        exit(-1);
    }
    return (uint8_t)val;
}

Point Secp256K1::Negation(Point &p)
{
    Point Q;
    Q.x.Set(&p.x);
    Q.y.Set(&this->P);
    Q.y.Sub(&p.y);
    Q.z.SetInt32(1);
    return Q;
}
bool Secp256K1::ParsePublicKeyHex(const char *str, Point &ret, bool &isCompressed)
{
    int len = strlen(str);
    ret.Clear();

    printf("Debug ParsePublicKeyHex: Input string = %s\n", str);
    printf("Debug ParsePublicKeyHex: String length = %d\n", len);

    if (len < 2)
    {
        fprintf(stderr, "ParsePublicKeyHex: Error invalid public key specified (66 or 130 character length)\n");
        return false;
    }

    uint8_t type = GetByte(str, 0);
    printf("Debug ParsePublicKeyHex: Key type = 0x%02X\n", type);

    switch (type)
    {
    case 0x02:
    case 0x03:
        if (len != 66)
        {
            fprintf(stderr, "ParsePublicKeyHex: Error invalid public key specified (66 character length)\n");
            return false;
        }
        for (int i = 0; i < 32; i++)
            ret.x.SetByte(31 - i, GetByte(str, i + 1));
        printf("Debug ParsePublicKeyHex: Parsed x = %s\n", ret.x.GetBase16());
        ret.y = GetY(ret.x, type == 0x02);
        printf("Debug ParsePublicKeyHex: Calculated y = %s\n", ret.y.GetBase16());
        isCompressed = true;
        break;

    case 0x04:
        if (len != 130)
        {
            fprintf(stderr, "ParsePublicKeyHex: Error invalid public key specified (130 character length)\n");
            return false;
        }
        for (int i = 0; i < 32; i++)
        {
            ret.x.SetByte(31 - i, GetByte(str, i + 1));
            ret.y.SetByte(31 - i, GetByte(str, i + 33));
        }
        printf("Debug ParsePublicKeyHex: Parsed x = %s\n", ret.x.GetBase16());
        printf("Debug ParsePublicKeyHex: Parsed y = %s\n", ret.y.GetBase16());
        isCompressed = false;
        break;

    default:
        fprintf(stderr, "ParsePublicKeyHex: Error invalid public key specified (Unexpected prefix (only 02,03 or 04 allowed)\n");
        return false;
    }

    ret.z.SetInt32(1);

    printf("Debug ParsePublicKeyHex: Before EC check: x = %s, y = %s\n", ret.x.GetBase16(), ret.y.GetBase16());

    if (!EC(ret))
    {
        fprintf(stderr, "ParsePublicKeyHex: Error invalid public key specified (Not lie on elliptic curve)\n");
        printf("Debug ParsePublicKeyHex: EC check failed\n");
        return false;
    }

    printf("Debug ParsePublicKeyHex: EC check passed\n");
    printf("Debug ParsePublicKeyHex: Final parsed point: x = %s, y = %s\n", ret.x.GetBase16(), ret.y.GetBase16());
    printf("Debug ParsePublicKeyHex: Is compressed: %s\n", isCompressed ? "true" : "false");

    return true;
}

char *Secp256K1::GetPublicKeyHex(bool compressed, Point &pubKey)
{
    unsigned char publicKeyBytes[65];
    char *ret = nullptr;
    if (!compressed)
    {
        publicKeyBytes[0] = 0x4;
        pubKey.x.Get32Bytes(publicKeyBytes + 1);
        pubKey.y.Get32Bytes(publicKeyBytes + 33);
        ret = (char *)tohex((char *)publicKeyBytes, 65);
    }
    else
    {
        publicKeyBytes[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
        pubKey.x.Get32Bytes(publicKeyBytes + 1);
        ret = (char *)tohex((char *)publicKeyBytes, 33);
    }
    return ret;
}

void Secp256K1::GetPublicKeyHex(bool compressed, Point &pubKey, char *dst)
{
    unsigned char publicKeyBytes[65];
    if (!compressed)
    {
        publicKeyBytes[0] = 0x4;
        pubKey.x.Get32Bytes(publicKeyBytes + 1);
        pubKey.y.Get32Bytes(publicKeyBytes + 33);
        tohex_dst((char *)publicKeyBytes, 65, dst);
    }
    else
    {
        publicKeyBytes[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
        pubKey.x.Get32Bytes(publicKeyBytes + 1);
        tohex_dst((char *)publicKeyBytes, 33, dst);
    }
}

char *Secp256K1::GetPublicKeyRaw(bool compressed, Point &pubKey)
{
    char *ret = (char *)malloc(compressed ? 33 : 65);
    if (ret == nullptr)
    {
        fprintf(stderr, "Can't alloc memory\n");
        exit(0);
    }
    if (!compressed)
    {
        ret[0] = 0x4;
        pubKey.x.Get32Bytes((unsigned char *)(ret + 1));
        pubKey.y.Get32Bytes((unsigned char *)(ret + 33));
    }
    else
    {
        ret[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
        pubKey.x.Get32Bytes((unsigned char *)(ret + 1));
    }
    return ret;
}

void Secp256K1::GetPublicKeyRaw(bool compressed, Point &pubKey, char *dst)
{
    if (!compressed)
    {
        dst[0] = 0x4;
        pubKey.x.Get32Bytes((unsigned char *)(dst + 1));
        pubKey.y.Get32Bytes((unsigned char *)(dst + 33));
    }
    else
    {
        dst[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
        pubKey.x.Get32Bytes((unsigned char *)(dst + 1));
    }
}

Point Secp256K1::AddDirect(Point &p1, Point &p2)
{
    Int _s, _p, dy, dx;
    Point r;
    r.z.SetInt32(1);

    dy.ModSub(&p2.y, &p1.y);
    dx.ModSub(&p2.x, &p1.x);
    dx.ModInv();
    _s.ModMulK1(&dy, &dx); // s = (p2.y - p1.y) * inverse(p2.x - p1.x);

    _p.ModSquareK1(&_s); // _p = pow2(s)

    r.x.ModSub(&_p, &p1.x);
    r.x.ModSub(&p2.x); // rx = pow2(s) - p1.x - p2.x;

    r.y.ModSub(&p2.x, &r.x);
    r.y.ModMulK1(&_s);
    r.y.ModSub(&p2.y); // ry = -p2.y - s * (ret.x - p2.x);

    return r;
}

Point Secp256K1::Add2(Point &p1, Point &p2)
{
    Int u, v, u1, v1, vs2, vs3, us2, a, us2w, vs2v2, vs3u2, _2vs2v2;
    Point r;
    u1.ModMulK1(&p2.y, &p1.z);
    v1.ModMulK1(&p2.x, &p1.z);
    u.ModSub(&u1, &p1.y);
    v.ModSub(&v1, &p1.x);
    us2.ModSquareK1(&u);
    vs2.ModSquareK1(&v);
    vs3.ModMulK1(&vs2, &v);
    us2w.ModMulK1(&us2, &p1.z);
    vs2v2.ModMulK1(&vs2, &p1.x);
    _2vs2v2.ModAdd(&vs2v2, &vs2v2);
    a.ModSub(&us2w, &vs3);
    a.ModSub(&_2vs2v2);

    r.x.ModMulK1(&v, &a);

    vs3u2.ModMulK1(&vs3, &p1.y);
    r.y.ModSub(&vs2v2, &a);
    r.y.ModMulK1(&r.y, &u);
    r.y.ModSub(&vs3u2);

    r.z.ModMulK1(&vs3, &p1.z);

    return r;
}

Point Secp256K1::Add(Point &p1, Point &p2)
{
    Int u, v, u1, u2, v1, v2, vs2, vs3, us2, w, a, us2w, vs2v2, vs3u2, _2vs2v2, x3, vs3y1;
    Point r;

    u1.ModMulK1(&p2.y, &p1.z);
    u2.ModMulK1(&p1.y, &p2.z);
    v1.ModMulK1(&p2.x, &p1.z);
    v2.ModMulK1(&p1.x, &p2.z);
    u.ModSub(&u1, &u2);
    v.ModSub(&v1, &v2);
    w.ModMulK1(&p1.z, &p2.z);
    us2.ModSquareK1(&u);
    vs2.ModSquareK1(&v);
    vs3.ModMulK1(&vs2, &v);
    us2w.ModMulK1(&us2, &w);
    vs2v2.ModMulK1(&vs2, &v2);
    _2vs2v2.ModAdd(&vs2v2, &vs2v2);
    a.ModSub(&us2w, &vs3);
    a.ModSub(&_2vs2v2);

    r.x.ModMulK1(&v, &a);

    vs3u2.ModMulK1(&vs3, &u2);
    r.y.ModSub(&vs2v2, &a);
    r.y.ModMulK1(&r.y, &u);
    r.y.ModSub(&vs3u2);

    r.z.ModMulK1(&vs3, &w);

    return r;
}

Point Secp256K1::DoubleDirect(Point &p)
{
    Int _s, _p, a;
    Point r;
    r.z.SetInt32(1);
    _s.ModMulK1(&p.x, &p.x);
    _p.ModAdd(&_s, &_s);
    _p.ModAdd(&_s);

    a.ModAdd(&p.y, &p.y);
    a.ModInv();
    _s.ModMulK1(&_p, &a); // s = (3 * pow2(p.x)) * inverse(2 * p.y);

    _p.ModMulK1(&_s, &_s);
    a.ModAdd(&p.x, &p.x);
    a.ModNeg();
    r.x.ModAdd(&a, &_p); // rx = pow2(s) + neg(2 * p.x);

    a.ModSub(&r.x, &p.x);

    _p.ModMulK1(&a, &_s);
    r.y.ModAdd(&_p, &p.y);
    r.y.ModNeg(); // ry = neg(p.y + s * (ret.x + neg(p.x)));
    return r;
}

Point Secp256K1::Double(Point &p)
{
    Int z2, x2, _3x2, w, s, s2, b, _8b, _8y2s2, y2, h;
    Point r;
    z2.ModSquareK1(&p.z);
    x2.ModSquareK1(&p.x);
    _3x2.ModAdd(&x2, &x2);
    _3x2.ModAdd(&x2);
    w.ModAdd(&z2, &_3x2);
    s.ModMulK1(&p.y, &p.z);
    b.ModMulK1(&p.y, &s);
    b.ModMulK1(&p.x);
    h.ModSquareK1(&w);
    _8b.ModAdd(&b, &b);
    _8b.ModDouble();
    _8b.ModDouble();
    h.ModSub(&_8b);
    r.x.ModMulK1(&h, &s);
    r.x.ModAdd(&r.x);
    s2.ModSquareK1(&s);
    y2.ModSquareK1(&p.y);
    _8y2s2.ModMulK1(&y2, &s2);
    _8y2s2.ModDouble();
    _8y2s2.ModDouble();
    _8y2s2.ModDouble();
    r.y.ModAdd(&b, &b);
    r.y.ModAdd(&b, &b);
    r.y.ModAdd(&r.y, &r.y);
    r.y.ModSub(&h);
    r.y.ModMulK1(&w);
    r.y.ModSub(&_8y2s2);
    r.z.ModMulK1(&s2, &s);
    r.z.ModDouble();
    r.z.ModDouble();
    r.z.ModDouble();
    return r;
}

Int Secp256K1::GetY(Int x, bool isEven)
{
    Int y, right;

    // Calculate right side: x^3 + 7 (mod p)
    right.Set(&x);
    right.ModSquareK1(&right); // x^2
    right.ModMulK1(&x);        // x^3
    right.Add(&_7);
    right.Mod(&P); // x^3 + 7 (mod p)

    printf("Debug GetY: x = %s\n", x.GetBase16());
    printf("Debug GetY: x^3 + 7 = %s\n", right.GetBase16());

    // Calculate y = (x^3 + 7)^((p+1)/4) mod p
    Int exp(P);
    exp.AddOne();
    exp.ShiftR(2); // (p+1)/4

    // Use ModExp correctly
    y.Set(&right);
    y.ModExp(&exp);
    y.Mod(&P);

    printf("Debug GetY: calculated y = %s\n", y.GetBase16());
    printf("Debug GetY: isEven parameter = %s\n", isEven ? "true" : "false");
    printf("Debug GetY: calculated y is even = %s\n", y.IsEven() ? "true" : "false");

    // Ensure the correct parity
    if (y.IsEven() != isEven)
    {
        y.Neg();
        y.Mod(&P);
        printf("Debug GetY: y negated\n");
    }

    // Verify that y^2 mod p == x^3 + 7 mod p
    Int check(y);
    check.ModSquareK1(&check);
    if (!check.IsEqual(&right))
    {
        printf("Debug GetY: y^2 != x^3 + 7 (mod p)\n");
    }
    else
    {
        printf("Debug GetY: y^2 == x^3 + 7 (mod p) - Valid solution found\n");
    }

    printf("Debug GetY: final y = %s\n", y.GetBase16());
    return y;
}

bool Secp256K1::EC(Point &p)
{
    Int right, left;

    // Calculate right side: x^3 + 7 (mod p)
    right.Set(&p.x);
    right.ModSquareK1(&right); // x^2
    right.ModMulK1(&p.x);      // x^3
    right.Add(&_7);
    right.Mod(&P); // (x^3 + 7) mod p

    // Calculate left side: y^2 (mod p)
    left.Set(&p.y);
    left.ModSquareK1(&left);
    left.Mod(&P);

    printf("Debug EC: x = %s\n", p.x.GetBase16());
    printf("Debug EC: y = %s\n", p.y.GetBase16());
    printf("Debug EC: x^3 + 7 mod p = %s\n", right.GetBase16());
    printf("Debug EC: y^2 mod p = %s\n", left.GetBase16());

    bool result = left.IsEqual(&right);
    printf("Debug EC: y^2 mod p %s x^3 + 7 mod p\n", result ? "==" : "!=");

    return result;
}

Point Secp256K1::ScalarMultiplication(Point &P, Int *scalar)
{
    Point R, Q;
    int bitLength = scalar->GetBitLength();
    R.Clear();
    R.z.SetInt32(1);
    if (!scalar->IsZero())
    {
        Q.Set(P);
        for (int i = 0; i < bitLength; i++)
        {
            if (scalar->GetBit(i))
            {
                R = Add(R, Q);
            }
            Q = Double(Q);
        }
    }
    R.Reduce();
    return R;
}

#define KEYBUFFCOMP(buff, p)                                                    \
    (buff)[0] = ((p).x.bits[7] >> 8) | ((uint32_t)(0x2 + (p).y.IsOdd()) << 24); \
    (buff)[1] = ((p).x.bits[6] >> 8) | ((p).x.bits[7] << 24);                   \
    (buff)[2] = ((p).x.bits[5] >> 8) | ((p).x.bits[6] << 24);                   \
    (buff)[3] = ((p).x.bits[4] >> 8) | ((p).x.bits[5] << 24);                   \
    (buff)[4] = ((p).x.bits[3] >> 8) | ((p).x.bits[4] << 24);                   \
    (buff)[5] = ((p).x.bits[2] >> 8) | ((p).x.bits[3] << 24);                   \
    (buff)[6] = ((p).x.bits[1] >> 8) | ((p).x.bits[2] << 24);                   \
    (buff)[7] = ((p).x.bits[0] >> 8) | ((p).x.bits[1] << 24);                   \
    (buff)[8] = 0x00800000 | ((p).x.bits[0] << 24);                             \
    (buff)[9] = 0;                                                              \
    (buff)[10] = 0;                                                             \
    (buff)[11] = 0;                                                             \
    (buff)[12] = 0;                                                             \
    (buff)[13] = 0;                                                             \
    (buff)[14] = 0;                                                             \
    (buff)[15] = 0x108;

#define KEYBUFFUNCOMP(buff, p)                                 \
    (buff)[0] = ((p).x.bits[7] >> 8) | 0x04000000;             \
    (buff)[1] = ((p).x.bits[6] >> 8) | ((p).x.bits[7] << 24);  \
    (buff)[2] = ((p).x.bits[5] >> 8) | ((p).x.bits[6] << 24);  \
    (buff)[3] = ((p).x.bits[4] >> 8) | ((p).x.bits[5] << 24);  \
    (buff)[4] = ((p).x.bits[3] >> 8) | ((p).x.bits[4] << 24);  \
    (buff)[5] = ((p).x.bits[2] >> 8) | ((p).x.bits[3] << 24);  \
    (buff)[6] = ((p).x.bits[1] >> 8) | ((p).x.bits[2] << 24);  \
    (buff)[7] = ((p).x.bits[0] >> 8) | ((p).x.bits[1] << 24);  \
    (buff)[8] = ((p).y.bits[7] >> 8) | ((p).x.bits[0] << 24);  \
    (buff)[9] = ((p).y.bits[6] >> 8) | ((p).y.bits[7] << 24);  \
    (buff)[10] = ((p).y.bits[5] >> 8) | ((p).y.bits[6] << 24); \
    (buff)[11] = ((p).y.bits[4] >> 8) | ((p).y.bits[5] << 24); \
    (buff)[12] = ((p).y.bits[3] >> 8) | ((p).y.bits[4] << 24); \
    (buff)[13] = ((p).y.bits[2] >> 8) | ((p).y.bits[3] << 24); \
    (buff)[14] = ((p).y.bits[1] >> 8) | ((p).y.bits[2] << 24); \
    (buff)[15] = ((p).y.bits[0] >> 8) | ((p).y.bits[1] << 24); \
    (buff)[16] = 0x00800000 | ((p).y.bits[0] << 24);           \
    (buff)[17] = 0;                                            \
    (buff)[18] = 0;                                            \
    (buff)[19] = 0;                                            \
    (buff)[20] = 0;                                            \
    (buff)[21] = 0;                                            \
    (buff)[22] = 0;                                            \
    (buff)[23] = 0;                                            \
    (buff)[24] = 0;                                            \
    (buff)[25] = 0;                                            \
    (buff)[26] = 0;                                            \
    (buff)[27] = 0;                                            \
    (buff)[28] = 0;                                            \
    (buff)[29] = 0;                                            \
    (buff)[30] = 0;                                            \
    (buff)[31] = 0x208;

#define KEYBUFFSCRIPT(buff, h)                                                                          \
    (buff)[0] = 0x00140000 | (uint32_t)h[0] << 8 | (uint32_t)h[1];                                      \
    (buff)[1] = (uint32_t)h[2] << 24 | (uint32_t)h[3] << 16 | (uint32_t)h[4] << 8 | (uint32_t)h[5];     \
    (buff)[2] = (uint32_t)h[6] << 24 | (uint32_t)h[7] << 16 | (uint32_t)h[8] << 8 | (uint32_t)h[9];     \
    (buff)[3] = (uint32_t)h[10] << 24 | (uint32_t)h[11] << 16 | (uint32_t)h[12] << 8 | (uint32_t)h[13]; \
    (buff)[4] = (uint32_t)h[14] << 24 | (uint32_t)h[15] << 16 | (uint32_t)h[16] << 8 | (uint32_t)h[17]; \
    (buff)[5] = (uint32_t)h[18] << 24 | (uint32_t)h[19] << 16 | 0x8000;                                 \
    (buff)[6] = 0;                                                                                      \
    (buff)[7] = 0;                                                                                      \
    (buff)[8] = 0;                                                                                      \
    (buff)[9] = 0;                                                                                      \
    (buff)[10] = 0;                                                                                     \
    (buff)[11] = 0;                                                                                     \
    (buff)[12] = 0;                                                                                     \
    (buff)[13] = 0;                                                                                     \
    (buff)[14] = 0;                                                                                     \
    (buff)[15] = 0xB0;

void Secp256K1::GetHash160(int type, bool compressed,
                           Point &k0, Point &k1, Point &k2, Point &k3,
                           uint8_t *h0, uint8_t *h1, uint8_t *h2, uint8_t *h3)
{

    ALIGN16 unsigned char sh0[64];
    ALIGN16 unsigned char sh1[64];
    ALIGN16 unsigned char sh2[64];
    ALIGN16 unsigned char sh3[64];

    switch (type)
    {
    case P2PKH:
    case BECH32:
    {
        if (!compressed)
        {
            uint32_t b0[32], b1[32], b2[32], b3[32];
            KEYBUFFUNCOMP(b0, k0);
            KEYBUFFUNCOMP(b1, k1);
            KEYBUFFUNCOMP(b2, k2);
            KEYBUFFUNCOMP(b3, k3);
            sha256sse_2B(b0, b1, b2, b3, sh0, sh1, sh2, sh3);
        }
        else
        {
            uint32_t b0[16], b1[16], b2[16], b3[16];
            KEYBUFFCOMP(b0, k0);
            KEYBUFFCOMP(b1, k1);
            KEYBUFFCOMP(b2, k2);
            KEYBUFFCOMP(b3, k3);
            sha256sse_1B(b0, b1, b2, b3, sh0, sh1, sh2, sh3);
        }
        ripemd160sse_32(sh0, sh1, sh2, sh3, h0, h1, h2, h3);
    }
    break;

    case P2SH:
    {
        unsigned char kh0[20], kh1[20], kh2[20], kh3[20];
        GetHash160(P2PKH, compressed, k0, k1, k2, k3, kh0, kh1, kh2, kh3);

        uint32_t b0[16], b1[16], b2[16], b3[16];
        KEYBUFFSCRIPT(b0, kh0);
        KEYBUFFSCRIPT(b1, kh1);
        KEYBUFFSCRIPT(b2, kh2);
        KEYBUFFSCRIPT(b3, kh3);

        sha256sse_1B(b0, b1, b2, b3, sh0, sh1, sh2, sh3);
        ripemd160sse_32(sh0, sh1, sh2, sh3, h0, h1, h2, h3);
    }
    break;
    }
}

void Secp256K1::GetHash160(int type, bool compressed, Point &pubKey, unsigned char *hash)
{
    unsigned char shapk[64];

    switch (type)
    {
    case P2PKH:
    case BECH32:
    {
        unsigned char publicKeyBytes[65];
        if (!compressed)
        {
            publicKeyBytes[0] = 0x4;
            pubKey.x.Get32Bytes(publicKeyBytes + 1);
            pubKey.y.Get32Bytes(publicKeyBytes + 33);
            sha256_65(publicKeyBytes, shapk);
        }
        else
        {
            publicKeyBytes[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
            pubKey.x.Get32Bytes(publicKeyBytes + 1);
            sha256_33(publicKeyBytes, shapk);
        }
        ripemd160_32(shapk, hash);
    }
    break;

    case P2SH:
    {
        unsigned char script[64];
        script[0] = 0x00; // OP_0
        script[1] = 0x14; // PUSH 20 bytes
        GetHash160(P2PKH, compressed, pubKey, script + 2);
        sha256(script, 22, shapk);
        ripemd160_32(shapk, hash);
    }
    break;
    }
}

#define KEYBUFFPREFIX(buff, k, fix)                          \
    (buff)[0] = (k->bits[7] >> 8) | ((uint32_t)(fix) << 24); \
    (buff)[1] = (k->bits[6] >> 8) | (k->bits[7] << 24);      \
    (buff)[2] = (k->bits[5] >> 8) | (k->bits[6] << 24);      \
    (buff)[3] = (k->bits[4] >> 8) | (k->bits[5] << 24);      \
    (buff)[4] = (k->bits[3] >> 8) | (k->bits[4] << 24);      \
    (buff)[5] = (k->bits[2] >> 8) | (k->bits[3] << 24);      \
    (buff)[6] = (k->bits[1] >> 8) | (k->bits[2] << 24);      \
    (buff)[7] = (k->bits[0] >> 8) | (k->bits[1] << 24);      \
    (buff)[8] = 0x00800000 | (k->bits[0] << 24);             \
    (buff)[9] = 0;                                           \
    (buff)[10] = 0;                                          \
    (buff)[11] = 0;                                          \
    (buff)[12] = 0;                                          \
    (buff)[13] = 0;                                          \
    (buff)[14] = 0;                                          \
    (buff)[15] = 0x108;

void Secp256K1::GetHash160_fromX(int type, unsigned char prefix,
                                 Int *k0, Int *k1, Int *k2, Int *k3,
                                 uint8_t *h0, uint8_t *h1, uint8_t *h2, uint8_t *h3)
{

    ALIGN16 unsigned char sh0[64];
    ALIGN16 unsigned char sh1[64];
    ALIGN16 unsigned char sh2[64];
    ALIGN16 unsigned char sh3[64];

    switch (type)
    {
    case P2PKH:
    {
        uint32_t b0[16], b1[16], b2[16], b3[16];
        KEYBUFFPREFIX(b0, k0, prefix);
        KEYBUFFPREFIX(b1, k1, prefix);
        KEYBUFFPREFIX(b2, k2, prefix);
        KEYBUFFPREFIX(b3, k3, prefix);

        sha256sse_1B(b0, b1, b2, b3, sh0, sh1, sh2, sh3);
        ripemd160sse_32(sh0, sh1, sh2, sh3, h0, h1, h2, h3);
    }
    break;

    case P2SH:
    {
        fprintf(stderr, "[E] Fixme unsupported case");
        exit(0);
    }
    break;
    }
}