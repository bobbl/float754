/* Half precision floating point emulation
 * Licenced under the ISC license (similar to the MIT/Expat license)
 *
 * Copyright (c) 2011-2012 Jörg Mische <bobbl@gmx.de>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
 * OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


#define clz_64(x) (__builtin_clzl(x))
#define clz_32(x) (__builtin_clzl(x)-32)
#define clz_16(x) (__builtin_clzl(x)-48)



typedef uint16_t float16;
typedef float float32;
typedef double float64;

typedef union {
    uint16_t u;
    float16 f;
} uf16_t;

typedef union {
    uint32_t u;
    float f;
} uf32_t;

typedef union {
    uint64_t u;
    double f;
} uf64_t;



#define EXC_INVALID_OP		1
#define EXC_INEXACT		2


//#define F16_SIGN_MASK		0x8000
//#define F16_EXP_MASK		0x7c00
//#define F16_MANT_MASK		0x03ff
//#define F16_SIGNALLING_MASK	0x0200
//#define F16_NAN			0xfe00
#define F16_MANT_WIDTH		10
//#define F16_EXP_MAX		0x1f
//#define F16_BIAS		15L

//#define F32_SIGN_MASK		0x80000000
//#define F32_EXP_MASK		0x7f800000
//#define F32_MANT_MASK		0x007fffff
//#define F32_SIGNALLING_MASK	0x00400000
//#define F32_NAN			0xffc00000
#define F32_MANT_WIDTH		23
//#define F32_EXP_MAX		0xff
//#define F32_BIAS		127L

//#define F64_SIGN_MASK		0x8000000000000000
//#define F64_EXP_MASK		0x7ff0000000000000
//#define F64_MANT_MASK		0x000fffffffffffff
//#define F64_SIGNALLING_MASK	0x0008000000000000
//#define F64_NAN			0xfff8000000000000
#define F64_MANT_WIDTH		52
//#define F64_EXP_MAX		0x7ff
//#define F64_BIAS		1023L

#define F128_MANT_WIDTH		112
//#define F128_BIAS		16383L

#define _MANT_WIDTH(width)	F ## width ## _MANT_WIDTH
#define _UINT_FAST(width)	uint_fast ## width ## _t
#define _FLOAT(width)		float ## width
#define _UF(width)		uf ## width ## _t
#define _CLZ(width)		clz_ ## width
#define _UMUL_PPMM(width)	UMUL_PPMM_ ## width
#define MANT_WIDTH(width)	_MANT_WIDTH(width)
#define UINT_FAST(width)	_UINT_FAST(width)
#define FLOAT(width)		_FLOAT(width)
#define UF(width)		_UF(width)
#define CLZ(width)		_CLZ(width)
#define UMUL_PPMM(width)	_UMUL_PPMM(width)

#define SIGN_MASK(width)	(1L<<(width-1))
#define EXP_MASK(width)		((1L<<(width-1))-(1L<<MANT_WIDTH(width)))
#define MANT_MASK(width)	((1L<<MANT_WIDTH(width))-1)
#define SIGNALLING_MASK(width)	(1L<<(MANT_WIDTH(width)-1))
#define NOTANUM(width)		((1L<<(width-1))|(((1L<<(width-MANT_WIDTH(width)))-1)<<(MANT_WIDTH(width)-1)))
#define EXP_MAX(width)		((1L<<(width-1-MANT_WIDTH(width)))-1)
#define BIAS(width)		((1L<<(width-2-MANT_WIDTH(width)))-1)
#define EXTRACT_EXP(width, f)	(((f)>>MANT_WIDTH(width))&EXP_MAX(width))
#define EXTRACT_MANT(width, f)	((f)&MANT_MASK(width))

#define EXP_TYPE(width)		int16_t




#define UMUL_PPMM_16(h, l, u, v) 					\
  do {									\
    uint_fast32_t p = (uint_fast32_t)(u) * (uint_fast32_t)(v);		\
    l = p & 0xffff;							\
    h = p >> 16;							\
  } while (0);

#define UMUL_PPMM_32(h, l, u, v) 					\
  do {									\
    uint_fast64_t p = (uint_fast64_t)(u) * (uint_fast64_t)(v);		\
    l = p & 0xffffffff;							\
    h = p >> 32;							\
  } while (0);

#define UMUL_PPMM_64(h, l, u, v) 					\
  do {									\
    uint_fast64_t ul = u & 0xffffffff;					\
    uint_fast64_t uh = u >> 32;						\
    uint_fast64_t vl = v & 0xffffffff;					\
    uint_fast64_t vh = v >> 32;						\
    uint_fast64_t plh = ul*vh;						\
    uint_fast64_t phl = uh*vl;						\
    plh += phl;								\
    uint_fast64_t phh = uh*vh + ((plh<phl) ? (uint_fast64_t)1<<32 : 0);	\
    uint_fast64_t carry = (plh<<32);					\
    uint_fast64_t pll = ul*vl + carry;					\
    h = phh + (plh>>32) + ((pll<carry) ? 1 : 0);			\
    l = pll;								\
  } while (0);



#define ROUND_1(shifted, remains) \
    (((shifted&1) == 0) \
	? (shifted>>1) \
	: ((remains == 0) \
	    ? (((shifted>>1)+1) & (~1)) \
	    : ((shifted>>1) + 1)))

#define ROUND_2(shifted, remains) \
    ((((shifted&1)!=0) && remains == 0) \
	? (((shifted>>1) + 1) & (~1)) \
	: ((shifted>>1) + (shifted&1)))

#define ROUND_3(shifted, remains) \
    (((shifted>>1) + (shifted&1)) & ((remains==0 && ((shifted&1)!=0)) ? (~1) : (~0)))

#define ROUND_4(shifted, remains) \
    (((shifted&2)==0) \
	? ((shifted>>1) + ((remains!=0)&(shifted&1))) \
	: ((shifted>>1) + (shifted&1)))

#define ROUND_5(shifted, remains) \
    ((shifted>>1) +  \
	(((shifted&2)==0) \
	    ? ((remains!=0)&(shifted&1)) \
	    : (shifted&1)))

#define ROUND_6(shifted, remains) \
    ((shifted>>1) + \
	(shifted & \
	    (((shifted&2)==0) \
		? (remains!=0) \
		: 1)))

#define ROUND_7(shifted, remains) \
    ((shifted>>1) + ((shifted&1) & ((shifted>>1) | (remains!=0))))

#define ROUND_8(shifted, remains) \
    ((shifted>>1) + ((remains!=0) \
	? (shifted&1) \
	: (shifted & (shifted>>1) & 1)))

// probably the fastest alternative
#define ROUND_9(shifted, remains) \
    ((shifted>>1) + (shifted & ((remains==0) ? ((shifted>>1)&1) : 1)))

#define ROUND ROUND_9



void raise_exception(int e)
{
}


bool float16_is_nan(float16 a)
{
    return (((a&0x7c00)==0x7c00) && ((a&0x03ff)!=0));
}




#define PLATFORM_WIDTH 64


#define WIDTH 16
#include "operations.inc.c"
#undef WIDTH

#define WIDTH 32
#include "operations.inc.c"
#undef WIDTH

#define WIDTH 64
#include "operations.inc.c"
#undef WIDTH


#define SMALL 16
#define BIG 32
#include "convert.inc.c"
#undef SMALL
#undef BIG

#define SMALL 16
#define BIG 64
#include "convert.inc.c"
#undef SMALL
#undef BIG

#define SMALL 32
#define BIG 64
#include "convert.inc.c"
#undef SMALL
#undef BIG





int64_t int64_from_float16_zero(float16 f)
{
    uint_fast16_t exp = (f>>MANT_WIDTH(16))&EXP_MAX(16);

    if (exp<BIAS(16))
	return 0;
    if (exp==EXP_MAX(16))
	return INT64_MIN; // too big or NaN

    int_fast32_t abs = (((f&MANT_MASK(16))|(MANT_MASK(16)+1))
	<< (exp-BIAS(16))) >> MANT_WIDTH(16);
    return (f&SIGN_MASK(16)) ? -abs : abs;
}


int64_t int64_from_float32_zero(float32 ff)
{
    uf32_t f;
    f.f = ff;
    uint_fast16_t exp = (f.u>>MANT_WIDTH(32))&EXP_MAX(32);
    int_fast64_t abs;
    if (exp<BIAS(32))
	return 0;
    if (exp>BIAS(32)+63)
	return INT64_MIN; // too big or NaN

    if (exp<BIAS(32)+64-MANT_WIDTH(32))
	abs = (((f.u&MANT_MASK(32))|(MANT_MASK(32)+1))
	    << (exp-BIAS(32))) >> MANT_WIDTH(32);
    else
	abs = ((f.u&MANT_MASK(32))|(MANT_MASK(32)+1))
	    << (exp-BIAS(32)-MANT_WIDTH(32));
    return (f.u&SIGN_MASK(32)) ? -abs : abs;
}


// shift right the value v by s and round to nearest or even
#define SHR_NEAREST_EVEN(v, s) (((v) + (1L<<((s)-1))-1 + (((v)>>((s)))&1))>>(s))

float16 float16_from_int64(int64_t i)
{
    uf16_t r; r.u=0;
    
    if (i<0) {
	if (i <= -(1L<<(EXP_MAX(16)-BIAS(16)))) {
	    r.u = SIGN_MASK(16) | EXP_MASK(16); // - infinity
	    return r.f;
	}
	r.u = SIGN_MASK(16);
	i = -i;
    } else if (i==0) {
	return r.f;
    } else if (i >= (1L<<(EXP_MAX(16)-BIAS(16)))) {
    	r.u = EXP_MASK(16); // + infinity
    	return r.f;
    }

    uint_fast16_t shift = clz_16(i);
	// highest representable value is smaller than 2^15, therefore
	// clz_16() is sufficient, other cases are catched in the if-clause
	// and return infinity

    r.u |= ((BIAS(16)+16-1-shift)<<MANT_WIDTH(16))
	+ SHR_NEAREST_EVEN((i<<shift)&0x7fff, 16-1-MANT_WIDTH(16));
    return r.f;
	// the + instead of a | guarantees that in the case of an overflow
	// by the rounding the exponent is increased by 1
}


float32 float32_from_int64(int64_t i)
{
    uf32_t r; r.u = 0;

    if (i<0) {
	r.u = SIGN_MASK(32);
	i = -i;
    } 
    if (i!=0) {
	uint_fast64_t shift = clz_64(i);
	r.u |= ((BIAS(32)+64-1-shift)<<MANT_WIDTH(32)) 
	    + SHR_NEAREST_EVEN((i<<shift)&0x7fffffffffffffff, 64-1-MANT_WIDTH(32));
	// the + instead of a | guarantees that in the case of an overflow
	// by the rounding the exponent is increased by 1
    }
    return r.f;
}






uint64_t rand64_seed = 0;

uint64_t rand64()
{
    rand64_seed = 6364136223846793005L*rand64_seed + 1442695040888963407L;
    return rand64_seed;
}


int main()
{
    unsigned i, j;
    uf16_t a16, b16, c16;
    uf32_t a32, b32, c32;
    uf64_t a64, b64, c64, d64;
    int64_t i64;
    
    for (i=0; i<0x10000; i++)
    {
	a16.u = i;
	
	// float16 <-> int64

	i64 = int64_from_float16_zero(a16.f);
	if (i64 != (int64_t)float64_from_float16(a16.f))
	{
	    printf("int64_from_float16(%04x %g) = %ld != %ld\n",
		a16.u, float64_from_float16(a16.f), 
		i64, (int64_t)float64_from_float16(a16.f));
	    return 1;
	}

	b16.f = float16_from_int64(i64);
	c16.f = float16_from_float64((double)i64);
	if (b16.u != c16.u)
	{
	    printf("float16_from_int64(%ld) = %04x %g != %04x %g\n",
		i64, b16.u, float64_from_float16(b16.f),
		     c16.u, float64_from_float16(c16.f));
	    return 1;
	}
	
	i64 = i ^ (i<<3);
	b16.f = float16_from_int64(i64);
	c16.f = float16_from_float64((double)i64);
	if (b16.u != c16.u)
	{
	    printf("float16_from_int64(%ld) = %04x %g != %04x %g\n",
		i64, b16.u, float64_from_float16(b16.f),
		     c16.u, float64_from_float16(c16.f));
	    return 1;
	}

	// float32 <-> int64	

	i64 = ((float)rand64());
	b32.f = float32_from_int64(i64);
	c32.f = (float)i64;
	if (b32.u != c32.u)
	{
	    printf("float16_from_int64(%ld) = %04x %g != %04x %g\n",
		i64, b32.u, b16.f, c32.u, c32.f);
	    return 1;
	}

	// float16 <-> float64

	a64.f = float64_from_float16(a16.f);
	b16.f = float16_from_float64(a64.f);
	if (a16.u!=b16.u && b16.u!=NOTANUM(16))
        {
	    printf("float64_from_float16 %04x (%016lx %g) != %04x (%g)\n",
		    a16.u, a64.u, float64_from_float16(a16.f),
		    b16.u, float64_from_float16(b16.f));
	    return 1;
	}
	
	// float32 <-> float64
	
	a64.u = rand64();
	a32.f = a64.f;
	a64.f = float64_from_float32(a32.f);
	b64.f = a32.f;
//	if (float32_from_float(double_from_float64(a64)) != float32_from_float64(a64))
	if (!(isnan(a32.f) && isnan(a64.f)) && a64.u != b64.u)
        {
	    printf("float64_from_float32(%08x %g) = %016lx %g != %016lx %g\n",
		a32.u, a32.f, a64.u, a64.f, b64.u, b64.f);
	    return 1;
	}
	

#define CORRECT_OP(x, y) 	(x * y)
#define NEW_OP(x, y) 		float16_mul(x, y)
//#define CORRECT_OP(x, y) 	(x + y)
//#define NEW_OP(x, y) 		float16_add(x, y)
	for (j=0; j<0x10000; j++)
	{
	    b16.u = j;
	    uf16_t test, correct;
	    test.f = NEW_OP(a16.f, b16.f);

/*
	    double sum = CORRECT_OP(double_from_float16(a), double_from_float16(b));
	    correct = float16_from_float64(sum);
	    if (correct!=test.u) // works only, because NaN is always default value
	    {
		printf("%04x (%g %08x) ° %04x (%g %08x) = "
		    "%04x (%g %08x exact: %g %08x) but: %04x (%g %08x)\n",
		    a, float_from_float16(a), hex_from_float(float_from_float16(a)),
		    b, float_from_float16(b), hex_from_float(float_from_float16(b)),
		    correct, float_from_float16(correct), hex_from_float(float_from_float16(correct)),
		    sum, hex_from_float(sum),
		    test, float_from_float16(test), hex_from_float(float_from_float16(test)));
		return 1;
	    }
*/

	    double sumd = CORRECT_OP(float64_from_float16(a16.f), float64_from_float16(b16.f));
	    c16.f = float16_from_float64(sumd);
	    if (c16.u!=test.u) // works only, because NaN is always default value
	    {
		a64.f = float64_from_float16(a16.f);
		b64.f = float64_from_float16(b16.f);
		c64.f = float64_from_float16(c16.f);
		printf("%04x (%g %016lx) ° %04x (%g %016lx) = "
		    "%04x (%g %016lx exact: %g) but: %04x (%g)\n",
		    a16.u, a64.f, a64.u,
		    b16.u, b64.f, b64.u,
		    c16.u, c64.f, c64.u,
		    sumd, //hex_from_float(sumd),
		    test.u, float64_from_float16(test.f));
		return 1;
	    }
	}
    }




//#define MAX_ITER 1000000000
#define MAX_ITER 100


    for (i=0; i<MAX_ITER; i++)
    {
	a64.u = rand64();
	b64.u = rand64();
	c64.f = a64.f + b64.f;
	d64.f = float64_add(a64.f, b64.f);
	if ((c64.u != d64.u) && !(isnan(c64.f) && isnan(d64.f)))
	    printf("%d: %016lx (%g) + %016lx (%g) = %016lx (%g) but: %016lx (%g)\n",
	    i, a64.u, a64.f, b64.u, b64.f, c64.u, c64.f, d64.u, d64.f);
    }



    return 0;
}


