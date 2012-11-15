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

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


#define clz_64(x) (__builtin_clzl(x))
#define clz_32(x) (__builtin_clzl(x)-32)
#define clz_16(x) (__builtin_clzl(x)-48)



typedef uint16_t float16;
typedef uint32_t float32;
typedef uint64_t float64;

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
#define _CLZ(width)		clz_ ## width
#define MANT_WIDTH(width)	_MANT_WIDTH(width)
#define UINT_FAST(width)	_UINT_FAST(width)
#define FLOAT(width)		_FLOAT(width)
#define CLZ(width)		_CLZ(width)

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



void raise_exception(int e)
{
}


bool float16_is_nan(float16 a)
{
    return (((a&0x7c00)==0x7c00) && ((a&0x03ff)!=0));
}


#define STOP 0x80011


#define WIDTH 16
#include "addsub.inc.c"
#undef WIDTH

#define WIDTH 32
#include "addsub.inc.c"
#undef WIDTH

#define WIDTH 64
#include "addsub.inc.c"
#undef WIDTH


#define SMALL 16
#define BIG 32
#include "convf2f.inc.c"
#undef SMALL
#undef BIG

#define SMALL 16
#define BIG 64
#include "convf2f.inc.c"
#undef SMALL
#undef BIG

#define SMALL 32
#define BIG 64
#include "convf2f.inc.c"
#undef SMALL
#undef BIG



#define _FUNC_MUL(W) float ## W ## _mul
#define FUNC_MUL(W) _FUNC_MUL(W)


FLOAT(WIDTH) FUNC_MUL(WIDTH) (FLOAT(WIDTH) a, FLOAT(WIDTH) b)
{
    EXP_TYPE(WIDTH) a_exp = EXTRACT_EXP(WIDTH, a);
    EXP_TYPE(WIDTH) b_exp = EXTRACT_EXP(WIDTH, b);
    UINT_FAST(WIDTH) a_mant = a & MANT_MASK(WIDTH);
    UINT_FAST(WIDTH) b_mant = b & MANT_MASK(WIDTH);

    if (a_exp==EXP_MAX(WIDTH)) {
/*
	if (a_mant==0) {
	    if (b_exp==EXP_MAX(WIDTH)) {
		if (b_mant==0)
		    return a ^ (b&SIGN_MASK(WIDTH)); // inf*inf = inf
		else
		    return b; // inf*nan = nan
	    } else if (b_exp==0)
		return NAN(WIDTH); // inf*0 = nan
	    else
		return a ^ (b&SIGN_MASK(WIDTH)); // inf*real = inf
	}
	else return NAN(WIDTH); // nan*x = nan
*/
	return ((a_mant!=0) || (b_exp==0 || (b&(SIGN_MASK(WIDTH)-1))>EXP_MASK_WIDTH))
	    ? NAN(WIDTH) // inf*(nan|0) = nan
	    : a ^ (b&SIGN_MASK(WIDTH)); // inf*(inf|real) = inf
    } else if (a_exp==0) {
	if (a_mant==0)
	    return  (b_exp==EXP_MAX(WIDTH)) 
		? NAN(WIDTH); // 0*(inf|nan) = nan
		: ((a ^ b) & SIGN_MASK(WIDTH)); // 0*(0|real) = 0

	// shift denorms into position and adjust exponent
	while (a_mant<=MANT_MASK(WIDTH)) {
	    a_mant <<= 1;
	    a_exp++;
	}
	// not yet complete ?
    }
	    
    // TODO: multiplication	
}




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


int64_t int64_from_float32_zero(float32 f)
{
    uint_fast16_t exp = (f>>MANT_WIDTH(32))&EXP_MAX(32);
    int_fast64_t abs;
    if (exp<BIAS(32))
	return 0;
    if (exp>BIAS(32)+63)
	return INT64_MIN; // too big or NaN

    if (exp<BIAS(32)+64-MANT_WIDTH(32))
	abs = (((f&MANT_MASK(32))|(MANT_MASK(32)+1))
	    << (exp-BIAS(32))) >> MANT_WIDTH(32);
    else
	abs = ((f&MANT_MASK(32))|(MANT_MASK(32)+1))
	    << (exp-BIAS(32)-MANT_WIDTH(32));
    return (f&SIGN_MASK(32)) ? -abs : abs;
}


// shift right the value v by s and round to nearest or even
#define SHR_NEAREST_EVEN(v, s) (((v) + (1L<<((s)-1))-1 + (((v)>>((s)))&1))>>(s))

float16 float16_from_int64(int64_t i)
{
    uint16_t sign=0;
    
    if (i<0)
    {
	if (i <= -(1L<<(EXP_MAX(16)-BIAS(16))))
	    return SIGN_MASK(16) | EXP_MASK(16); // - infinity
	sign = SIGN_MASK(16);
	i = -i;
    }
    else
    {
	if (i==0)
	    return 0;
	if (i >= (1L<<(EXP_MAX(16)-BIAS(16))))
    	    return EXP_MASK(16); // + infinity
    }

    uint_fast16_t shift = clz_16(i);
    // highest representable value is smaller than 2^15, therefore
    // clz_16() is sufficient, other cases are catched in the if-clause
    // and return infinity
    uint_fast16_t mant = (i<<shift)&0x7fff;

    // rounding
    mant = SHR_NEAREST_EVEN(mant, 16-1-MANT_WIDTH(16));

    return sign | (((BIAS(16)+16-1-shift)<<MANT_WIDTH(16)) + mant);
    // the + instead of a | guarantees that in the case of an overflow
    // by the rounding the exponent is increased by 1
}


float32 float32_from_int64(int64_t i)
{
    float32 sign=0;
    
    if (i<0)
    {
	sign = SIGN_MASK(32);
	i = -i;
    }
    else if (i==0)
	return 0;

    uint_fast64_t shift = clz_64(i);
    uint_fast64_t mant = (i<<shift)&0x7fffffffffffffff;

    // rounding
    mant = SHR_NEAREST_EVEN(mant, 64-1-MANT_WIDTH(32));

    return sign | (((BIAS(32)+64-1-shift)<<MANT_WIDTH(32)) + mant);
    // the + instead of a | guarantees that in the case of an overflow
    // by the rounding the exponent is increased by 1
}






typedef union {
    float32 u;
    float f;
} uf32_t;

typedef union {
    float64 u;
    double f;
} uf64_t;

float float_from_float16(float16 f)
{
    uf32_t r;
    r.u = float32_from_float16(f);
    return r.f;
}

float16 float16_from_float(float f)
{
    uf32_t a;
    a.f = f;
    return float16_from_float32(a.u);
}

unsigned hex_from_float(float f)
{
    uf32_t r;
    r.f = f;
    return r.u;
}

double double_from_float16(float16 f)
{
    uf64_t r;
    r.u = float64_from_float16(f);
    return r.f;
}


float16 float16_from_double(double f)
{
    uf64_t a;
    a.f = f;
    return float16_from_float64(a.u);
}

float32 float32_from_float(float f)
{
    uf32_t r;
    r.f = f;
    return r.u;
}

float float_from_float32(float32 f)
{
    uf32_t r;
    r.u = f;
    return r.f;
}

float64 float64_from_double(double f)
{
    uf64_t r;
    r.f = f;
    return r.u;
}

double double_from_float64(float64 f)
{
    uf64_t r;
    r.u = f;
    return r.f;
}


uint64_t rand64_seed = 0;

uint64_t rand64()
{
    rand64_seed = 6364136223846793005L*rand64_seed + 1442695040888963407L;
    return rand64_seed;
}

double rand_double()
{
    uf64_t	a;
    a.u = rand64();
    return a.f;
}




int main()
{
    unsigned i, j;
    float16 a, b;
    //float16 b, correct, test;
    int64_t i64;
    float32 a32, b32;
    float a32f, b32f;
    float64 a64, b64;
//    double a64f, b64f;



    for (i=0; i<0x10000; i++)
    {
	a = i;
	
	// float16 <-> int64

	i64 = int64_from_float16_zero(a);
	if (i64 != (int64_t)double_from_float16(a))
	{
	    printf("int64_from_float16(%04x %g) = %ld != %ld\n",
		a, float_from_float16(a), 
		i64, (int64_t)double_from_float16(a));
	    return 1;
	}

	if (float16_from_int64(i64) != float16_from_double((double)i64))
	{
	    printf("float16_from_int64(%ld) = %04x %g != %04x %g\n",
		i64, float16_from_int64(i64), double_from_float16(float16_from_int64(i64)),
		float16_from_double((double)i64),
		double_from_float16(float16_from_double((double)i64)));
	    return 1;
	}
	
	i64 = a ^ (a<<3);
	if (float16_from_int64(i64) != float16_from_double((double)i64))
	{
	    printf("float16_from_int64(%ld) = %04x %g != %04x %g\n",
		i64, float16_from_int64(i64), double_from_float16(float16_from_int64(i64)),
		float16_from_double((double)i64),
		double_from_float16(float16_from_double((double)i64)));
	    return 1;
	}

	// float32 <-> int64	

	i64 = ((float)rand64());
	if (float32_from_int64(i64) != float32_from_float((float)i64))
	{
	    printf("float32_from_int64(%ld) = %08x %g != %08x %g\n",
		i64, float32_from_int64(i64), float_from_float32(float32_from_int64(i64)),
		float32_from_float((float)i64),
		(float)i64);
	    return 1;
	}

	// float16 <-> float64

	b = float16_from_float64(float64_from_float16(a));
	if (a!=b && b!=NOTANUM(16))
        {
	    printf("%04x (%g float:%g) != %04x (%g)\n",
		    a, double_from_float16(a), float_from_float16(a),
		    b, double_from_float16(b));
	    return 1;
	}
	
	// float32 <-> float64
	
	a64 = rand64();
	a32 = a64;
	a32f = float_from_float32(a32);
	a64 = float64_from_float32(a32);
//	if (float32_from_float(double_from_float64(a64)) != float32_from_float64(a64))
	if (!(isnan(a32f) && isnan(double_from_float64(a64))) &&
	    float64_from_double(a32f) != a64)
        {
	    printf("float64_from_float32(%08x %g) = %016lx %g != %016lx %g\n",
		a32, a32f,
		float64_from_float32(a32), double_from_float64(float64_from_float32(a32)),
		float64_from_double(a32f), float_from_float32(a32));
	    return 1;
	}
	

	for (j=0; j<0x10000; j++)
	{
	    b=j;
	    float16 test = float16_add(a, b);

	    float sum = float_from_float16(a) + float_from_float16(b);
	    float16 correct = float16_from_float(sum);
	    if (correct!=test) // works only, as NaN is always default value
	    {
		printf("%04x (%g %08x) + %04x (%g %08x) = "
		    "%04x (%g %08x exact: %g %08x) but: %04x (%g %08x)\n",
		    a, float_from_float16(a), hex_from_float(float_from_float16(a)),
		    b, float_from_float16(b), hex_from_float(float_from_float16(b)),
		    correct, float_from_float16(correct), hex_from_float(float_from_float16(correct)),
		    sum, hex_from_float(sum),
		    test, float_from_float16(test), hex_from_float(float_from_float16(test)));
		return 1;
	    }

	    double sumd = double_from_float16(a) + double_from_float16(b);
	    correct = float16_from_double(sumd);
	    if (correct!=test) // works only, as NaN is always default value
	    {
		printf("%04x (%g %016lx) + %04x (%g %016lx) = "
		    "%04x (%g %016lx exact: %g) but: %04x (%g)\n",
		    a, double_from_float16(a), float64_from_float16(a),
		    b, double_from_float16(b), float64_from_float16(b),
		    correct, float_from_float16(correct), float64_from_float16(correct),
		    sumd, //hex_from_float(sumd),
		    test, float_from_float16(test));
		return 1;
	    }
	}
    }




#define MAX_ITER 1000000000



    uint64_t *pc, *pd;
    double dc, dd;
    
    pc = (uint64_t *)&dc;
    pd = (uint64_t *)&dd;
    for (i=0; i<MAX_ITER; i++)
    {
	a64 = float64_from_double(rand_double());
	b64 = float64_from_double(rand_double());
	
	dc = double_from_float64(a64)+double_from_float64(b64);
	dd = double_from_float64(float64_add(a64, b64));

	if ((*pc != *pd) && !(isnan(dc) && isnan(dd)))
	    printf("%d: %016lx (%g) + %016lx (%g) = %016lx (%g) but: %016lx (%g)\n",
	    i, a64, double_from_float64(a64), b64, double_from_float64(b64), *pc, dc, *pd, dd);
/*
	else
	    printf("%d: %g + %g = %g\n",
	    iter, a, b, c);
*/
    }



    return 0;
}


