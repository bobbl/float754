/* Half precision floating point emulation
 * Licenced under the ISC license (similar to the MIT/Expat license)
 *
 * Copyright (c) 2011 Jörg Mische <bobbl@gmx.de>
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
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

    unsigned i, j;

/*

uint_fast16_t clz_16(uint_fast16_t x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);

    // count ones (population count)
    x -= ((x >> 1) & 0x5555);
    x = ((x >> 2) & 0x3333) + (x & 0x3333);
    x = ((x >> 4) + x) & 0x0f0f;
    x += (x >> 8);

    return (16 - (x & 0x7f));
}

uint_fast16_t clz_32(uint_fast32_t x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    // count ones (population count)
    x -= ((x >> 1) & 0x55555555);
    x = ((x >> 2) & 0x33333333) + (x & 0x33333333);
    x = ((x >> 4) + x) & 0x0f0f0f0f;
    x += (x >> 8);
    x += (x >> 16);

    return (32 - (x & 0x7f));
}

uint_fast16_t clz_64(uint_fast64_t x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);

    // count ones (population count)
    x -= ((x >> 1) & 0x5555555555555555);
    x = ((x >> 2) & 0x3333333333333333) + (x & 0x3333333333333333);
    x = ((x >> 4) + x) & 0x0f0f0f0f0f0f0f0f;
    x += (x >> 8);
    x += (x >> 16);
    x += (x >> 32);

    return (64 - (x & 0x7f));
}

uint_fast16_t _clz_64(uint_fast64_t x)
{

    uint_fast8_t count = 0x40;
    uint_fast64_t b;
    b = (x<0x100000000) ? 0 : 32; count -= b; x >>= b;
    b = (x<0x10000) ? 0 : 16; count -= b; x >>= b;
    b = (x<0x100) ? 0 : 8; count -= b; x >>= b;
    b = (x<0x10) ? 0 : 4; count -= b; x >>= b;
    b = (x<0x4) ? 0 : 2; count -= b; x >>= b;
    return count - (x<0x2 ? x : 2);
}


*/
#define clz_64(x) (__builtin_clzl(x))
#define clz_32(x) (__builtin_clzl(x)-32)
#define clz_16(x) (__builtin_clzl(x)-48)



typedef uint16_t float16;
typedef uint32_t float32;
typedef uint64_t float64;

#define EXC_INVALID_OP	1
#define EXC_INEXACT		2


#define F16_SIGN_MASK		0x8000
#define F16_EXP_MASK		0x7c00
#define F16_MANT_MASK		0x03ff
#define F16_SIGNALLING_MASK	0x0200
#define F16_NAN			0xfe00
#define F16_MANT_WIDTH		10
#define F16_SPECIAL_EXP		0x1f
#define F16_BIAS		15
#define F16_EXTRACT_EXP(f)	(((f)>>F16_MANT_WIDTH)&F16_SPECIAL_EXP)
#define F16_EXTRACT_MANT(f)	((f)&F16_MANT_MASK)

#define F32_SIGN_MASK		0x80000000
#define F32_EXP_MASK		0x7f800000
#define F32_MANT_MASK		0x007fffff
#define F32_SIGNALLING_MASK	0x00400000
#define F32_NAN			0xffc00000
#define F32_MANT_WIDTH		23
#define F32_SPECIAL_EXP		0xff
#define F32_BIAS		127
#define F32_EXTRACT_EXP(f)	(((f)>>F32_MANT_WIDTH)&F32_SPECIAL_EXP)
#define F32_EXTRACT_MANT(f)	((f)&F32_MANT_MASK)

#define F64_SIGN_MASK		0x8000000000000000
#define F64_EXP_MASK		0x7ff0000000000000
#define F64_MANT_MASK		0x000fffffffffffff
#define F64_SIGNALLING_MASK	0x0008000000000000
#define F64_NAN			0xfff8000000000000
#define F64_MANT_WIDTH		52
#define F64_SPECIAL_EXP		0x7ff
#define F64_BIAS		1023
#define F64_EXTRACT_EXP(f)	(((f)>>F64_MANT_WIDTH)&F64_SPECIAL_EXP)
#define F64_EXTRACT_MANT(f)	((f)&F64_MANT_MASK)



void raise_exception(int e)
{
}


bool float16_is_nan(float16 a)
{
    return (((a&0x7c00)==0x7c00) && ((a&0x03ff)!=0));
}


#define STOP 0x80011

float16 float16_add(float16 a, float16 b)
{
    int_fast8_t lo_exp = F16_EXTRACT_EXP(a);
    int_fast8_t hi_exp = F16_EXTRACT_EXP(b);
    int_fast16_t diff_exp = lo_exp - hi_exp;
    float16 hi_sign;
    int_fast16_t sum;

    if (diff_exp == 0)
    {
	if (hi_exp == F16_SPECIAL_EXP)
	{
	    // infinity or NaN
	    if ( ((a|b) & F16_MANT_MASK) == 0)
	    {
		// two times infinity (mantissa==0)
		if (a != b) // already equal except the sign
		{
		    // infinity minus infinity raises an error
		    raise_exception(EXC_INVALID_OP);
		    return F16_NAN;
		}
		else
		    return a;
	    }
	    if ( ((a&b) & F16_SIGNALLING_MASK) == 0)
		// one of the two NaNs is signalling 
		// (highest bit of the mantissa==0)
		raise_exception(EXC_INVALID_OP);
	    return F16_NAN;
	}

	// same exponent
	uint64_t a_mant = F16_EXTRACT_MANT(a);
	uint64_t b_mant = F16_EXTRACT_MANT(b);
	hi_sign = a & F16_SIGN_MASK;

	if (((a ^ b) & F16_SIGN_MASK) == 0)
	{
	    // both numbers have the same sign -> addition
	    sum = a_mant + b_mant;

	    if (hi_exp==(F16_SPECIAL_EXP-1))
		// by the addition the exponent is increased by 1
		// and then it is too big 
		return hi_sign | F16_EXP_MASK;
	    if (hi_exp == 0)
		// both numbers are subnormal. If there is an overflow, i.e. the
		// sum is normal, bit 52 of the result is set, which means that
		// the exponent is 1 and not 0. Hence the conversion from subnormal
		// to normal works automatically, only the sign has to be added.
		return hi_sign | sum;

	    sum = (sum + ((sum>>1)&1)) >> 1; // round to nearest or even
	    hi_exp++;
		// because of the two implicit ones the result has one digit
		// to much, i.e. the mantissa has to be shifted right by 1
		// and the exponent must be increased by 1
	}
	else
	{
	    // numbers have different signs -> subtraction
	    sum = a_mant - b_mant;

	    if (sum == 0)
	        return 0;
	    
	    if (sum < 0)
	    {
	        hi_sign ^= F16_SIGN_MASK;
		sum = -sum;
	    }

	    // very easy for subnormal numbers
	    if (hi_exp == 0)
		return hi_sign | sum;

	    // normalise mantissa
	    uint_fast8_t zeros = clz_16(sum) - 5;
	    sum <<= zeros;
    	    hi_exp = hi_exp - 1 - zeros;

	    // subnormal result
	    if (hi_exp < 0)
	        return hi_sign | (sum >> (-hi_exp));
	}
    }
    else 
    {
	uint64_t lo_mant, hi_mant;
	if (diff_exp > 0)
	{
	    hi_exp += diff_exp;
	    lo_exp -= diff_exp;
		// exchange hi_exp and lo_exp
	    lo_mant = F16_EXTRACT_MANT(b);
	    hi_mant = F16_EXTRACT_MANT(a);
	    hi_sign = a;

	}
	else // (diff_exp < 0)
	{
	    diff_exp = -diff_exp;
	    lo_mant = F16_EXTRACT_MANT(a);
	    hi_mant = F16_EXTRACT_MANT(b);
	    hi_sign = b;
	}
	// temporarily the lower bits of hi_sign are not masked out,
	// in order to return the correct value in the case of an exception

	if (hi_exp == F16_SPECIAL_EXP)
	{
	    if (hi_mant==0)
		return hi_sign;				// infinity
	    if ((hi_mant & 0x0200)==0)
	        raise_exception(EXC_INVALID_OP);        // signalling NaN
	    return F16_NAN;			// NaN
	}

	if (diff_exp > F16_MANT_WIDTH+2)		// lo much too low
	{
	    if (lo_exp!=0 || lo_mant!=0)
		raise_exception(EXC_INEXACT);
    	    return hi_sign;
	}

	// from now on, only the sign of hi_sign is needed
	hi_sign &= F16_SIGN_MASK;

	if (lo_exp == 0)
	    // lo subnormal
    	    diff_exp--;
	else
	    // insert implicit 1
    	    lo_mant |= 1<<F16_MANT_WIDTH;


#define REM_BITS (F16_MANT_WIDTH+2)
#define REM_MASK ((1<<REM_BITS)-1)
#define REM_HALF (1<<(REM_BITS-1))
	uint64_t rem = (lo_mant << (REM_BITS-diff_exp)) & REM_MASK;

	lo_mant >>= diff_exp;
	hi_mant |= 1<<F16_MANT_WIDTH;

	if (((a ^ b) & F16_SIGN_MASK) == 0)
	{
	    // same sign
	    // simplifies the normalisation, but overflow checks are needed
	    sum = hi_mant + lo_mant;

	    if (sum < (2<<F16_MANT_WIDTH))		// no bit overflow
	    {
		if (rem == REM_HALF)
		    sum = (sum+1) & (~1);    		// round to even
		else if (rem > REM_HALF)
		    sum++;				// round up
		// else round down (nothing to do)

		hi_exp--; 
	    }
	    else					// 1 bit overflow
	    {
		if (hi_exp>=F16_SPECIAL_EXP-1)
		    return hi_sign | F16_EXP_MASK;	// infinity

		if ((sum&1) == 0)
		    sum = sum>>1;			// round down
		else if (rem == 0)
		    sum = ((sum>>1)+1) & (~1);    	// round to even
		else
		    sum = (sum>>1) + 1;			// round up

// alternative formulation:
//		if (((sum&1)!=0) && rem == 0)
//		    sum = ((sum>>1) + 1) & (~1);	// round up to even
//		else
//		    sum = (sum>>1) + (sum&1);		// round to nearest

// another alternative formulation:
//		sum = ((sum>>1) + (sum&1)) & ((rem==0 && ((sum&1)!=0)) ? (~1) : (~0));
	    }
	}
	else
	{
	    // not the same sign
	    sum = hi_mant - lo_mant;

	    if (rem != 0)
		sum--;

	    uint_fast8_t zeros = clz_16(sum) - 5;

	    sum = (sum << zeros) | (((-rem) & REM_MASK) >> (REM_BITS-zeros));
	    rem = (rem << zeros) & REM_MASK;
	    if (rem == REM_HALF)
    		sum = (sum+1) & (~1);    		// round to even
	    else if (rem != 0 && rem < REM_HALF)
		sum++;

	    hi_exp = hi_exp - 1 - zeros;
	    if (hi_exp < 0)
		return hi_sign | (sum >> (-hi_exp));
	}
    }
    // Now the leading 1 must be masked out. But it is more efficient to
    // decrement the exponent by 1 and then add the implicit 1.
    // The decrementation as already been done earlier.
    return hi_sign | (((uint64_t)hi_exp<<F16_MANT_WIDTH) + sum);
}



float32 float32_from_float16(float16 f)
{
    uint_fast32_t sign = ((uint_fast32_t)f&F16_SIGN_MASK)<<16;
    uint_fast32_t abs = f&(F16_EXP_MASK|F16_MANT_MASK);
    
    if (abs<=F16_MANT_MASK)
    {
	if (abs==0) return sign;		// zero

	uint_fast32_t shift = clz_16(abs);	// subnormal
	return sign | (((F32_BIAS-F16_BIAS-shift+16-F16_MANT_WIDTH-1)<<F32_MANT_WIDTH) 
	    + (abs<<(shift+F32_MANT_WIDTH-16+1)));
    }
    
    if (abs>=F16_EXP_MASK)			// infinity with sign or NaN
	return (abs==F16_EXP_MASK) ? (sign|F32_EXP_MASK) : F32_NAN;

    return (sign | (abs<<(F32_MANT_WIDTH-F16_MANT_WIDTH)))
	+ ((uint_fast32_t)(F32_BIAS-F16_BIAS)<<F32_MANT_WIDTH);	// normal
}


float64 float64_from_float16(float16 f)
{
    uint_fast64_t sign = ((uint_fast64_t)f&F16_SIGN_MASK)<<48;
    uint_fast16_t abs = f&(F16_EXP_MASK|F16_MANT_MASK);
    
    if (abs<=F16_MANT_MASK)
    {
	if (abs==0) return sign;		// zero

	uint_fast64_t shift = clz_16(abs);	// subnormal
	return sign | (((F64_BIAS-F16_BIAS-shift+16-F16_MANT_WIDTH-1)<<F64_MANT_WIDTH) 
	    + (abs<<(shift+F64_MANT_WIDTH-16+1)));
    }
    
    if (abs>=F16_EXP_MASK)			// infinity with sign or NaN
	return (abs==F16_EXP_MASK) ? (sign|F64_EXP_MASK) : F64_NAN;

    return (sign | (abs<<(F64_MANT_WIDTH-F16_MANT_WIDTH)))
	+ ((uint_fast64_t)(F64_BIAS-F16_BIAS)<<F64_MANT_WIDTH);	// normal

}


float16 float16_from_float32(float32 f)
{
    uint_fast32_t sign = (f>>16)&F16_SIGN_MASK;
    uint_fast32_t abs = f&0x7fffffff;

    if (abs >= ((F32_BIAS+31-F16_BIAS)<<F32_MANT_WIDTH))
	// already Nan or infinity or too big => infinity
	return (abs > F32_EXP_MASK) ? F16_NAN : (sign | F16_EXP_MASK);

    if (abs < ((uint_fast32_t)(F32_BIAS-F32_MANT_WIDTH-1)<<F32_MANT_WIDTH))
	return sign;			// too small => zero

    uint_fast32_t mant = (f&F32_MANT_MASK) 
	+ F32_MANT_MASK + ((uint32_t)1<<(F32_MANT_WIDTH-F16_MANT_WIDTH-1))
	+ ((f>>(F32_MANT_WIDTH-F16_MANT_WIDTH)) & 1);
    uint_fast32_t exp = abs>>F32_MANT_WIDTH;

    if (mant >= ((uint_fast32_t)2<<F32_MANT_WIDTH))
    {
	if (exp >= F32_BIAS-F16_BIAS)
	    return sign | (((exp-(F32_BIAS-F16_BIAS))<<10) 
		+ (mant>>(F32_MANT_WIDTH-F16_MANT_WIDTH+1))); 	// normal
    }
    else if (exp > F32_BIAS-F16_BIAS)
	return sign | (((exp-(F32_BIAS-F16_BIAS+1))<<10) 
	    + (mant>>(F32_MANT_WIDTH-F16_MANT_WIDTH))); 	// normal

    return sign | (mant >> (F32_BIAS-exp-F16_BIAS+F32_MANT_WIDTH-F16_MANT_WIDTH+1));	// subnormal
}


float16 float16_from_float64(float64 f)
{
    uint_fast64_t sign = (f>>48)&F16_SIGN_MASK;
    uint_fast64_t abs = f&0x7fffffffffffffff;

//printf("abs=%lx border=%lx\n", abs, ((uint_fast64_t)(F64_BIAS-F64_MANT_WIDTH-1)<<F64_MANT_WIDTH)); 
    if (abs >= ((uint64_t)(F64_BIAS+31-F16_BIAS)<<F64_MANT_WIDTH))
	// already Nan or infinity or too big => infinity
	return (abs > F64_EXP_MASK) ? F16_NAN : (sign | F16_EXP_MASK);

    if (abs < ((uint_fast64_t)(F64_BIAS-F64_MANT_WIDTH-1)<<F64_MANT_WIDTH))
	return sign;			// too small => zero

    uint_fast64_t mant = (f&F64_MANT_MASK) 
	+ F64_MANT_MASK + ((uint64_t)1<<(F64_MANT_WIDTH-F16_MANT_WIDTH-1))
	+ ((f>>(F64_MANT_WIDTH-F16_MANT_WIDTH)) & 1);
    uint_fast64_t exp = abs>>F64_MANT_WIDTH;


    if (mant >= ((uint_fast64_t)2<<F64_MANT_WIDTH))
    {
//printf("A");
	if (exp >= F64_BIAS-F16_BIAS)
	    return sign | (((exp-(F64_BIAS-F16_BIAS))<<10) 
		+ (mant>>(F64_MANT_WIDTH-F16_MANT_WIDTH+1))); 	// normal
    }
    else if (exp > F64_BIAS-F16_BIAS)
    {
//printf("B");
	return sign | (((exp-(F64_BIAS-F16_BIAS+1))<<10) 
	    + (mant>>(F64_MANT_WIDTH-F16_MANT_WIDTH))); 	// normal
    }
//printf("C");
//printf("abs=%lx mant=%lx exp=%lx shifted=%lx\n", abs, mant, exp, 
//(mant >> (F64_BIAS-exp-F16_BIAS+F64_MANT_WIDTH-F16_MANT_WIDTH+1-1)));
//    return sign | (mant >> (F64_BIAS-1-exp+28));	// subnormal
    return sign | (mant >> (F64_BIAS-exp-F16_BIAS+F64_MANT_WIDTH-F16_MANT_WIDTH+1));	// subnormal
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


#define MAX_ITER 100000000


double rand_double()
{
    uint64_t	u, v;
    double	*p = (double *)&v;

    u = rand();
    v = (u<<32) | rand();
    return *p;
}


int main()
{
    float16 a, b, correct, test;


    for (i=0; i<0x10000; i++)
    {
	a = i;

	b = float16_from_float64(float64_from_float16(a));
	if (a!=b && b!=F16_NAN)
        {
	    printf("%04x (%g float:%g) != %04x (%g)\n",
		    a, double_from_float16(a), float_from_float16(a),
		    b, double_from_float16(b));
	    return 1;
	}

/*
	for (j=0; j<0x10000; j++)
	{
	    b=j;
	    float sum = float_from_float16(a) + float_from_float16(b);
	    correct = float16_from_float(sum);
	    test = float16_add(a, b);
	    
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
	}
*/

	for (j=0; j<0x10000; j++)
	{
	    b=j;
	    double sum = double_from_float16(a) + double_from_float16(b);
	    correct = float16_from_double(sum);
	    test = float16_add(a, b);
	    
	    if (correct!=test) // works only, as NaN is always default value
	    {
		printf("%04x (%g %016lx) + %04x (%g %016lx) = "
		    "%04x (%g %016lx exact: %g) but: %04x (%g)\n",
		    a, double_from_float16(a), float64_from_float16(a),
		    b, double_from_float16(b), float64_from_float16(b),
		    correct, float_from_float16(correct), float64_from_float16(correct),
		    sum, //hex_from_float(sum),
		    test, float_from_float16(test));
		return 1;
	    }
	}

    }
    return 0;
}
