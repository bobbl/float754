/* This is the ISC license (similar to the MIT/Expat license)
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



uint_fast16_t count_leading_zeros(uint_fast64_t x)
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

    return (64 - (x & 0x7f));
}


typedef uint64_t float64;

#define F64_EXC_INVALID_OP	1
#define F64_EXC_INEXACT		2

#define F64_DEFAULT_NAN		0xffffffffffffffff

#define F64_EXTRACT_EXP(f)	(((f)>>52)&0x7ff)
#define F64_EXTRACT_MANT(f)	((f)&0x000fffffffffffff)


void raise_exception(int e)
{
}


float64 float64_add(float64 a, float64 b)
{
    int16_t lo_exp = F64_EXTRACT_EXP(a);
    int16_t hi_exp = F64_EXTRACT_EXP(b);
    int16_t diff_exp = lo_exp - hi_exp;
    float64 hi_sign;
    int64_t sum;

    if (diff_exp == 0)
    {
	if (hi_exp == 0x7ff)
	{
	    // infinity or NaN
	    if ( ((a|b) & 0x000fffffffffffff) == 0)
	    {
		// two times infinity (mantissa==0)
		if (a != b) // already equal except the sign
		{
		    // infinity minus infinity raises an error
		    raise_exception(F64_EXC_INVALID_OP);
		    return F64_DEFAULT_NAN;
		}
		else
		    return a;
	    }
	    if ( ((a&b) & 0x0008000000000000) == 0)
		// one of the two NaNs is signalling 
		// (highest bit of the mantissa==0)
		raise_exception(F64_EXC_INVALID_OP);
	    return (a>b) ? a : b;
		// return greater mantissa
	}

	// same exponent
	uint64_t a_mant = F64_EXTRACT_MANT(a);
	uint64_t b_mant = F64_EXTRACT_MANT(b);
	hi_sign = a & 0x8000000000000000;

	if (((a ^ b) & 0x8000000000000000) == 0)
	{
	    // both numbers have the same sign -> addition
	    sum = a_mant + b_mant;

	    if (hi_exp==0x7fe)
		// by the addition the exponent is increased by 1
		// and then it is too big 
		return hi_sign | 0x7ff0000000000000;
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

	    if (sum < 0)
	    {
	        hi_sign ^= 0x8000000000000000;
		sum = -sum;
	    }
	    else if (sum == 0)
	    {
//	    	if (float64_rounding_mode == F64_ROUND_DOWN)
//		    return 0x8000000000000000;
//	   	else
		    return 0;
	    }

	    // very easy for subnormal numbers
	    if (hi_exp == 0)
		return hi_sign | sum;

	    // normalise mantissa
	    uint_fast8_t zeros = count_leading_zeros(sum) - 11;
	    sum <<= zeros;
    	    hi_exp = hi_exp - 1 - zeros;

	    // subnormal result
	    if (hi_exp < 0)
	        return hi_sign | (sum >> (-hi_exp));
/*
	    // only for subbormal numbers the exponent need not be decreased
	    if (hi_exp > 0)
		hi_exp--;

	    // normalise mantissa
	    uint_fast64_t border = 0x0000000000200000;
	    uint_fast16_t step;
	    for (step=32; step>0;)
	    {
		if (sum < border)
		{
		    hi_exp -= step;
		    sum <<= step;
		}
		step >>= 1;
		border <<= step;
	    }
	    if (hi_exp < 0)
		return hi_sign | (sum << (-hi_exp));
*/
/*
	    while (sum < 0x0010000000000000)
	    {
		if (hi_exp == 0) 
		    // result is subnormal
		    return hi_sign | sum;
		sum <<= 1;
		hi_exp--;
	    }
*/
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
	    lo_mant = F64_EXTRACT_MANT(b);
	    hi_mant = F64_EXTRACT_MANT(a);
	    hi_sign = a;

	}
	else // (diff_exp < 0)
	{
	    diff_exp = -diff_exp;
	    lo_mant = F64_EXTRACT_MANT(a);
	    hi_mant = F64_EXTRACT_MANT(b);
	    hi_sign = b;
	}
	// temporarily the lower bits of hi_sign are not masked out,
	// in order to return the correct value in the case of an exception

	if (hi_exp == 0x7ff)	        	// infinity or NaN
	{
	    if (((hi_mant & 0x8000000000000)==0) && (hi_mant!=0))
	        // signalling NaN
	        raise_exception(F64_EXC_INVALID_OP);
	    return hi_sign;
	}
	if (diff_exp > 53)			// lo much too low
	{
	    if (lo_exp!=0 || lo_mant!=0)
		raise_exception(F64_EXC_INEXACT);
    	    return hi_sign;
	}

	// from now on, only the sign of hi_sign is needed
	hi_sign &= 0x8000000000000000;

	if (lo_exp == 0)
	    // lo subnormal
    	    diff_exp--;
	else
	    // insert implicit 1
    	    lo_mant |= 0x10000000000000;

	uint64_t rem = (diff_exp==0) ? 0 : lo_mant << (64-diff_exp);
	    // correction, if diff_exp==0 (only possible if lo_exp==0)
	    // because (lo_mant << 64) == (lo_mant << 0) != 0 on intel platforms

	lo_mant >>= diff_exp;
	hi_mant |= 0x10000000000000;

	if (((a ^ b) & 0x8000000000000000) == 0)
	{
	    // same sign
	    // simplifies the normalisation, but overflow checks are needed
	    sum = hi_mant + lo_mant;

	    if (sum < 0x20000000000000)			// no bit overflow
	    {
		if (rem == 0x8000000000000000)
		    sum = (sum+1) & (~1);    		// round to even
		else if (rem > 0x8000000000000000)
		    sum++;				// round up
		// else round down (nothing to do)

		hi_exp--; 
	    }
	    else					// 1 bit overflow
	    {
		if (hi_exp>=0x7fe)
		    return hi_sign | 0x7ff0000000000000;	// infinity

		if ((sum&1) == 0)
		    sum = sum>>1;			// round down
		else if (rem == 0)
		    sum = ((sum>>1)+1) & (~1);    	// round to even
		else
		    sum = (sum>>1) + 1;			// round up

// alternate formulation:
//		if (((sum&1)!=0) && rem == 0)
//		    sum = ((sum>>1) + 1) & (~1);		// round up to even
//		else
//		    sum = (sum>>1) + (sum&1);			// round to nearest

// another alternate formulation:
//		sum = ((sum>>1) + (sum&1)) & ((rem==0 && ((sum&1)!=0)) ? (~1) : (~0));		// round to nearest or even
	    }
	}
	else
	{
/*
	    // not the same sign
	    sum = hi_mant - lo_mant;
	    hi_exp--; 

	    // normalise mantissa
	    uint_fast64_t border = 0x0000000000200000;
	    uint_fast16_t step;
	    for (step=32; step>0;)
	    {
		if (sum < border)
		{
		    hi_exp -= step;
		    sum = (sum << step) | (rem >> (64-step));
		    rem <<= step;
		}
		step >>= 1;
		border <<= step;
	    }
	    if (hi_exp < 0)
		return hi_sign | (sum << (-hi_exp));
*/
#define REM_MASK ((1L<<(52+2))-1)
#define REM_HALF (1L<<(52+1))
	    // not the same sign
	    sum = hi_mant - lo_mant;

	    if (rem != 0)
		sum--;

	    uint_fast8_t zeros = count_leading_zeros(sum) - (64-1-52);

	    sum = (sum << zeros) | (((-rem) & ((1L<<54)-1))
		>> ((52+2)-zeros));
	    rem = (rem << zeros) & ((1L<<54)-1);
	    if (rem == REM_HALF)
    		sum = (sum+1) & (~1);    		// round to even
	    else if (rem != 0 && rem < REM_HALF)
		sum++;

	    hi_exp = hi_exp - 1 - zeros;
	    if (hi_exp < 0)
		return hi_sign | (sum >> (-hi_exp));
/*
	    while (sum < 0x10000000000000)
	    {
		if (hi_exp==0) 
		    // result is subnormal
		    return hi_sign | sum;
		sum = (sum<<1) - (rem>>63);
		rem <<= 1;
		hi_exp--;
	    }
*/

	    if (rem == 0x8000000000000000)
    		sum = sum & (~1);    		// round to even
	    else if (rem > 0x8000000000000000)
    		sum--;				// round down
	    // else round up (nothing to do)
	}
    }
    // Now the leading 1 must be masked out. But it is more efficient to
    // decrement the exponent by 1 and then add the implicit 1.
    // The decrementation as already been done earlier.
    return hi_sign | (((uint64_t)hi_exp<<52) + sum);
}













#define clz_64(x) (__builtin_clzl(x))
#define clz_32(x) (__builtin_clzl(x)-32)
#define clz_16(x) (__builtin_clzl(x)-48)


#undef REM_HALF
typedef uint16_t float16;
typedef uint32_t float32;

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

#define WIDTH 64

FLOAT(WIDTH) _float64_add(FLOAT(WIDTH) a, FLOAT(WIDTH) b)
{
    EXP_TYPE(WIDTH) lo_exp = EXTRACT_EXP(WIDTH, a);
    EXP_TYPE(WIDTH) hi_exp = EXTRACT_EXP(WIDTH, b);
    EXP_TYPE(WIDTH) diff_exp = lo_exp - hi_exp;
    FLOAT(WIDTH) hi_sign;
    int_fast16_t sum;

    if (diff_exp == 0)
    {
	if (hi_exp == EXP_MAX(WIDTH))
	{
	    // infinity or NaN
	    if ( ((a|b) & MANT_MASK(WIDTH)) == 0)
	    {
		// two times infinity (mantissa==0)
		if (a != b) // already equal except the sign
		{
		    // infinity minus infinity raises an error
		    raise_exception(EXC_INVALID_OP);
		    return NOTANUM(WIDTH);
		}
		else
		    return a;
	    }
	    if ( ((a&b) & SIGNALLING_MASK(WIDTH)) == 0)
		// one of the two NaNs is signalling 
		// (highest bit of the mantissa==0)
		raise_exception(EXC_INVALID_OP);
	    return NOTANUM(WIDTH);
	}

	// same exponent
	uint64_t a_mant = EXTRACT_MANT(WIDTH, a);
	uint64_t b_mant = EXTRACT_MANT(WIDTH, b);
	hi_sign = a & SIGN_MASK(WIDTH);

	if (((a ^ b) & SIGN_MASK(WIDTH)) == 0)
	{
	    // both numbers have the same sign -> addition
	    sum = a_mant + b_mant;

	    if (hi_exp == EXP_MAX(WIDTH)-1)
		// by the addition the exponent is increased by 1
		// and then it is too big 
		return hi_sign | EXP_MASK(WIDTH);
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

	    if (sum == 0) return 0;
	    if (sum < 0) {
	        hi_sign ^= SIGN_MASK(WIDTH);
		sum = -sum;
	    }

	    // very easy for subnormal numbers
	    if (hi_exp == 0)
		return hi_sign | sum;

	    // normalise mantissa
	    uint_fast8_t zeros = CLZ(WIDTH)(sum) - (WIDTH-1-MANT_WIDTH(WIDTH));
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
	    lo_mant = EXTRACT_MANT(WIDTH, b);
	    hi_mant = EXTRACT_MANT(WIDTH, a);
	    hi_sign = a;

	}
	else // (diff_exp < 0)
	{
	    diff_exp = -diff_exp;
	    lo_mant = EXTRACT_MANT(WIDTH, a);
	    hi_mant = EXTRACT_MANT(WIDTH, b);
	    hi_sign = b;
	}
	// temporarily the lower bits of hi_sign are not masked out,
	// in order to return the correct value in the case of an exception

	if (hi_exp == EXP_MAX(WIDTH))
	{
	    if (hi_mant==0)
		return hi_sign;				// infinity
	    if ((hi_mant & SIGNALLING_MASK(WIDTH))==0)
	        raise_exception(EXC_INVALID_OP);        // signalling NaN
	    return NOTANUM(WIDTH);			// NaN
	}

	if (diff_exp > MANT_WIDTH(WIDTH)+2)		// lo much too low
// +1 für float64?
	{
	    if (lo_exp!=0 || lo_mant!=0)
		raise_exception(EXC_INEXACT);
    	    return hi_sign;
	}

	// from now on, only the sign of hi_sign is needed
	hi_sign &= SIGN_MASK(WIDTH);

	if (lo_exp == 0)
	    // lo subnormal
    	    diff_exp--;
	else
	    // insert implicit 1
    	    lo_mant |= 1L<<MANT_WIDTH(WIDTH);


//#define REM_MASK ((1<<(MANT_WIDTH(WIDTH)+2))-1)
//#define REM_HALF (1<<(MANT_WIDTH(WIDTH)+1))
//#define REM_MASK ((MANT_MASK(WIDTH)<<2)|3)
//#define REM_HALF (1<<(MANT_WIDTH(WIDTH)+1))
//#define REM_HALF ((MANT_WIDTH(WIDTH)<<1)+2)
#define REM_HALF (2L<<MANT_WIDTH(WIDTH))
	uint64_t rem = (lo_mant << ((MANT_WIDTH(WIDTH)+2)-diff_exp)) 
	    & ((MANT_MASK(WIDTH)<<2)|3);

	lo_mant >>= diff_exp;
	hi_mant |= 1L<<MANT_WIDTH(WIDTH);

	if (((a ^ b) & SIGN_MASK(WIDTH)) == 0)
	{
	    // same sign
	    // simplifies the normalisation, but overflow checks are needed
	    sum = hi_mant + lo_mant;

	    if (sum < (2L<<MANT_WIDTH(WIDTH)))		// no bit overflow
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
		if (hi_exp>=EXP_MAX(WIDTH)-1)
		    return hi_sign | EXP_MASK(WIDTH);	// infinity

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

	    uint_fast8_t zeros = CLZ(WIDTH)(sum) - (WIDTH-1-MANT_WIDTH(WIDTH));

	    sum = (sum << zeros) | (((-rem) & ((MANT_MASK(WIDTH)<<2)|3))
		>> ((MANT_WIDTH(WIDTH)+2)-zeros));
	    rem = (rem << zeros) & ((MANT_MASK(WIDTH)<<2)|3);
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
    return hi_sign | (((UINT_FAST(WIDTH))hi_exp<<MANT_WIDTH(WIDTH)) + sum);
}




























static inline double float64_to_double(float64 f)
{
    return *(double *)&f;
}

static inline float64 double_to_float64(double d)
{
    return *(uint64_t *)&d;
}






//#define MAX_ITER 100000000
#define MAX_ITER 10


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
    int iter;
    uint64_t *pa, *pb, *pc, *pd, *pe;
    double a, b, c, d, e;
    
    pa = (uint64_t *)&a;
    pb = (uint64_t *)&b;
    pc = (uint64_t *)&c;
    pd = (uint64_t *)&d;
    pe = (uint64_t *)&e;
    for (iter=0; iter<MAX_ITER; iter++)
    {
	a = rand_double();
	b = rand_double();
	c = a+b;
	d = float64_to_double(float64_add(double_to_float64(a), double_to_float64(b)));
	e = float64_to_double(_float64_add(double_to_float64(a), double_to_float64(b)));

	if ((*pc != *pd) && !(isnan(c) && isnan(d)))
	    printf("%d: %016lx (%g) + %016lx (%g) = %016lx (%g) but: %016lx (%g)\n",
	    iter, *pa, a, *pb, b, *pc, c, *pd, d);

	if ((*pc != *pe) && !(isnan(c) && isnan(e)))
	    printf("X%d: %016lx (%g) + %016lx (%g) = %016lx (%g) but: %016lx (%g)\n",
	    iter, *pa, a, *pb, b, *pc, c, *pe, e);
	else
	    printf("%d:  %016lx (%g) +  %016lx (%g) =  %016lx (%g)\n",
	    iter, *pa, a, *pb, b, *pc, c);
    }
    return 0;
}
