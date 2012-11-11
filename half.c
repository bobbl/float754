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

#define EXC_INVALID_OP	1
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
//#define _BIAS(width)		F ## width ## _BIAS
#define _UINT_FAST(width)	uint_fast ## width ## _t
#define _FLOAT(width)		float ## width
#define _CLZ(width)		clz_ ## width
#define MANT_WIDTH(width)	_MANT_WIDTH(width)
//#define BIAS(width)		_BIAS(width)
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

float16 float16_add(float16 a, float16 b)
{
    int_fast8_t lo_exp = EXTRACT_EXP(16, a);
    int_fast8_t hi_exp = EXTRACT_EXP(16, b);
    int_fast16_t diff_exp = lo_exp - hi_exp;
    float16 hi_sign;
    int_fast16_t sum;

    if (diff_exp == 0)
    {
	if (hi_exp == EXP_MAX(16))
	{
	    // infinity or NaN
	    if ( ((a|b) & MANT_MASK(16)) == 0)
	    {
		// two times infinity (mantissa==0)
		if (a != b) // already equal except the sign
		{
		    // infinity minus infinity raises an error
		    raise_exception(EXC_INVALID_OP);
		    return NOTANUM(16);
		}
		else
		    return a;
	    }
	    if ( ((a&b) & SIGNALLING_MASK(16)) == 0)
		// one of the two NaNs is signalling 
		// (highest bit of the mantissa==0)
		raise_exception(EXC_INVALID_OP);
	    return NOTANUM(16);
	}

	// same exponent
	uint64_t a_mant = EXTRACT_MANT(16, a);
	uint64_t b_mant = EXTRACT_MANT(16, b);
	hi_sign = a & SIGN_MASK(16);

	if (((a ^ b) & SIGN_MASK(16)) == 0)
	{
	    // both numbers have the same sign -> addition
	    sum = a_mant + b_mant;

	    if (hi_exp==(EXP_MAX(16)-1))
		// by the addition the exponent is increased by 1
		// and then it is too big 
		return hi_sign | EXP_MASK(16);
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
	        hi_sign ^= SIGN_MASK(16);
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
	    lo_mant = EXTRACT_MANT(16, b);
	    hi_mant = EXTRACT_MANT(16, a);
	    hi_sign = a;

	}
	else // (diff_exp < 0)
	{
	    diff_exp = -diff_exp;
	    lo_mant = EXTRACT_MANT(16, a);
	    hi_mant = EXTRACT_MANT(16, b);
	    hi_sign = b;
	}
	// temporarily the lower bits of hi_sign are not masked out,
	// in order to return the correct value in the case of an exception

	if (hi_exp == EXP_MAX(16))
	{
	    if (hi_mant==0)
		return hi_sign;				// infinity
	    if ((hi_mant & 0x0200)==0)
	        raise_exception(EXC_INVALID_OP);        // signalling NaN
	    return NOTANUM(16);				// NaN
	}

	if (diff_exp > MANT_WIDTH(16)+2)		// lo much too low
	{
	    if (lo_exp!=0 || lo_mant!=0)
		raise_exception(EXC_INEXACT);
    	    return hi_sign;
	}

	// from now on, only the sign of hi_sign is needed
	hi_sign &= SIGN_MASK(16);

	if (lo_exp == 0)
	    // lo subnormal
    	    diff_exp--;
	else
	    // insert implicit 1
    	    lo_mant |= 1<<MANT_WIDTH(16);


#define REM_BITS (MANT_WIDTH(16)+2)
#define REM_MASK ((1<<REM_BITS)-1)
#define REM_HALF (1<<(REM_BITS-1))
	uint64_t rem = (lo_mant << (REM_BITS-diff_exp)) & REM_MASK;

	lo_mant >>= diff_exp;
	hi_mant |= 1<<MANT_WIDTH(16);

	if (((a ^ b) & SIGN_MASK(16)) == 0)
	{
	    // same sign
	    // simplifies the normalisation, but overflow checks are needed
	    sum = hi_mant + lo_mant;

	    if (sum < (2<<MANT_WIDTH(16)))		// no bit overflow
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
		if (hi_exp>=EXP_MAX(16)-1)
		    return hi_sign | EXP_MASK(16);	// infinity

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
    return hi_sign | (((uint64_t)hi_exp<<MANT_WIDTH(16)) + sum);
}

/*
#define D_BIAS (uint_fast32_t)(BIAS(32)-BIAS(16))
#define D_MANT (MANT_WIDTH(32)-MANT_WIDTH(16))

float32 float32_from_float16(float16 f)
{
    uint_fast32_t sign = ((uint_fast32_t)f&SIGN_MASK(16))<<16;
    uint_fast32_t abs = f&(EXP_MASK(16)|MANT_MASK(16));
    
    if (abs<=MANT_MASK(16))
    {
	if (abs==0) return sign;		// zero

	uint_fast32_t shift = clz_16(abs);	// subnormal
	return sign | (((D_BIAS-MANT_WIDTH(16)+16-1-shift)<<MANT_WIDTH(32)) 
	    + (abs<<(MANT_WIDTH(32)-16+1+shift)));
    }
    
    if (abs>=EXP_MASK(16))			// infinity with sign or NaN
	return (abs==EXP_MASK(16)) ? (sign|EXP_MASK(32)) : NOTANUM(32);

    return (sign | (abs<<D_MANT)) + (D_BIAS<<MANT_WIDTH(32));
}

float16 float16_from_float32(float32 f)
{
    uint_fast32_t sign = (f>>16)&SIGN_MASK(16);
    uint_fast32_t abs = f&(EXP_MASK(32)|MANT_MASK(32));

    if (abs >= ((D_BIAS+EXP_MAX(16))<<MANT_WIDTH(32)))
	// already Nan or infinity or too big => infinity
	return (abs > EXP_MASK(32)) ? NOTANUM(16) : (sign | EXP_MASK(16));

    if (abs < ((uint_fast32_t)(BIAS(32)-MANT_WIDTH(32)-1)<<MANT_WIDTH(32)))
	// too small => zero
	return sign;

    uint_fast32_t mant = (f&MANT_MASK(32)) + MANT_MASK(32)
	+ ((uint_fast32_t)1<<(D_MANT-1)) + ((f>>D_MANT) & 1);
    uint_fast32_t exp = abs>>MANT_WIDTH(32);

    if (mant >= ((uint_fast32_t)2<<MANT_WIDTH(32)))
    {
	if (exp >= D_BIAS)
	    return sign | (((exp-D_BIAS)<<MANT_WIDTH(16)) + (mant>>(D_MANT+1)));
    }
    else if (exp > D_BIAS)
	return sign | (((exp-D_BIAS-1)<<MANT_WIDTH(16)) + (mant>>D_MANT));

    // subnormal
    return sign | (mant >> (D_BIAS+D_MANT+1-exp));
}
#undef D_BIAS
#undef D_MANT
*/

#define SMALL 16
#define BIG 32

float32 float32_from_float16(float16 f)
{
    UINT_FAST(BIG) sign = ((UINT_FAST(BIG))f&SIGN_MASK(SMALL))<<(BIG-SMALL);
    UINT_FAST(SMALL) abs = f&(SIGN_MASK(SMALL)-1);
    
    if (abs<=MANT_MASK(SMALL))
    {
	if (abs==0) return sign;			// zero

	UINT_FAST(BIG) shift = SMALL-1-CLZ(SMALL)(abs);	// subnormal
	return sign | (((BIAS(BIG)-BIAS(SMALL)-MANT_WIDTH(SMALL)+shift)<<MANT_WIDTH(BIG)) 
	    + (abs<<(MANT_WIDTH(BIG)-shift)));
    }
    
    if (abs<EXP_MASK(SMALL))				// representable
	return (sign | (abs<<(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))))
	    + ((BIAS(BIG)-BIAS(SMALL))<<MANT_WIDTH(BIG));

    if (abs==EXP_MASK(SMALL))				// infinity with sign
	return (sign|EXP_MASK(BIG));
	
    return NOTANUM(BIG); 				// NaN
}

float16 float16_from_float32(float32 f)
{
    UINT_FAST(SMALL) sign = (f>>(BIG-SMALL))&SIGN_MASK(SMALL);
    EXP_TYPE(BIG) exp = ((f>>MANT_WIDTH(BIG))&EXP_MAX(BIG)) - BIAS(BIG) + BIAS(SMALL);

    if (exp >= EXP_MAX(SMALL))
	// return NaN if already Nan or infinity if infinity or too big
	return ((f&(SIGN_MASK(BIG)-1)) > EXP_MASK(BIG)) 
	    ? NOTANUM(SMALL) : (sign | EXP_MASK(SMALL));

    if (exp < BIAS(SMALL)-MANT_WIDTH(BIG)-1)
	// too small => zero
	return sign;

    UINT_FAST(BIG) mant = ((f&MANT_MASK(BIG)) + MANT_MASK(BIG)
	+ (MANT_MASK(BIG)>>(MANT_WIDTH(SMALL)+1)) + 1
	+ ((f>>(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))) & 1))
	>> (MANT_WIDTH(BIG)-MANT_WIDTH(SMALL));

    if (mant >= (MANT_MASK(SMALL)+1)<<1)
    {
	if (exp >= 0)
	    return sign | ((exp<<MANT_WIDTH(SMALL)) + (mant>>1));
    }
    else if (exp > 0)
	return sign | (((exp-1)<<MANT_WIDTH(SMALL)) + mant);

    // subnormal
    return sign | (mant >> (1-exp));
}


#undef SMALL
#undef BIG
#define SMALL 16
#define BIG 64

float64 float64_from_float16(float16 f)
{
    UINT_FAST(BIG) sign = ((UINT_FAST(BIG))f&SIGN_MASK(SMALL))<<(BIG-SMALL);
    UINT_FAST(SMALL) abs = f&(SIGN_MASK(SMALL)-1);
    
    if (abs<=MANT_MASK(SMALL))
    {
	if (abs==0) return sign;			// zero

	UINT_FAST(BIG) shift = SMALL-1-CLZ(SMALL)(abs);	// subnormal
	return sign | (((BIAS(BIG)-BIAS(SMALL)-MANT_WIDTH(SMALL)+shift)<<MANT_WIDTH(BIG)) 
	    + (abs<<(MANT_WIDTH(BIG)-shift)));
    }
    
    if (abs<EXP_MASK(SMALL))				// representable
	return (sign | (abs<<(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))))
	    + ((BIAS(BIG)-BIAS(SMALL))<<MANT_WIDTH(BIG));

    if (abs==EXP_MASK(SMALL))				// infinity with sign
	return (sign|EXP_MASK(BIG));
	
    return NOTANUM(BIG); 				// NaN
}

float16 float16_from_float64(float64 f)
{
    UINT_FAST(SMALL) sign = (f>>(BIG-SMALL))&SIGN_MASK(SMALL);
    EXP_TYPE(BIG) exp = ((f>>MANT_WIDTH(BIG))&EXP_MAX(BIG)) - BIAS(BIG) + BIAS(SMALL);

    if (exp >= EXP_MAX(SMALL))
	// return NaN if already Nan or infinity if infinity or too big
	return ((f&(SIGN_MASK(BIG)-1)) > EXP_MASK(BIG)) 
	    ? NOTANUM(SMALL) : (sign | EXP_MASK(SMALL));

    if (exp < BIAS(SMALL)-MANT_WIDTH(BIG)-1)
	// too small => zero
	return sign;

    UINT_FAST(BIG) mant = ((f&MANT_MASK(BIG)) + MANT_MASK(BIG)
	+ (MANT_MASK(BIG)>>(MANT_WIDTH(SMALL)+1)) + 1
	+ ((f>>(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))) & 1))
	>> (MANT_WIDTH(BIG)-MANT_WIDTH(SMALL));

    if (mant >= (MANT_MASK(SMALL)+1)<<1)
    {
	if (exp >= 0)
	    return sign | ((exp<<MANT_WIDTH(SMALL)) + (mant>>1));
    }
    else if (exp > 0)
	return sign | (((exp-1)<<MANT_WIDTH(SMALL)) + mant);

    // subnormal
    return sign | (mant >> (1-exp));
}


#undef SMALL
#undef BIG
#define SMALL 32
#define BIG 64

float64 float64_from_float32(float32 f)
{
    UINT_FAST(BIG) sign = ((UINT_FAST(BIG))f&SIGN_MASK(SMALL))<<(BIG-SMALL);
    UINT_FAST(SMALL) abs = f&(SIGN_MASK(SMALL)-1);
    
    if (abs<=MANT_MASK(SMALL))
    {
	if (abs==0) return sign;			// zero

	UINT_FAST(BIG) shift = SMALL-1-CLZ(SMALL)(abs);	// subnormal
	return sign | (((BIAS(BIG)-BIAS(SMALL)-MANT_WIDTH(SMALL)+shift)<<MANT_WIDTH(BIG)) 
	    + (abs<<(MANT_WIDTH(BIG)-shift)));
    }
    
    if (abs<EXP_MASK(SMALL))				// representable
	return (sign | (abs<<(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))))
	    + ((BIAS(BIG)-BIAS(SMALL))<<MANT_WIDTH(BIG));

    if (abs==EXP_MASK(SMALL))				// infinity with sign
	return (sign|EXP_MASK(BIG));
	
    return NOTANUM(BIG); 				// NaN
}

float32 float32_from_float64(float64 f)
{
    UINT_FAST(SMALL) sign = (f>>(BIG-SMALL))&SIGN_MASK(SMALL);
    EXP_TYPE(BIG) exp = ((f>>MANT_WIDTH(BIG))&EXP_MAX(BIG)) - BIAS(BIG) + BIAS(SMALL);

    if (exp >= EXP_MAX(SMALL))
	// return NaN if already Nan or infinity if infinity or too big
	return ((f&(SIGN_MASK(BIG)-1)) > EXP_MASK(BIG)) 
	    ? NOTANUM(SMALL) : (sign | EXP_MASK(SMALL));

    if (exp < BIAS(SMALL)-MANT_WIDTH(BIG)-1)
	// too small => zero
	return sign;

    UINT_FAST(BIG) mant = ((f&MANT_MASK(BIG)) + MANT_MASK(BIG)
	+ (MANT_MASK(BIG)>>(MANT_WIDTH(SMALL)+1)) + 1
	+ ((f>>(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))) & 1))
	>> (MANT_WIDTH(BIG)-MANT_WIDTH(SMALL));

    if (mant >= (MANT_MASK(SMALL)+1)<<1)
    {
	if (exp >= 0)
	    return sign | ((exp<<MANT_WIDTH(SMALL)) + (mant>>1));
    }
    else if (exp > 0)
	return sign | (((exp-1)<<MANT_WIDTH(SMALL)) + mant);

    // subnormal
    return sign | (mant >> (1-exp));
}

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
	if (i <= -(1<<(EXP_MAX(16)-BIAS(16))))
	    return SIGN_MASK(16) | EXP_MASK(16); // - infinity
	sign = SIGN_MASK(16);
	i = -i;
    }
    else
    {
	if (i==0)
	    return 0;
	if (i >= (1<<(EXP_MAX(16)-BIAS(16))))
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

float double_from_float64(float64 f)
{
    uf64_t r;
    r.u = f;
    return r.f;
}


#define MAX_ITER 100000000

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
	

/*
	for (j=0; j<0x10000; j++)
	{
	    b=j;
	    test = float16_add(a, b);

	    float sum = float_from_float16(a) + float_from_float16(b);
	    correct = float16_from_float(sum);
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
*/
    }
    return 0;
}
