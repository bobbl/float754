// Modified version, that look more like hardware with multiple stages
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>



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
    uint64_t rem;
    uint64_t lo_mant, hi_mant;


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
	hi_mant = F64_EXTRACT_MANT(a);
	lo_mant = F64_EXTRACT_MANT(b);
	hi_sign = a & 0x8000000000000000;

	if (((a ^ b) & 0x8000000000000000) == 0)
	{
	    // both numbers have the same sign -> addition

	    if (hi_exp==0x7fe)
		// by the addition the exponent is increased by 1
		// and then it is too big 
		return hi_sign | 0x7ff0000000000000;

	    // nothing to do if both numbers are subnormal
	    if (hi_exp != 0)
	    {
		rem = (hi_mant ^ lo_mant) << 63;
		hi_mant = ((hi_mant>>1) | 0x0010000000000000) + (hi_mant & lo_mant & 1);
		lo_mant = lo_mant >> 1;
		hi_exp++;
	    }
	}
	else
	{
	    // not the same sign
	}
    }
    else 
    {
	if (diff_exp > 0)
	{
	    int16_t h=hi_exp; hi_exp=lo_exp; lo_exp=h;
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

	rem = (diff_exp==0) ? 0 : lo_mant << (64-diff_exp);
	    // correction, if diff_exp==0 (only possible if lo_exp==0)
	    // because (lo_mant << 64) == (lo_mant << 0) != 0 on intel platforms

	lo_mant >>= diff_exp;
	hi_mant |= 0x10000000000000;
    }




    if (((a ^ b) & 0x8000000000000000) == 0)
    {
	// same sign
	// simplifies the normalisation, but overflow checks are needed
	sum = hi_mant + lo_mant;

	if (hi_exp != 0)
	{
	    // both numbers are subnormal. Therefore simple adding is enough.
	    // If there is an overflow, i.e. the sum is normal, bit 52 of the
	    // result is set, which means that the exponent is 1 and not 0.
	    // Hence the conversion from subnormal to normal works
	    // automatically, only the sign has to be added.
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
    }
    else
    {
	    // not the same sign
	    sum = hi_mant - lo_mant;

	    // only possible, if exponents are equal
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

	    // Only for subbormal numbers the exponent need not be decreased.
	    // This can alos only happen, if exponents are equal.
	    if (hi_exp > 0)
		hi_exp--;

	    // normalise mantissa
	    while (sum < 0x10000000000000)
	    {
		if (hi_exp==0) 
		    // result is subnormal
		    return hi_sign | sum;
		sum = (sum<<1) - (rem>>63);
		rem <<= 1;
		hi_exp--;
	    }

	    if (rem == 0x8000000000000000)
    		sum = sum & (~1);    		// round to even
	    else if (rem > 0x8000000000000000)
    		sum--;				// round down
	    // else round up (nothing to do)
    }
    // Now the leading 1 must be masked out. But it is more efficient to
    // decrement the exponent by 1 and then add the implicit 1.
    // The decrementation as already been done earlier.
    return hi_sign | (((uint64_t)hi_exp<<52) + sum);
}

inline double float64_to_double(float64 f)
{
    return *(double *)&f;
}

inline float64 double_to_float64(double d)
{
    return *(uint64_t *)&d;
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
    int iter;
    uint64_t *pa, *pb, *pc, *pd;
    double a, b, c, d;
    
    pa = (uint64_t *)&a;
    pb = (uint64_t *)&b;
    pc = (uint64_t *)&c;
    pd = (uint64_t *)&d;
    for (iter=0; iter<MAX_ITER; iter++)
    {
	a = rand_double();
	b = rand_double();
	c = a+b;
	d = float64_to_double(float64_add(double_to_float64(a), double_to_float64(b)));

	if ((*pc != *pd) && !(isnan(c) && isnan(d)))
	    printf("%d: %016lx (%g) + %016lx (%g) = %016lx (%g) but: %016lx (%g)\n",
	    iter, *pa, a, *pb, b, *pc, c, *pd, d);
/*
	else
	    printf("%d: %g + %g = %g\n",
	    iter, a, b, c);
*/
    }
    return 0;
}



