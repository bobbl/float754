/* Floating point addition and subtraction
 * Licenced under the ISC license (similar to the MIT/Expat license)
 *
 * Copyright (c) 2012 Jörg Mische <bobbl@gmx.de>
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

/* Include file that just works when included into a proper C source file.
   Set the macro width to the bit width of the types. Then the functions
   floatWIDTH_add() and floatWIDTH_sub() are generated. */

   
#define _FUNC_ADD(W) float ## W ## _add
#define FUNC_ADD(W) _FUNC_ADD(W)
#define _FUNC_SUB(W) float ## W ## _sub
#define FUNC_SUB(W) _FUNC_SUB(W)
#define _FUNC_MUL(W) float ## W ## _mul
#define FUNC_MUL(W) _FUNC_MUL(W)

#define BIGONE ((UINT_FAST(WIDTH))1)

FLOAT(WIDTH) FUNC_ADD(WIDTH) (FLOAT(WIDTH) af, FLOAT(WIDTH) bf)
{
    UF(WIDTH) a, b, r;
    a.f = af;
    b.f = bf;
    EXP_TYPE(WIDTH) lo_exp = EXTRACT_EXP(WIDTH, a.u);
    EXP_TYPE(WIDTH) hi_exp = EXTRACT_EXP(WIDTH, b.u);
    EXP_TYPE(WIDTH) diff_exp = lo_exp - hi_exp;

    if (diff_exp == 0) {
	// same exponent
	if (hi_exp == EXP_MAX(WIDTH)) {
	    // infinity or NaN
	    r.u = (a.u==b.u && (a.u & MANT_MASK(WIDTH))==0)
		? a.u // inf+inf = inf | -inf-inf = -inf
		: NOTANUM(WIDTH);

	} else if (((a.u ^ b.u) & SIGN_MASK(WIDTH)) == 0) {
	    // both numbers have the same sign -> addition
	    if (hi_exp == EXP_MAX(WIDTH)-1) {
		// by the addition the exponent is increased by 1
		// and then it is too big 
		r.u = (a.u & SIGN_MASK(WIDTH)) | EXP_MASK(WIDTH);
	    } else if (hi_exp == 0) {
		// both numbers are subnormal. If there is an overflow, i.e. the
		// sum is normal, bit 52 of the result is set, which means that
		// the exponent is 1 and not 0. Hence the conversion from subnormal
		// to normal works automatically, only the sign has to be added.
		r.u = a.u + (b.u & MANT_MASK(WIDTH));
	    } else {
		// because of the two implicit ones the result has one digit
		// to much, i.e. the mantissa has to be shifted right by 1
		// and the exponent must be increased by 1
		INT_FAST(WIDTH) sum 
		    = (a.u&MANT_MASK(WIDTH)) + (b.u&MANT_MASK(WIDTH));
		r.u = (a.u|MANT_MASK(WIDTH)) + 1 + ROUND(sum, 0);
	    }
	} else {
	    // numbers have different signs -> subtraction
	    INT_FAST(WIDTH) sum = EXTRACT_MANT(WIDTH, a.u) - EXTRACT_MANT(WIDTH, b.u);

	    if (sum == 0) {
		r.u = 0;
	    } else {
		if (sum < 0) {
	    	    r.u = b.u & SIGN_MASK(WIDTH);
		    sum = -sum;
		} else {
		    r.u = a.u & SIGN_MASK(WIDTH);
		}

		// very easy for subnormal numbers
		if (hi_exp == 0) {
		    r.u |= sum;
		} else {
		    // normalise mantissa
		    uint_fast8_t zeros = CLZ(WIDTH)(sum) - (WIDTH-1-MANT_WIDTH(WIDTH));
		    sum <<= zeros;
    		    hi_exp -= zeros;
		    r.u |= (hi_exp<1)
			? (sum >> (1-hi_exp)) // subnormal result
			: (((UINT_FAST(WIDTH))(hi_exp-1)<<MANT_WIDTH(WIDTH)) + sum);
			    // Now the leading 1 must be masked out. But it is more efficient to
			    // decrement the exponent by 1 and then add the implicit 1.
		}
	    }
	}

    } else {
	UINT_FAST(WIDTH) lo_mant, hi_mant;
	
	if (diff_exp > 0) {
	    hi_exp += diff_exp;
	    lo_exp -= diff_exp;
		// exchange hi_exp and lo_exp
	    lo_mant = EXTRACT_MANT(WIDTH, b.u);
	    hi_mant = EXTRACT_MANT(WIDTH, a.u);
	    r.u = a.u;
	} else {
	    diff_exp = -diff_exp;
	    lo_mant = EXTRACT_MANT(WIDTH, a.u);
	    hi_mant = EXTRACT_MANT(WIDTH, b.u);
	    r.u = b.u;
	}
	// temporarily the lower bits of r.u are not masked out,
	// in order to return the correct value in the case of an exception

	if (hi_exp == EXP_MAX(WIDTH)) {
	    if (hi_mant!=0) {
		r.u = NOTANUM(WIDTH); // x+nan = nan
	    } // x+inf = inf

	} else if (diff_exp <= MANT_WIDTH(WIDTH)+2) { // lo not too low
	    // from now on, only the sign of r.u is needed
	    r.u &= SIGN_MASK(WIDTH);

	    if (lo_exp == 0)
    		diff_exp--; // lo subnormal
	    else
    		lo_mant |= BIGONE<<MANT_WIDTH(WIDTH); // insert implicit 1

	    if (((a.u ^ b.u) & SIGN_MASK(WIDTH)) == 0) {
		// same sign
		// simplifies the normalisation, but overflow checks are needed
		INT_FAST(WIDTH) sum = (lo_mant >> diff_exp) + hi_mant;
		if (sum < (BIGONE<<MANT_WIDTH(WIDTH))) {		// no bit overflow
		    lo_mant <<= 1;
		    sum = (lo_mant >> diff_exp) + 2*hi_mant;
		} else if (hi_exp<EXP_MAX(WIDTH)-1) {		// 1 bit overflow
		    sum += (BIGONE<<MANT_WIDTH(WIDTH));
		} else {					
		    r.u |= EXP_MASK(WIDTH); // infinity
		    return r.f;
		}
		UINT_FAST(WIDTH) rem = lo_mant & ((BIGONE<<(diff_exp))-1);
		r.u |= ((UINT_FAST(WIDTH))hi_exp<<MANT_WIDTH(WIDTH)) + ROUND(sum, rem);
		    // Now the leading 1 must be masked out. But it is more efficient
		    // to decrement the exponent by 1 and then add the implicit 1.

	    } else {
		// not the same sign
		UINT_FAST(WIDTH) rem;
		INT_FAST(WIDTH) sum = (hi_mant | (BIGONE<<MANT_WIDTH(WIDTH)))
		    - (lo_mant>>diff_exp)
		    - (((lo_mant & ((BIGONE<<(diff_exp))-1))!=0) ? 1 : 0);

		if (diff_exp > 1) {
		    if (sum < (BIGONE<<MANT_WIDTH(WIDTH))) {
			sum = (sum << 2) | (((-lo_mant) >> (diff_exp-2)) & 3);
			rem = lo_mant & ((BIGONE<<(diff_exp-2))-1);
			hi_exp = hi_exp - 1;
		    } else {
		        sum = (sum << 1) | (((-lo_mant) >> (diff_exp-1)) & 1);
			rem = lo_mant & ((BIGONE<<(diff_exp-1))-1);
		    }
		} else {
		    uint_fast8_t zeros = (sum==0)
			? MANT_WIDTH(WIDTH)+2
			: CLZ(WIDTH)(sum) + MANT_WIDTH(WIDTH) - WIDTH + 2;
		    // in LLVM clz(0)=clz(1)+1, but in gcc it is undefined,
		    // hence we must deal with this special case

		    sum = (sum << zeros)
			| (((-lo_mant) << (zeros-diff_exp)) & ((BIGONE<<zeros)-1));
		    rem = 0;
		    hi_exp = hi_exp - zeros + 1;
		}
		r.u |= (hi_exp < 1)
		    ? (ROUND(sum, rem) >> (1-hi_exp))
		    : (((UINT_FAST(WIDTH))(hi_exp-1)<<MANT_WIDTH(WIDTH)) + ROUND(sum, rem));
			// Now the leading 1 must be masked out. But it is more efficient to
			// decrement the exponent by 1 and then add the implicit 1.
	    }
	} // else lo much to low, return hi
    }
    return r.f;

//if (a.u==0x3000 && b.u==0xafff)
//    printf("sum=%lx rem=%lx hi_mant=%lx lo_mant=%lx hi_exp=%d lo_exp=%d\n",
//	sum, rem, hi_mant, lo_mant, hi_exp, lo_exp);
}


FLOAT(WIDTH) FUNC_SUB(WIDTH) (FLOAT(WIDTH) af, FLOAT(WIDTH) bf)
{
    UF(WIDTH) b;
    b.f = bf;
    b.u ^= SIGN_MASK(WIDTH);
    return FUNC_SUB(WIDTH) (af, b.f);
}


FLOAT(WIDTH) FUNC_MUL(WIDTH) (FLOAT(WIDTH) af, FLOAT(WIDTH) bf)
{
    UF(WIDTH) a, b, r;
    a.f = af;
    b.f = bf;
    r.u = (a.u^b.u) & SIGN_MASK(WIDTH); 
    UINT_FAST(WIDTH) a_mant = a.u & (SIGN_MASK(WIDTH)-1);
    UINT_FAST(WIDTH) b_mant = b.u & (SIGN_MASK(WIDTH)-1);

    if (a_mant>=EXP_MASK(WIDTH)) {
	r.u = (a_mant==EXP_MASK(WIDTH) && b_mant>0 && b_mant<=EXP_MASK(WIDTH))
	    ? (r.u | EXP_MASK(WIDTH)) // inf*inf = inf*real = inf
	    : NOTANUM(WIDTH); // inf*nan = inf*0 = nan*x = nan

    } else if (b_mant>=EXP_MASK(WIDTH)) {
	r.u = (b_mant>EXP_MASK(WIDTH) || a_mant==0)
	    ? NOTANUM(WIDTH) // 0*inf = 0*nan = real*nan = nan
	    : (r.u | EXP_MASK(WIDTH)); // real*inf = inf

    } else if (a_mant!=0 && b_mant!=0) { // else 0*0 = 0*real = real*0 = 0
	EXP_TYPE(WIDTH) a_exp = EXTRACT_EXP(WIDTH, a.u);
	EXP_TYPE(WIDTH) b_exp = EXTRACT_EXP(WIDTH, b.u);
	EXP_TYPE(WIDTH) r_exp = a_exp + b_exp - BIAS(WIDTH);
	a_mant &= MANT_MASK(WIDTH);
	b_mant &= MANT_MASK(WIDTH);

	if (a_exp==0) {
	    // shift subnormal into position and adjust exponent
	    EXP_TYPE(WIDTH) shift = CLZ(WIDTH)(a_mant)-WIDTH+MANT_WIDTH(WIDTH)+1;
	    r_exp -= shift-1;
	    a_mant <<= shift;
		// If b is also subnormal, the result is 0 anyway.
		// This is recognised later.
	} else if (b_exp==0) {
	    // shift subnormal into position and adjust exponent
	    EXP_TYPE(WIDTH) shift = CLZ(WIDTH)(b_mant)-WIDTH+MANT_WIDTH(WIDTH)+1;
	    r_exp -= shift-1;
	    b_mant <<= shift;
	} 
	a_mant |= MANT_MASK(WIDTH)+1;
	b_mant |= MANT_MASK(WIDTH)+1;

	// here a lot of platform dependent optimisation is possible
#if PLATFORM_WIDTH > WIDTH
	UINT_FAST(PLATFORM_WIDTH) product
	    = (UINT_FAST(PLATFORM_WIDTH))a_mant
	    * (UINT_FAST(PLATFORM_WIDTH))b_mant;
	if (product >= ((UINT_FAST(PLATFORM_WIDTH))2<<(2*MANT_WIDTH(WIDTH)))) {
	    product = (product >> 1) | (product & 1); // keep overflow bits
	    r_exp++;
	}
	UINT_FAST(WIDTH) r_mant = product >> (MANT_WIDTH(WIDTH)-1);
	UINT_FAST(WIDTH) remains = product & ((BIGONE << (MANT_WIDTH(WIDTH)-1))-1);
#elif PLATFORM_WIDTH == WIDTH
	UINT_FAST(WIDTH) ph, pl;
	UMUL_PPMM(WIDTH)(ph, pl, a_mant, b_mant);
	
	UINT_FAST(WIDTH) r_mant, remains;
	if (ph >= (2L<<(2*MANT_WIDTH(WIDTH)-WIDTH))) {
	    r_mant = (ph << (WIDTH-MANT_WIDTH(WIDTH))) | (pl >> (MANT_WIDTH(WIDTH)));
	    remains = pl & ((BIGONE << (MANT_WIDTH(WIDTH)))-1);
	    r_exp++;
	} else {
	    r_mant = (ph << (WIDTH-MANT_WIDTH(WIDTH)+1)) | (pl >> (MANT_WIDTH(WIDTH)-1));
	    remains = pl & ((BIGONE << (MANT_WIDTH(WIDTH)-1))-1);
	}
	
#else
#error "word width of platform is to small"
#endif

	if (r_exp >= EXP_MAX(WIDTH)) {
	    // overflow => +/- infinity
	    r.u |= EXP_MASK(WIDTH);
	
	} else if (r_exp <= 0) {
	    if (r_exp >= -MANT_WIDTH(WIDTH)-1) {
		// subnormal
		remains |= r_mant & ((BIGONE << (1-r_exp))-1);
		r_mant >>= (1-r_exp);
		r.u |= ROUND(r_mant, remains);
	    } // else +/-0
	} else {
	    // round to nearest or even
	    r_mant = ROUND(r_mant, remains);
	    if (r_mant>=(2*MANT_MASK(WIDTH)+2)) { // only == is possible
		r_mant <<= 1;
		r_exp++;
	    }

	    r.u |= ((UINT_FAST(WIDTH))r_exp<<MANT_WIDTH(WIDTH))
		| (r_mant&MANT_MASK(WIDTH));
	}
    }
    return r.f;
}





#undef _FUNC_ADD
#undef FUNC_ADD
#undef _FUNC_SUB
#undef FUNC_SUB
#undef _FUNC_MUL
#undef FUNC_MUL

#undef BIGONE
