/* Conversion between floating point types
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
   Set the macros BIG and SMALL to the bit width of the types. Two
   functions floatSMALL_from_floatBIG() and floatBIG_from_floatSMALL() are
   generated. */
   
#define _FUNC_HEADER(T, F) FLOAT(T) float ## T ## _from_float ## F (FLOAT(F) ff)
#define FUNC_HEADER(T, F) _FUNC_HEADER(T, F)

// big from small
FUNC_HEADER(BIG, SMALL)
{
    UF(SMALL) s; s.f = ff;
    UF(BIG) b; b.u = ((UINT_FAST(BIG))s.u&SIGN_MASK(SMALL))<<(BIG-SMALL);
    UINT_FAST(SMALL) abs = s.u&(SIGN_MASK(SMALL)-1);
    
    if (abs<=MANT_MASK(SMALL))
    {
	if (abs!=0) { // else  zero
	    UINT_FAST(BIG) shift = SMALL-1-CLZ(SMALL)(abs); // subnormal
	    b.u |= (((BIAS(BIG)-BIAS(SMALL)-MANT_WIDTH(SMALL)+shift)<<MANT_WIDTH(BIG)) 
		+ (abs<<(MANT_WIDTH(BIG)-shift)));
	}
    } else if (abs<EXP_MASK(SMALL)) {			// representable
	b.u = (b.u | (abs<<(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))))
	    + ((BIAS(BIG)-BIAS(SMALL))<<MANT_WIDTH(BIG));
    } else if (abs==EXP_MASK(SMALL)) {			// infinity with sign
	b.u |= EXP_MASK(BIG);
    } else
	b.u = NOTANUM(BIG); 				// NaN
    return b.f;
}

FUNC_HEADER(SMALL, BIG)
{
    UF(BIG) b; b.f = ff;
    UF(SMALL) s; s.u = (b.u>>(BIG-SMALL))&SIGN_MASK(SMALL);
    EXP_TYPE(BIG) exp = ((b.u>>MANT_WIDTH(BIG))&EXP_MAX(BIG)) - BIAS(BIG) + BIAS(SMALL);

    if (exp >= EXP_MAX(SMALL)) {
	// return NaN if already Nan or infinity if infinity or too big
	s.u = ((b.u&(SIGN_MASK(BIG)-1)) > EXP_MASK(BIG)) 
	    ? NOTANUM(SMALL)
	    : (s.u | EXP_MASK(SMALL));

    } else if (exp >= BIAS(SMALL)-MANT_WIDTH(BIG)-2) {
	UINT_FAST(BIG) mant = ((b.u&MANT_MASK(BIG)) + MANT_MASK(BIG)
	    + (MANT_MASK(BIG)>>(MANT_WIDTH(SMALL)+1)) + 1
	    + ((b.u>>(MANT_WIDTH(BIG)-MANT_WIDTH(SMALL))) & 1))
	    >> (MANT_WIDTH(BIG)-MANT_WIDTH(SMALL));

	if (mant >= (MANT_MASK(SMALL)+1)<<1) {
	    if (exp >= 0) {
		s.u |= ((exp<<MANT_WIDTH(SMALL)) + (mant>>1));
		return s.f;
	    }
	} else if (exp > 0) {
	    s.u |= (((exp-1)<<MANT_WIDTH(SMALL)) + mant);
	    return s.f;
	}

	// subnormal
	mant = (b.u&MANT_MASK(BIG)) + MANT_MASK(BIG) + 1; // add implicit 1
	UINT_FAST(BIG) corr = (mant << (1+MANT_WIDTH(SMALL)-1+exp)) & (2*MANT_MASK(BIG)+1);
	mant >>= MANT_WIDTH(BIG) - MANT_WIDTH(SMALL) + 1 - exp;
	corr = (corr > (MANT_MASK(BIG)+1)) 
	    | ((corr == (MANT_MASK(BIG)+1)) & (mant & 1));
	s.u |= (mant + corr);
    } // else too small -> zero
    return s.f;
}

#undef _FUNC_HEADER
#undef FUNC_HEADER
