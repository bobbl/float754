/* Conversion between floating point types
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

/* Include file that just works when included into a proper C source file.
   Set the macros BIG and SMALL to the bit width of the types. Two
   functions floatSMALL_from_floatBIG() and floatBIG_from_floatSMALL() are
   generated. */
   
#define _FUNC_HEADER(T, F) FLOAT(T) float ## T ## _from_float ## F (FLOAT(F) f)
#define FUNC_HEADER(T, F) _FUNC_HEADER(T, F)


FUNC_HEADER(BIG, SMALL)
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

FUNC_HEADER(SMALL, BIG)
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

#undef _FUNC_HEADER
#undef FUNC_HEADER
