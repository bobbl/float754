IEEE 754-2008 Floating point emulation and hardware
===================================================

This project has two objectives

 1. Fast C code for emulating IEEE 754-2008 floating point operations
 2. Based on 1. an efficient floating point hardware in VHDL

Roadmap
-------
 1. C implementation of half precision: add, multiply, divide
 2. Complete brute force checking of half precision arithmetics (2^32 cases)
 3. C implementation of single precision
 4. C implementation of double precision
 5. C implementation of quad precision
 6. Rearrangement of C code to better fit hardware demands
 7. Translation to VHDL

License
-------

Licensed under the ISC licence (similar to the MIT/Expat license).
