# Modular Exponentiation with Precomputation with GMP

In ElGamal, we need to compute the modular exponentiation for the fixed base (e.g., the generator, the public key).

The functions in this repo allow precomputation -- a look-up table is generated.

It reduces the number of modular multiplication.

The idea of spending more space for such a table was mentioned in the following paper in EUROCRYPT'92:
*Fast Exponentiation with Precomputation*
Ernest F. Brickell, Daniel M. Gordon, Kevin S. McCurley, and DavidB. Wilson
https://www.ccrwest.org/gordon/fast.pdf

If you are looking for precomputation of modular exponentiation for bilinear mapping, I would recommend the following library by Alin Tomescu. @alinush
https://github.com/alinush/libbilinear
