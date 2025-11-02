# AffChab1
Sage code for computing ℤ[ζ₃]-points on a hyperelliptic genus 2 curve via the Affine Chabauty method

The file `ExImag.sage` contains the computations for Theorem 7.2 of [LL25]. Using the Affine Chabauty method developed in that paper it is verified that the affine hyperelliptic curve defined by

$$ y^2 = x^6 - 4x^5 + 2x^4 + 6x^3 + x^2 - 10x + 1 $$

has no integral points other than $(0, \pm 1), (2, \pm 1), (1, \pm\sqrt{-3})$. The curve has genus 2 and Mordell–Weil rank 2 over ℚ(ζ₃), so the classical Chabauty method does not apply. The code has been tested with Sage 10.7.

## References
- [LL25] Marius Leonhardt and Martin Lüdtke, "Affine Chabauty I"

## Authors

- Marius Leonhardt
- Martin Lüdtke
