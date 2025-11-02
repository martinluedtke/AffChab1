r"""
Sage code for Theorem 7.2 of the paper "Affine Chabauty I" [LL25].
It contains the computation of ℤ[ζ₃]-points on the affine hyperelliptic curve
  y² = x⁶ - 4x⁵ + 2x⁴ + 6x³ + x² - 10x + 1
with LMFDB label 1549.a.1549.1. The curve has genus 2 and Mordell–Weil rank 2 over ℚ(ζ₃).

REFERENCES:

 - [LL25] Marius Leonhardt and Martin Lüdtke, "Affine Chabauty I"

AUTHORS:
 - Marius Leonhardt
 - Martin Lüdtke
 
"""

p = 7
K = Qp(p,15)
x = polygen(K)
f = x^6 - 4*x^5 + 2*x^4 + 6*x^3 + x^2 - 10*x + 1
X = HyperellipticCurve(f)

α = sqrt(K(-3))

# known solutions (up to hyperelliptic involution)
P0 = X(0,1)
P1 = X(2,1)
P2 = X(1,α)

# hyperelliptic involution
def ι(P):
    return X(P[0],-P[1])

known_sols = [P0,P1,P2, ι(P0),ι(P1),ι(P2)]
known_sol_names = ["P₀", "P₁", "P₂", "ιP₀", "ιP₁", "ιP₂"]

print("Known solutions to y² = x⁶ - 4x⁵ + 2x⁴ + 6x³ + x² - 10x + 1 over ℤ[ζ₃]")
print("    (0,±1),  (2,±1),  (1,±√-3)")
print("")

print("Let P₀ = (0,1), P₁ = (2,1), P₂ = (1,√-3).")
print("Basis of log differentials H^0(X, Ω^1(D)):")
print("    ω₀ = dx/y,  ω₁ = x dx/y,  ω₂ = x^2 dx/y")
print(f"Auxiliary prime: p = (2-√-3) over rational prime 7")
print("")
print("Check that [P₁] - [P₀] and [P₂] - [P₀] are linearly independent in the Jacobian")
print("by computing the determinant of the 2x2 matrix (∫_{P₀}^{Pᵢ} ωⱼ):")
det = matrix([2*X.coleman_integrals_on_basis(P0,Q)[:2] for Q in [P1,P2]]).determinant()
print(f"    {det}")
print("From Magma's RankBound we know the rank of the Jacobian over K is at most 2,")
print("so [P₁] - [P₀] and [P₂] - [P₀] span a subgroup of full rank.")
print("")

print(f"The only bad prime for the minimal model is 1549 and there is only one component.")
print(f"So all O_K-points have the same reduction type.")
print("")

print("The Affine Chabauty Condition is satisfied, so there exists a nontrivial log differential η")
print("such that ∫_P^Q η vanishes for all pairs of O_K-points P,Q.")
print("Hence the function ρ(P) := det(∫_{P₀}^{Pᵢ} ωⱼ)_{1 ≤ i,j ≤ 3} with P₃ = P vanishes for all O_K-points P.")
print("")

chabauty_function = lambda P: matrix([2*X.coleman_integrals_on_basis(P0,Q)[:3] for Q in [P1,P2,P]]).determinant()

print("Check that the function vanishes on all known points:")
for i in range(len(known_sols)):
    print(f"    ρ({known_sol_names[i]}) = {chabauty_function(known_sols[i])}")
print("")

print(f"The function vanishes on integral points but not necessarily on rational points:")
print(f"    ρ((33/20, 4073/8000)) = {chabauty_function(X(33/20, 4073/8000))}")
print("")

print("Find coefficients αⱼ for the annihilating differential η = α₀ω₀ + α₁ω₁ + α₂ω₂ as 2x2 minors in the determinant equation:")

coeffs = [
    matrix([[2*X.coleman_integrals_on_basis(P0,Q)[j] for j in [1,2]] for Q in [P1,P2]]).determinant(),
    -matrix([[2*X.coleman_integrals_on_basis(P0,Q)[j] for j in [0,2]] for Q in [P1,P2]]).determinant(),
    matrix([[2*X.coleman_integrals_on_basis(P0,Q)[j] for j in [0,1]] for Q in [P1,P2]]).determinant()
]

print(f"    α₀ = {coeffs[0]}")
print(f"    α₁ = {coeffs[1]}")
print(f"    α₂ = {coeffs[2]}")
print("")


print(f"For each known point P, compute the power series expansion of ∫_P^{{Q(t)}} η in the parameter t with x(t) = x(P) + {p}t.")
print("Valuations of coefficients:")
R.<t> = PowerSeriesRing(K,'t')
for i in range(len(known_sols)):
    P = known_sols[i]
    xt = P[0] + p*t
    yt = sqrt(f(xt))
    if yt(0) != P[1]:
        yt = -yt
    g = sum(coeffs[j] * p*xt^j/yt for j in range(3)).integral()
    vals = [g[k].valuation() for k in range(15)]
    print(f"    P = {known_sol_names[i]}:\t {vals}")

print("In each case, -∞ is the only nonpositive slope, corresponding to t=0, so the known points are the only O_K-points in their residue disc.")
print("")

print("Check that known points cover all residue discs.")
x = polygen(GF(p))
f = x^6 - 4*x^5 + 2*x^4 + 6*x^3 + x^2 - 10*x + 1
modp_sols = []
for a in GF(p):
    if f(a) == 0:
        modp_sols.append((a,0))
    elif f(a).is_square():
        b = sqrt(f(a))
        modp_sols.append((a,b))
        modp_sols.append((a,-b))
print(f"𝔽₇-points: {modp_sols}")

known_sol_modp = [(0,1), (2,1), (1,2), (0,6), (2,6), (1,5)]
for i in range(len(known_sols)):
    print(f"    {known_sol_names[i]} reduces to {known_sol_modp[i]} modulo (2-√-3)")

print("All residue discs are covered, so we found all O_K-points.")
