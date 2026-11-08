# Verification record

This file lists every check of the theorem that we can document publicly, with dates, tools and authors. No
professional mathematical peer review has been carried out. The full Erdős Problem #203 remains open.

## What has to be checked

1. **The pool is complete and correct.** Every listed factor is prime, the product of the listed prime powers equals
   G_N = gcd(2^N − 1, 3^N − 1) for N = 12,252,240, and therefore the listed primes are exactly the primes p > 3
   with e_p | N.
2. **The arithmetic of the bound.** For every pool prime, ord_p(2), ord_p(3), ord_p(6), e_p and d_p are correct.
   The weights are nonnegative, sum to 1 and avoid the classes 0 (mod 3) and 0 (mod 4). Table 1 and
   B = 3,961,943 / 4,084,080 < 1 follow exactly.
3. **The argument.** Lemma 1 (line density), Lemma 2 (weighted obstruction), the CRT normalisation of 5 and 7 and
   the reduction to the finite pool are correct. The extension to the 480 divisors of N is immediate, because a
   period dividing a divisor of N also divides N.
4. **Optional:** the dual certificate of Appendix A is feasible and matches B, so the weights are optimal within the
   fixed cut. The theorem does not depend on this.

## Up to 1 November 2026: checks within the author's workflow

The article refers to an initial frozen report, to a GMP computation and to an informal review (Sections 3 and 5).
These were made within the author's private AI-assisted workflow. They are not documented here, and the article's
word "independent" for them is too strong.

## 1 November 2026: review by Erdős Problem a Day

At the author's request, Erdős Problem a Day reviewed the submission with its own code and reported the limited claim
as verified. Erdős Problem a Day is run by Patrick White with Claude, an AI. The review was sent in private
correspondence, and it is theirs to publish; it is not reproduced here. Its suggestions are reflected in article v1.0
or in the README of this repository:

- the corrected digit count of G_N;
- the monotonicity argument for the divisors of N;
- definitions of "pure" and "covering period";
- statements of the torus reduction and of the line lemma;
- attribution of the sources.

## 2 November 2026: confirmation of the corrected article

These checks were made by ChatGPT/Codex within the author's workflow. ChatGPT/Codex wrote both checkers and
reconstructed `optional_lp_dual.json` from the tables of the article. Status: `CONFIRMED_IN_STATED_SCOPE`.

- **Exact standard-library checker** (the original of `verify/verify_dual_stdlib.py`): PASS. It covers checks 1
  (primality and orders), 2 and 4. It also enumerates the projection directly for 10 representatives and all 12 CRT
  normalisation pairs. Its record reports that a wrong factor product, a deleted dual entry and a wrong projection
  were rejected. Of these, only the deleted dual entry is fed back through a check in the code; the other two follow
  from the checks that pass.
- **GMP checker under WSL (GMP 6.3.0):** PASS, with the direct gcd equal to the product of the factors. The SHA-256 of
  the 732-digit decimal string of G_N, without a trailing newline, is
  `1befdff4d20ed0ab7e5371250ba9408e735c966de9fee71131dd2caae0b7daa5`. The file `data/G_N_decimal.txt` ends with a
  newline, so its own SHA-256 is the one in `SHA256SUMS`. `tr -d '\n' < data/G_N_decimal.txt | sha256sum` gives
  `1befdff4…b7daa5`.
- **Written review by ChatGPT** of the argument of Sections 1–4 of the article, covering Lemmas 1–2, the CRT
  normalisation, the pool reduction, the divisor extension and Proposition 1: no error found.

## 8 November 2026: public verifiers in this repository

These checks were made by Claude (Anthropic) as research operator, on the files exactly as published here.

| Check | Platform | Result | Output |
|---|---|---|---|
| `verify/verify_dual_stdlib.py` | Windows 11, Python 3.12 | PASS, B = 3961943/4084080, gap 122137/4084080, 192 primes, 480 divisors | `EXACT_CHECK_windows.json` |
| `verify/verify_dual_stdlib.py` | WSL2 Ubuntu, Python 3.12.3 | PASS (same values) | `EXACT_CHECK_wsl.json` |
| `verify/verify_pool_gmp.py` | Windows 11, GMP 6.3.0 (`libgmp-10.dll` from Git for Windows) | PASS: gcd = product = `G_N_decimal.txt`, 2,432 bits, 732 digits; perturbed data rejected | `GMP_CHECK_windows.json` |
| `verify/verify_pool_gmp.py` | WSL2 Ubuntu, GMP 6.3.0 (`libgmp.so.10`) | PASS (same values) | `GMP_CHECK_wsl.json` |
| `verify/independent_check.py` | Windows 11, GMP 6.3.0 | PASS on all checks | `INDEPENDENT_CHECK_windows.json` |
| `verify/independent_check.py` | WSL2 Ubuntu, GMP 6.3.0 | PASS on all checks | `INDEPENDENT_CHECK_wsl.json` |

The decimal G_N has the same SHA-256 as on 2 November (`1befdff4…b7daa5`, the digits without newline).

`verify/verify_dual_stdlib.py` differs from the 2 November checker only in its header, file paths and output
handling, as the header lists; the mathematical checks are unchanged. `verify/verify_pool_gmp.py` is a public rewrite
of the 2 November GMP checker. It adds a gmpy2 backend, finds the library on Windows and macOS, compares the result
with `data/G_N_decimal.txt`, and adds perturbation controls.

## 8 November 2026: three from-scratch checks by separate AI agents

Three separate Claude agent instances each received the claim and the article and wrote their own code. None had
access to the earlier checkers or to the other agents' work, and each was told to look for errors rather than to
confirm the claim. These are independent implementations, run within the author's workflow; they are not an
independent human verification.

**1. Pool completeness and arithmetic. Verdict: CONFIRMED.**
- G_N was computed with GMP 6.3.0 on Windows and under WSL, each in two ways: directly, and as 3^N mod (2^N − 1)
  followed by a gcd. It was also computed a third time with plain CPython integers (Mersenne folding plus `math.gcd`,
  129 s), which rules out a GMP-specific error. All runs agree: 2,432 bits, 732 decimal digits, and SHA-256 of the
  digits `1befdff4…b7daa5`. `mpz_sizeinbase(G, 10)` returns 733, an overestimate by one.
- All 192 listed factors are prime by trial division. The exponents equal min(v_p(2^N − 1), v_p(3^N − 1)), and the
  product equals G_N exactly.
- As a separate check without GMP, a segmented sieve tested all 203,280,219 primes 3 < p < 2^32 with exact 64-bit
  modular powers. Exactly 192 satisfy 2^N ≡ 3^N ≡ 1 (mod p), and they are the listed ones.

**2. Adversarial review of the argument. Verdict: CONFIRMED_WITH_CAVEATS; the caveats are expository, not errors.**
- Lemma 1, the averaging modulo 60, the reduction to the torus, the CRT normalisation, Lemma 2, the extension to the
  480 divisors and the weak-duality argument of Proposition 1 were all checked and found correct.
- Lemma 1 was brute-forced for every coset of all 38 pool primes with e_p ≤ 260.
- The pool was found by brute force up to 3·10^7, and G_N divided by the pool primes leaves cofactor 1.
- B, Table 1 and the dual were recomputed independently. All 30 dual inequalities are tight.
- The theorem and abstract are correctly scoped, and the prior work [3] is identified and credited.
- With uniform weights on the 30 admissible classes, the method of [3] reproduces that post's published values
  (3965/4752 at period 55,440; 11255/12012 at 720,720). For the present pool, uniform weights give
  1557209/1531530 > 1, so the non-uniform weights are needed.
- Caveats: three steps are left implicit in the text (the primes 2 and 3, an unused 5 or 7, and averaging over
  infinitely many lines); see *Notes on the article* in the README. The article's word "independent" is too strong.
  Novelty should not be overstated: the new parts are the complete N-pool and explicit non-uniform weights for it.

**3. Exact bound. Verdict: CONFIRMED.**
- The pool was recomputed: GMP gcd, then trial division until the cofactor is 1. As a cross-check, every prime up to
  2·10^7 was tested with the criterion 2^gcd(p−1,N) ≡ 3^gcd(p−1,N) ≡ 1 (mod p). Both give the same 192 primes and
  valuations.
- ord_p(2), ord_p(3) and ord_p(6) were found by factoring p − 1. |⟨2, 3⟩| = e_p was confirmed by enumeration for the
  180 primes with e_p ≤ 2·10^5.
- B was computed exactly, grouped by d and also prime by prime: B = 3961943/4084080, N·B = 11,885,829. All 10 rows of
  Table 1 match. The primes per group d = 1, 2, 3, 4, 5, 6, 10, 12, 15, 20 number 94, 29, 24, 10, 10, 12, 3, 4, 3, 1.
- For every coset of all 78 pool primes with e_p ≤ 1500, a direct computation on the torus shows that the largest
  weighted contribution equals (d_p/e_p)·μ_{d_p} exactly. A separate run also confirmed this for the primes with
  d = 10 or 20 (2861, 4421, 19891 and 51481, with e_p up to 25,740).
- For every divisor T of N, the bound for T is at most B. All dual budget constraints are tight, all 30 support
  inequalities hold with L = B, and complementary slackness holds. All 2,000 random admissible weight vectors give a
  bound of at least L.

**Consolidated script.** `verify/independent_check.py` consolidates code from these checks into one script that starts
from N alone. It adds two controls: dropping the largest prime from the pool must leave cofactor 8,168,161, and uniform
weights must fail to give B < 1.

## 8 November 2026: which periods are new

`verify/new_periods.py` takes the earlier exclusions from the Erdős Problem a Day post of 14 October 2026: every
period up to 65,520 and every divisor of 110880, 166320, 221760, 332640 or 720720. It counts 411 of the 480
divisors of N as previously excluded and lists the other 69, from 67,320 to N, all multiples of 17. Its output is
`NEW_PERIODS.csv`. For each new period it computes the pool, the density Σ 1/e_p and the line bound with the
article's weights and with uniform weights; every article-weight bound is below 1. By the weakest sufficient argument:

- 50 periods have pool density below 1;
- 18 fall to uniform weights, with the largest bound 0.9092 at 942,480;
- only N itself needs the non-uniform weights.

## 8 November 2026: pre-publication review of this package

Five further Claude agent instances reviewed the package before publication, each from one angle: claims and scope,
privacy, reproducibility from a fresh clone, mathematical accuracy, and metadata. Every medium or high finding was
then checked by another instance that tried to refute it. The confirmed findings, for example wording on
independence, credit for prior work, the number of verifiers and the CFF reference types, were fixed before release.

## What has not been done

- No professional mathematical peer review.
- No exhaustive search of the literature for historical novelty.
- No check of the full Erdős Problem #203, which remains open.
