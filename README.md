# Erdős Problem #203: no pure finite prime covering with period dividing lcm(1,…,17) = 12,252,240

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24050485.svg)](https://doi.org/10.5281/zenodo.24050485)

This repository contains a partial result on [Erdős Problem #203](https://www.erdosproblems.com/203), its exact
certificate, three verifiers that anyone can run in seconds, and the list of the periods that are new.

**Theorem.** No pure finite prime covering for Erdős Problem #203 has period dividing
N = lcm(1,…,17) = 12,252,240. In particular, none has period equal to any of the 480 divisors of N. Most of these
periods (411 of the 480) had already been excluded by earlier work. The 69 new ones are listed in
`results/NEW_PERIODS.csv`; see *Prior work* and *What is new*.

**What this is not.** It is **not** a solution of Erdős Problem #203, which **remains open**. The theorem concerns
only *pure* coverings, built from prime-divisibility cosets with at most one coset per prime and nothing else. It says
nothing about coverings that use algebraic factorisations, about non-periodic mechanisms, or about whether a
suitable integer m exists. The article is a research draft; it has not been submitted to a journal or peer reviewed.

## The statement in brief

Erdős and Graham asked whether some integer m ≥ 1 with gcd(m, 6) = 1 makes 2^k·3^l·m + 1 composite for every
k, l ≥ 0. For a prime p > 3 let H_p = ⟨2, 3⟩ ⊆ F_p^× and e_p = |H_p| = lcm(ord_p 2, ord_p 3). The exponent pairs
(k, l) at which p divides 2^k·3^l·m + 1 form either the empty set or one coset of
L_p = {(k, l) : 2^k·3^l ≡ 1 (mod p)}, a subgroup of Z² of index e_p. A *pure finite prime covering* is a finite family
of such cosets, at most one per prime, whose union is Z². Its period is the lcm of the e_p used.

The proof has three ingredients:

1. **A finite prime pool.** A prime p > 3 can occur in a covering of period dividing N only if e_p | N, i.e. only if
   p divides G_N = gcd(2^N − 1, 3^N − 1). G_N has 2,432 bits and 732 decimal digits. Its complete factorisation has
   192 distinct primes, the largest being 8,168,161; 5, 7, 11, 13 and 17 occur squared.
2. **Line projection.** On each line k − l = τ, a coset of L_p covers a proportion 1/ord_p(6) of the points if τ lies
   in one fixed class a_p modulo i_p = e_p/ord_p(6), and none of them otherwise. Averaged over the lines
   τ ≡ t (mod 60), the proportion is d_p/e_p if t ≡ a_p (mod d_p) and 0 otherwise, where d_p = gcd(i_p, 60). Hence,
   with weights λ_t on the classes modulo 60 (step 3), p contributes at most (d_p/e_p)·μ_{d_p}, where
   μ_d = max_c Σ_{t ≡ c (mod d)} λ_t. This method comes from earlier work by Erdős Problem a Day (see *Prior work*).
3. **Explicit weights modulo 60.** Weights 1/20 on t ∈ {1, 26, 41, 46} and 1/10 on t ∈ {2, 5, 7, 10, 14, 19, 23, 58}
   avoid the classes used by 5 (mod 4) and 7 (mod 3) after a CRT translation, and give the exact bound

   B = 3,961,943 / 4,084,080 < 1   (gap 122,137 / 4,084,080).

   A covering would need B ≥ 1, a contradiction.

A rational dual certificate (article, Appendix A) also shows that these weights are optimal **within this fixed cut**:
direction (1, 1), modulus 60 and the 5/7 normalisation. This is not a claim of optimality among all obstructions, and
the theorem does not depend on it.

## Verify it yourself

All scripts need Python 3.9 or newer. The second, and the third unless it is run with `--no-gmp`, also need GMP
(see below). The commands below assume the repository root, but the scripts work from any directory. On Windows, use
`py -3` or `python` instead of `python3`.

```bash
python3 verify/verify_dual_stdlib.py
```

The first script uses only the standard library and finishes in about 2 s. It proves that all 192 listed factors are
prime (deterministic trial division) and computes ord_p(2), ord_p(3), ord_p(6) and e_p with minimality witnesses.
It recomputes Table 1 of the article and the bound B exactly with fractions, and checks the weights, all
25 dual entries, the 10 budget and 30 support inequalities and the 480 divisors. It also enumerates the line
projection directly for representative primes, and checks that a deliberately weakened dual certificate is rejected.
Expected: the printed JSON begins with `"status": "PASS_EXACT_ARITHMETIC_AND_FINITE_PROJECTION_CHECKS"`.

```bash
python3 verify/verify_pool_gmp.py
```

The second script computes gcd(2^N − 1, 3^N − 1) directly (2^N − 1 has 12,252,240 bits). It checks that the gcd
equals the product of the listed prime powers and the decimal number in `data/G_N_decimal.txt`, and that perturbed
data are rejected. Together with the primality check of the first script, this proves that the list is the
*complete* pool. It takes about 2 s and needs GMP: either `pip install gmpy2`, or a system GMP library
(`libgmp.so.10` on Linux, `libgmp.10.dylib` on macOS, `libgmp-10.dll` on Windows, for example the copy that ships with
Git for Windows). Pass `--gmp PATH` to name the library. Expected: `"status": "PASS"` with `"gcd_bits": 2432` and
`"gcd_exact_decimal_digits": 732`.

Any computer algebra system can repeat this check on its own. Compute `gcd(2^12252240 - 1, 3^12252240 - 1)` and
compare it with `data/G_N_decimal.txt`.

```bash
python3 verify/independent_check.py
```

The third script imports neither of the other two and does not read the factor list. Starting from N alone, it finds
the pool by brute force over all primes up to 3·10^7 and shows with GMP that G_N has no other prime factor
(cofactor 1). It recomputes the orders, Table 1 and B exactly. It verifies Lemma 1 by brute force for every coset of
the 38 pool primes with e_p ≤ 260, checks the 480 divisors and checks the dual certificate as transcribed from the
article. It takes about 15 s. With `--no-gmp` it runs without GMP, but then it covers only primes up to 3·10^7 and
reports a partial status.

```bash
python3 verify/new_periods.py
```

The fourth script does not verify the theorem; it determines what is new. It lists the divisors of N that the
earlier work did not exclude, and for each one gives its pool, the pool density and the bounds with the article's
weights and with uniform weights. It takes about 1 s, and its output is `results/NEW_PERIODS.csv`.

File integrity: `sha256sum -c SHA256SUMS` (Linux, macOS with coreutils, WSL or Git Bash).

Our outputs are in `results/`. `results/VERIFICATION.md` lists every check we can document publicly, and by whom.

## Contents

| Path | What it is |
|---|---|
| `paper/Erdos203_Article_Draft_2026-11-01.pdf` | The article, research draft v1.0 of 1 November 2026, unchanged (SHA-256 `fc77c1b9…3c661f`) |
| `data/pool_factors.json` | The factorisation of G_N: 192 pairs [p, v] |
| `data/pool_primes.csv` | For each pool prime: v, ord_p(2), ord_p(3), ord_p(6), e_p, i_p = e_p/ord_p(6), d_p = gcd(i_p, 60). A derived table for convenience; the verifiers recompute these values and do not read it. |
| `data/G_N_decimal.txt` | G_N = gcd(2^N − 1, 3^N − 1) in decimal (732 digits and a final newline) |
| `data/optional_lp_dual.json` | Weights, Table 1 and the dual certificate of Appendix A, as exact rationals |
| `verify/verify_dual_stdlib.py` | Exact checker of primality, orders, Table 1, B, the dual and the divisors (standard library only) |
| `verify/verify_pool_gmp.py` | Direct GMP computation of G_N and comparison with the factor list |
| `verify/independent_check.py` | Recomputation from N alone: pool, completeness, B, Lemma 1 by brute force, divisors, dual |
| `verify/new_periods.py` | Lists the divisors of N not excluded by earlier work, with pool, density and bounds for each |
| `results/NEW_PERIODS.csv` | The 69 new periods (output of `verify/new_periods.py`) |
| `results/` | Verifier outputs and the verification record |
| `SHA256SUMS` | Checksums of all files above |

## Prior work

The setup and the method come from earlier public work:

- **AnimishSharma**, on the [erdosproblems.com discussion thread for #203](https://www.erdosproblems.com/forum/thread/203)
  (21 June 2026): the formula |H_p| = lcm(ord_p 2, ord_p 3) and the density condition Σ 1/e_p ≥ 1 for finite
  coverings.
- **veljjanoski**, on the same thread (15 September 2026; code at https://github.com/veljjanoski/erdos203): the
  formulation as a covering by cosets of lattices of index e_p, and complete prime pools computed through gcds.
- **Erdős Problem a Day** (Patrick White and Claude),
  ["Erdős #203: along the lines of slope 1 the primes 5 and 7 leave half the lines to the rest of the pool"](https://erdosproblemaday.com/day/203-the-quotient-lemma-lines-close-the-open-pools),
  14 October 2026: the weighted condition on lines, the slope-one 5/7 cut with uniform weights modulo 60, and the
  exclusion of every period up to 65,520 and of every divisor of 720,720 = lcm(1,…,16). That post names the pool of
  lcm(1,…,17) as the next frontier. It credits the underlying coset-restriction lemma to classical work
  (B. H. Neumann 1954; Tomkinson; Z.-W. Sun).

## What is new

The new ingredients are the complete prime pool of N = lcm(1,…,17), and with it the pools of all divisors of N, plus
explicit non-uniform weights for N itself. Uniform weights give 1557209/1531530 > 1 for the N-pool, so they do not
suffice there.

Of the 480 divisors of N, the post above already excluded 411: every period up to 65,520 and every divisor of
720,720. The remaining 69 run from 67,320 up to N, and all of them are multiples of 17. As far as we know, none of
them had been excluded before. They are listed in `results/NEW_PERIODS.csv`, which `verify/new_periods.py`
reproduces. Once the pools are known, most of them need only the earlier tools:

| New periods | Weakest argument that excludes them |
|---|---|
| 50 | the pool density Σ 1/e_p is below 1 |
| 18 | the 5/7 cut with uniform weights gives a bound below 1 (largest 0.9092, at 942,480) |
| 1 (N itself) | only the non-uniform weights of the article give a bound below 1 |

A preliminary literature search is not a proof of historical novelty.

## Notes on the article

The article is included exactly as written on 1 November 2026. On 8 November 2026 a separate Claude agent instance
was asked to look for errors. This was not a peer review; see `results/VERIFICATION.md`. It found no mathematical
error, but noted three steps that the text leaves implicit:

- **The primes 2 and 3.** The article does not say why only primes p > 3 are considered. For gcd(m, 6) = 1, the
  prime 2 divides 2^k·3^l·m + 1 only when k = 0, and 3 only when l = 0. The cosets of primes p > 3 must therefore
  cover all pairs with k, l ≥ 1. All the cosets are periodic, so they then cover Z².
- **An unused 5 or 7.** The normalisation assumes that both 5 and 7 occur in the covering. If one of them does not
  occur, add an arbitrary coset of it first. This only enlarges the union, so a covering stays a covering.
- **Averaging over lines.** "Averaging over all lines τ ≡ t (mod 60)" means averaging over infinitely many lines. It
  is made rigorous by the remark at the end of the proof of Lemma 1: all densities can be computed on the finite
  torus (Z/NZ)², because every e_p and 60 divide N.

The word "independent" in the article (the abstract and Section 3) is too strong. Those checks, and the re-checks of
2 November, were all made within the author's own AI-assisted workflow. The from-scratch checks of 8 November were
written by separate Claude agent instances without access to the earlier checkers or to each other's code. They are
therefore independent implementations, but they are not an independent human verification.

Further remarks for readers of this repository:

- **Section 5** describes the supplement as containing the factor list "with original Git provenance" and "a GMP/C
  checker". The original experiment files and their Git history belong to the author's private research repository
  and are not included here. Correctness does not depend on them: the factor list is given as data, and its
  completeness is checked directly by `verify/verify_pool_gmp.py`, a Python script that calls GMP.
- **Appendix A** names `optional_lp_dual.json` and `verify_dual_stdlib.py`; they are in `data/` and `verify/`.
- **Section 3** already gives the corrected digit count of G_N: 732. An earlier internal report had 733, the output
  of GMP's `mpz_sizeinbase`, which may overestimate by one.
- **References.** The article credits the Erdős Problem a Day post but not the forum contributions listed under
  *Prior work*.

## Use of AI

The research was conceived and directed by the author. As the article states, the author used an AI-assisted research
workflow, including the author's private research system SAPIO MORPHEUS / COMPLEX and ChatGPT (OpenAI), for algorithm
design, computation, checking and preparation of the draft. The PDF was produced in ChatGPT and is published
unchanged, with its OpenAI Content Credentials (a C2PA manifest naming ChatGPT as the generator) intact. The
2026-11-02 checkers were written by ChatGPT/Codex. This repository was prepared, and the result re-verified, with
Claude (Anthropic) acting as a research operator under the author's direction. Three further from-scratch checks were
written by separate Claude agent instances that did not use the existing verifier code (see
`results/VERIFICATION.md`). No part of the COMPLEX system is needed to check the result, and none is included. The
author is responsible for the content.

## Citation

If you use this result, please cite the archived record:

> Majchrzak, F. (2026). *Erdős Problem #203: no pure finite prime covering with period dividing lcm(1,...,17) =
> 12,252,240* (Version 1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.24050485

The DOI above always resolves to the latest version; version 1.0.0 specifically is https://doi.org/10.5281/zenodo.24050486.
GitHub's "Cite this repository" uses `CITATION.cff`.

## License

Article, data and documentation: CC BY 4.0 (`LICENSE`). Verifier code in `verify/`: MIT (`LICENSE-CODE`).
