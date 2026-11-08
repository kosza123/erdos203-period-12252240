"""Which of the 480 divisors of N = lcm(1,...,17) are new exclusions, and why each one is excluded.

Prior exclusions are taken from the Erdős Problem a Day post of 14 October 2026
(https://erdosproblemaday.com/day/203-the-quotient-lemma-lines-close-the-open-pools), which states:
"No covering system has period N for any N <= 65,520, nor for any divisor of 110880, 166320, 221760, 332640 or
720720 = lcm(1, ..., 16)" (pure coverings).  A divisor T of N is counted as previously excluded if T <= 65520 or T
divides one of those five numbers; all other divisors are listed as new.

For every new T the script computes, from data/pool_factors.json (whose completeness is proved by
verify_pool_gmp.py):
  pool_size   number of primes p > 3 with e_p | T  (the pool of T, a subset of the pool of N)
  density     sum of 1/e_p over that pool (a covering needs density >= 1, AnimishSharma's condition)
  B_article   the line bound with the weights of the article, eq. (8)  (always <= B = 3961943/4084080 < 1)
  B_uniform   the same bound with uniform weights 1/30 on the 30 classes mod 60 avoiding 0 mod 3 and 0 mod 4
  reason      the weakest argument that already excludes T:
                density<1          the pool is too thin for any covering
                uniform_weights    the 5/7 cut of the 14 October post with uniform weights gives B < 1
                article_weights    only the non-uniform weights of the article give B < 1
Every new T is excluded by the theorem; the reason column only shows how much of the argument is needed.

usage:  python3 verify/new_periods.py [--output NEW_PERIODS.csv]
"""
import argparse, json, math
from fractions import Fraction as F
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
N = 12252240
PRIOR = (110880, 166320, 221760, 332640, 720720)
SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17)


def order(a, p):
    o = N
    for q in SMALL_PRIMES:
        while o % q == 0 and pow(a, o // q, p) == 1:
            o //= q
    assert pow(a, o, p) == 1
    return o


def factor_string(T):
    parts = []
    for q in SMALL_PRIMES:
        k = 0
        while T % q == 0:
            T //= q; k += 1
        if k:
            parts.append(f"{q}^{k}" if k > 1 else str(q))
    return "*".join(parts)


def main(args):
    if args.output and Path(args.output).exists():
        raise SystemExit(f"{args.output} exists; choose another name")
    assert N == math.lcm(*range(1, 18))
    pairs = json.loads((ROOT / "data" / "pool_factors.json").read_bytes())["factors"]
    info = {}
    for p, _ in pairs:
        assert p > 3 and pow(2, N, p) == 1 and pow(3, N, p) == 1
        e = math.lcm(order(2, p), order(3, p))
        info[p] = (e, gcd(e // order(6, p), 60))
    article = {t: F(1, 20) for t in (1, 26, 41, 46)} | {t: F(1, 10) for t in (2, 5, 7, 10, 14, 19, 23, 58)}
    uniform = {t: F(1, 30) for t in range(60) if t % 3 and t % 4}

    def mu(w, d):
        return max(sum((x for t, x in w.items() if t % d == c), F(0)) for c in range(d))

    def bound(w, T):
        return sum((F(d, e) * mu(w, d) for p, (e, d) in info.items() if p not in (5, 7) and T % e == 0), F(0))

    divisors = [T for T in range(1, N + 1) if N % T == 0]
    prior = [T for T in divisors if T <= 65520 or any(x % T == 0 for x in PRIOR)]
    new = [T for T in divisors if T not in set(prior)]
    rows, counts = [], {}
    for T in new:
        pool = [p for p, (e, _) in info.items() if T % e == 0]
        dens = sum((F(1, info[p][0]) for p in pool), F(0))
        ba, bu = bound(article, T), bound(uniform, T)
        reason = "density<1" if dens < 1 else ("uniform_weights" if bu < 1 else "article_weights")
        assert ba < 1
        counts[reason] = counts.get(reason, 0) + 1
        rows.append((T, factor_string(T), len(pool), dens, ba, bu, reason))
    print(f"divisors of N: {len(divisors)}; previously excluded: {len(prior)}; new: {len(new)} "
          f"(from {min(new)} to {max(new)}; all multiples of 17: {all(T % 17 == 0 for T in new)})")
    print("new periods by weakest sufficient argument:", counts)
    worst = max((r for r in rows if r[6] != "density<1" and r[0] != N), key=lambda r: r[5], default=None)
    if worst:
        print(f"largest uniform-weight bound among proper divisors with density >= 1: {float(worst[5]):.4f} at T = {worst[0]}")
    lines = ["T,factorisation,pool_size,density,density_decimal,B_article,B_article_decimal,B_uniform,B_uniform_decimal,reason"]
    for T, fs, n, dens, ba, bu, reason in rows:
        lines.append(f"{T},{fs},{n},{dens},{float(dens):.6f},{ba},{float(ba):.6f},{bu},{float(bu):.6f},{reason}")
    text = "\n".join(lines) + "\n"
    if args.output:
        with open(args.output, "x", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print("written", args.output)
    else:
        print(text, end="")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", default=None, help="write the CSV here (must not exist); default: print it")
    main(ap.parse_args())
