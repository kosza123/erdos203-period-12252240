"""Independent from-scratch check of the theorem (does not use verify_dual_stdlib.py or data/pool_factors.json).

Consolidated from the from-scratch checks written on 2026-11-08 by separate Claude agent instances that had no
access to the other verifiers in this repository.  It recomputes everything from N alone:

  1. the pool by brute force: all primes 3 < p <= 3*10^7 with 2^N = 3^N = 1 (mod p)   [no GMP needed]
  2. completeness: G_N = gcd(2^N - 1, 3^N - 1) with GMP, and G_N divided by the pool primes leaves cofactor 1,
     so no pool prime exceeds the sieve limit                                         [needs GMP, see below]
  3. multiplicative orders (by factoring p - 1), e_p, i_p, d_p, the weights of equation (8), R_d, mu_d,
     Table 1 and B = 3961943/4084080 < 1, exactly with fractions
  4. Lemma 1 by brute force: for every coset of every pool prime with e_p <= 260, the proportion on each line
     k - l = tau and the pooled proportion on each class t (mod 60)
  5. the 480 divisors of N (pool inclusion keeps the bound below 1)
  6. the dual certificate of Appendix A, transcribed here from the article and compared with
     data/optional_lp_dual.json (optimality within the fixed cut only; the theorem does not need it)

GMP: gmpy2 if installed, otherwise the system GMP library through ctypes (use --gmp PATH if it is not found).
With --no-gmp, step 2 is skipped; the result then only covers primes up to 3*10^7 and the status is
PASS_PARTIAL_PRIMES_UP_TO_3E7 instead of PASS.

usage:  python3 verify/independent_check.py [--gmp PATH] [--no-gmp] [--output RESULT.json]
Runs in about 15 seconds.
"""
import argparse, ctypes, ctypes.util, json, math, os, time
from fractions import Fraction as F
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
N = 12252240
LIM = 30_000_000
results = {}


def record(name, ok, detail=""):
    results[name] = {"pass": bool(ok), "detail": detail}
    print(("PASS " if ok else "FAIL ") + name + (f" :: {detail}" if detail != "" else ""))
    return ok


def gcd_with_gmp(explicit):
    if explicit is None:
        try:
            import gmpy2
            return int(gmpy2.gcd(gmpy2.mpz(2) ** N - 1, gmpy2.mpz(3) ** N - 1)), "gmpy2"
        except ImportError:
            pass
    names = [explicit] if explicit else [ctypes.util.find_library("gmp"), "libgmp.so.10", "libgmp.10.dylib",
             "/opt/homebrew/lib/libgmp.dylib", "/usr/local/lib/libgmp.dylib", "libgmp-10.dll",
             r"C:\Program Files\Git\mingw64\bin\libgmp-10.dll", r"C:\msys64\mingw64\bin\libgmp-10.dll"]
    lib = None
    for name in filter(None, names):
        try:
            if os.name == "nt" and os.path.isabs(name) and os.path.exists(name):
                os.add_dll_directory(os.path.dirname(name))
            lib = ctypes.CDLL(name); break
        except OSError as err:
            if explicit:
                raise SystemExit(f"cannot load GMP from {name}: {err}")
            continue
    if lib is None:
        raise SystemExit("GMP not found: pip install gmpy2, pass --gmp PATH, or use --no-gmp")

    class MPZ(ctypes.Structure):
        _fields_ = [("alloc", ctypes.c_int), ("size", ctypes.c_int), ("d", ctypes.c_void_p)]
    P = ctypes.POINTER(MPZ)
    for fn, args in (("init", [P]), ("clear", [P]), ("ui_pow_ui", [P, ctypes.c_ulong, ctypes.c_ulong]),
                     ("sub_ui", [P, P, ctypes.c_ulong]), ("gcd", [P, P, P]), ("get_str", [ctypes.c_char_p, ctypes.c_int, P])):
        getattr(lib, "__gmpz_" + fn).argtypes = args
    lib.__gmpz_sizeinbase.argtypes = [P, ctypes.c_int]; lib.__gmpz_sizeinbase.restype = ctypes.c_size_t
    lib.__gmpz_get_str.restype = ctypes.c_void_p
    a, b, g = MPZ(), MPZ(), MPZ()
    for z in (a, b, g):
        lib.__gmpz_init(ctypes.byref(z))
    lib.__gmpz_ui_pow_ui(ctypes.byref(a), 2, N); lib.__gmpz_sub_ui(ctypes.byref(a), ctypes.byref(a), 1)
    lib.__gmpz_ui_pow_ui(ctypes.byref(b), 3, N); lib.__gmpz_sub_ui(ctypes.byref(b), ctypes.byref(b), 1)
    lib.__gmpz_gcd(ctypes.byref(g), ctypes.byref(a), ctypes.byref(b))
    buf = ctypes.create_string_buffer(lib.__gmpz_sizeinbase(ctypes.byref(g), 16) + 2)
    lib.__gmpz_get_str(buf, 16, ctypes.byref(g))
    for z in (a, b, g):
        lib.__gmpz_clear(ctypes.byref(z))
    return int(buf.value.decode(), 16), "GMP via ctypes"


def order(x, p):
    n, t, fs, q = p - 1, p - 1, [], 2
    while q * q <= t:
        if t % q == 0:
            fs.append(q)
            while t % q == 0:
                t //= q
        q += 1
    if t > 1:
        fs.append(t)
    m = n
    for q in fs:
        while m % q == 0 and pow(x, m // q, p) == 1:
            m //= q
    return m


def main(args):
    if args.output and Path(args.output).exists():
        raise SystemExit(f"{args.output} exists; choose another name")
    t0 = time.time()
    record("N = lcm(1..17)", N == math.lcm(*range(1, 18)), N)

    # 1. pool by brute force
    sieve = bytearray([1]) * (LIM + 1); sieve[0] = sieve[1] = 0
    for i in range(2, int(LIM ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, LIM + 1, i)))
    pool = [p for p in range(5, LIM + 1) if sieve[p] and pow(2, N, p) == 1 and pow(3, N, p) == 1]
    record("pool primes up to 3e7: 192, largest 8168161", len(pool) == 192 and max(pool) == 8168161,
           f"{len(pool)} primes, largest {max(pool)}")

    # 2. completeness with GMP
    if args.no_gmp:
        print("NOTE completeness (GMP gcd) skipped by --no-gmp: the result covers primes <= 3e7 only")
        backend = None
    else:
        G, backend = gcd_with_gmp(args.gmp)
        co, val = G, {}
        for p in pool:
            v = 0
            while co % p == 0:
                co //= p; v += 1
            val[p] = v
        record("G_N has 2432 bits and 732 decimal digits", G.bit_length() == 2432 and len(str(G)) == 732,
               (G.bit_length(), len(str(G))))
        record("2^N = 3^N = 1 mod G_N", pow(2, N, G) == 1 and pow(3, N, G) == 1)
        record("every pool prime divides G_N", all(v >= 1 for v in val.values()))
        record("cofactor of G_N after the pool primes is 1 (pool complete)", co == 1, co)
        g2 = G
        for p in pool[:-1]:
            while g2 % p == 0:
                g2 //= p
        record("control: dropping the largest pool prime leaves cofactor 8168161", g2 == 8168161, g2)
        record("valuation 2 exactly for 5, 7, 11, 13, 17", {p: v for p, v in val.items() if v != 1} == {5: 2, 7: 2, 11: 2, 13: 2, 17: 2})

    # 3. orders, weights, Table 1, B
    info = {}
    for p in pool:
        o2, o3, o6 = order(2, p), order(3, p), order(6, p)
        e = math.lcm(o2, o3)
        info[p] = dict(e=e, h=o6, i=e // o6, d=gcd(e // o6, 60))
    record("e_p | N and ord_p(6) | e_p for all pool primes", all(N % v["e"] == 0 and v["e"] % v["h"] == 0 for v in info.values()))
    record("p=5: e=4, i=4, d=4; p=7: e=6, h=2, i=3, d=3",
           (info[5]["e"], info[5]["i"], info[5]["d"]) == (4, 4, 4) and (info[7]["e"], info[7]["h"], info[7]["i"], info[7]["d"]) == (6, 2, 3, 3))
    lam = {t: F(1, 20) for t in (1, 26, 41, 46)} | {t: F(1, 10) for t in (2, 5, 7, 10, 14, 19, 23, 58)}
    record("weights: nonnegative, sum 1, avoid 0 mod 3 and 0 mod 4",
           sum(lam.values()) == 1 and all(w >= 0 and t % 3 and t % 4 for t, w in lam.items()))
    data = json.loads((ROOT / "data" / "optional_lp_dual.json").read_bytes())
    record("weights equal those in data/optional_lp_dual.json", {int(t): F(w) for t, w in data["weights"].items()} == lam)
    R = {}
    for p in pool:
        if p not in (5, 7):
            R[info[p]["d"]] = R.get(info[p]["d"], F(0)) + F(1, info[p]["e"])

    def mu(d):
        return max(sum((w for t, w in lam.items() if t % d == c), F(0)) for c in range(d))
    B = sum(d * mu(d) * R[d] for d in R)
    table = {1: ("1", "319967/556920", 7039274), 2: ("1/2", "53617/291720", 2251914), 3: ("1/2", "120079/2042040", 1080711),
             4: ("1/2", "91723/4084080", 550338), 5: ("1/5", "92789/2450448", 463945), 6: ("1/4", "10343/875160", 217203),
             10: ("1/10", "16/36465", 5376), 12: ("1/4", "181/32760", 203082), 15: ("1/10", "61/16380", 68442), 20: ("1/10", "1/4420", 5544)}
    record("Table 1 of the article (mu_d, R_d, N d mu_d R_d) reproduced exactly",
           set(table) == set(R) and all(F(m) == mu(d) and F(r) == R[d] and c == N * d * mu(d) * R[d] for d, (m, r, c) in table.items()))
    record("B = 3961943/4084080 < 1", B == F(3961943, 4084080) and B < 1, f"{B} = {float(B):.6f}, gap {1 - B}")
    c5 = F(info[5]["d"], info[5]["e"]) * sum((w for t, w in lam.items() if t % info[5]["d"] == 0), F(0))
    c7 = F(info[7]["d"], info[7]["e"]) * sum((w for t, w in lam.items() if t % info[7]["d"] == 0), F(0))
    record("normalised 5 (a_5 = 0 mod 4) and 7 (a_7 = 0 mod 3) contribute 0", c5 == 0 and c7 == 0, (str(c5), str(c7)))
    uniform = {t: F(1, 30) for t in range(60) if t % 3 and t % 4}
    Bu = sum(d * max(sum((w for t, w in uniform.items() if t % d == c), F(0)) for c in range(d)) * R[d] for d in R)
    record("control: uniform weights on the 30 admissible classes do NOT give B < 1 for this pool", Bu >= 1, f"{Bu} = {float(Bu):.4f}")

    # 4. Lemma 1 by brute force
    small = [p for p in pool if info[p]["e"] <= 260]
    ok = True
    for p in small:
        e, h, i, d = (info[p][k] for k in ("e", "h", "i", "d"))
        M = math.lcm(e, 60)
        pow6 = [pow(6, u, p) for u in range(e)]
        counts = []
        for tau in range(M):
            c, base = {}, pow(2, tau, p)
            for u in range(e):
                v = base * pow6[u] % p
                c[v] = c.get(v, 0) + 1
            counts.append(c)
        H = {pow(2, k, p) * pow(3, l, p) % p for k in range(e) for l in range(e)}
        ok &= len(H) == e
        for g in H:
            prop = [F(counts[tau].get(g, 0), e) for tau in range(M)]
            a = next(tau for tau in range(M) if prop[tau]) % i
            ok &= all(prop[tau] == (F(i, e) if tau % i == a else 0) for tau in range(M))
            for t in range(60):
                lines = range(t, M, 60)
                ok &= sum((prop[tau] for tau in lines), F(0)) / len(lines) == (F(d, e) if t % d == a % d else 0)
    record(f"Lemma 1 by brute force for all cosets of the {len(small)} pool primes with e_p <= 260", ok)

    # 5. divisors
    divs = [x for x in range(1, N + 1) if N % x == 0]
    worst = max(sum((F(info[p]["d"], info[p]["e"]) * mu(info[p]["d"]) for p in pool if p not in (5, 7) and D % info[p]["e"] == 0), F(0)) for D in divs)
    record("480 divisors of N, each with bound <= B < 1", len(divs) == 480 and worst <= B, f"max bound {worst}")

    # 6. dual certificate, transcribed from Appendix A of the article
    dual = [(1,0,"319967/556920"),(2,0,"16519/157080"),(2,1,"31523/120120"),(3,1,"27697/291720"),(3,2,"83179/1021020"),
            (4,2,"91723/1021020"),(5,0,"89597/2450448"),(5,1,"7199/188496"),(5,2,"7199/188496"),(5,3,"7199/188496"),
            (5,4,"7199/188496"),(6,2,"28019/408408"),(6,5,"1569/680680"),(10,0,"161/58344"),(10,5,"95/58344"),
            (12,10,"181/2730"),(15,2,"61/5460"),(15,5,"61/5460"),(15,8,"61/5460"),(15,11,"61/5460"),(15,14,"61/5460"),
            (20,2,"1/884"),(20,6,"1/884"),(20,14,"1/884"),(20,18,"1/884")]
    record("Appendix A equals data/optional_lp_dual.json",
           sorted((a, b, F(c)) for a, b, c in data["dual"]) == sorted((a, b, F(c)) for a, b, c in dual))
    q = {(d, c): F(v) for d, c, v in dual}
    S = [t for t in range(60) if t % 3 and t % 4]
    record("dual feasible: q >= 0, budgets sum_c q_dc <= d R_d, all 30 inequalities >= B",
           all(v >= 0 for v in q.values())
           and all(sum((v for (dd, c), v in q.items() if dd == d), F(0)) <= d * R.get(d, F(0)) for d in range(1, 61) if 60 % d == 0)
           and len(S) == 30 and all(sum((v for (d, c), v in q.items() if t % d == c), F(0)) >= B for t in S))

    status = "FAIL" if not all(r["pass"] for r in results.values()) else ("PASS_PARTIAL_PRIMES_UP_TO_3E7" if args.no_gmp else "PASS")
    out = {"status": status, "scope": "pool <= 3e7 only" if args.no_gmp else "complete", "gmp_backend": backend,
           "B": str(B), "checks": results, "elapsed_seconds": round(time.time() - t0, 1)}
    print(status, f"({out['elapsed_seconds']} s)")
    if args.output:
        with open(args.output, "x", encoding="utf-8", newline="\n") as f:
            json.dump(out, f, indent=2); f.write("\n")
    if status == "FAIL":
        raise SystemExit(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gmp", default=None, help="path or name of the GMP shared library (skips gmpy2)")
    ap.add_argument("--no-gmp", action="store_true", help="skip the gcd; checks primes up to 3*10^7 only")
    ap.add_argument("--output", default=None, help="write JSON here (must not exist)")
    main(ap.parse_args())
