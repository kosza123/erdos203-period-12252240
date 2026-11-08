"""Completeness of the prime pool: direct gcd(2^N - 1, 3^N - 1) with GMP, compared with the factor list.

For N = lcm(1,...,17) = 12252240 a prime p > 3 satisfies e_p | N exactly when 2^N = 3^N = 1 (mod p), i.e. when
p divides G_N = gcd(2^N - 1, 3^N - 1).  This script computes G_N directly (2^N - 1 has 12,252,240 bits) and checks
that it equals the product of p^v over data/pool_factors.json, and that it equals data/G_N_decimal.txt.
Together with the primality of the listed factors (checked by verify_dual_stdlib.py) this proves that the list
is the complete pool.  Perturbation controls: dropping a listed prime, raising a valuation, or changing one digit of
the decimal must each break the equality.

GMP backends, tried in order: gmpy2 (pip install gmpy2), then the system GMP library through ctypes
(libgmp.so.10 on Linux, libgmp.10.dylib on macOS, libgmp-10.dll on Windows, e.g. from Git for Windows' mingw64
or MSYS2).  Pass --gmp PATH to name the library file explicitly.

Based on the 2026-11-02 GMP checker (ChatGPT/Codex); this public version adds the gmpy2 backend, library
discovery on Windows and macOS, and the decimal-file comparison.

usage:  python3 verify/verify_pool_gmp.py [--gmp PATH] [--output RESULT.json]
"""
import argparse
import ctypes
import ctypes.util
import hashlib
import json
import os
import platform
import time
from math import lcm, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CANDIDATES = ["libgmp.so.10", "libgmp.so", "libgmp.10.dylib", "/opt/homebrew/lib/libgmp.dylib",
              "/usr/local/lib/libgmp.dylib", "libgmp-10.dll", r"C:\Program Files\Git\mingw64\bin\libgmp-10.dll",
              r"C:\msys64\mingw64\bin\libgmp-10.dll", r"C:\msys64\ucrt64\bin\libgmp-10.dll"]


def gcd_gmpy2(N):
    import gmpy2
    a = gmpy2.mpz(2) ** N - 1
    b = gmpy2.mpz(3) ** N - 1
    return int(gmpy2.gcd(a, b)), "gmpy2 " + gmpy2.version() + " / GMP " + gmpy2.mp_version()


class MPZ(ctypes.Structure):
    _fields_ = [("alloc", ctypes.c_int), ("size", ctypes.c_int), ("limbs", ctypes.c_void_p)]


def load_gmp(explicit):
    names = [explicit] if explicit else ([ctypes.util.find_library("gmp")] + CANDIDATES)
    for name in names:
        if not name:
            continue
        try:
            if os.name == "nt" and os.path.isabs(name) and os.path.exists(name):
                os.add_dll_directory(os.path.dirname(name))
            return ctypes.CDLL(name), name
        except OSError as err:
            if explicit:
                raise SystemExit(f"cannot load GMP from {name}: {err}")
            continue
    raise SystemExit("GMP not found: install gmpy2 (pip install gmpy2) or pass --gmp PATH-TO-LIBGMP")


def gcd_ctypes(N, explicit):
    gmp, name = load_gmp(explicit)
    ptr = ctypes.POINTER(MPZ)

    def bind(fn, args, res=None):
        f = getattr(gmp, "__gmpz_" + fn)
        f.argtypes, f.restype = args, res
        return f

    init, clear = bind("init", [ptr]), bind("clear", [ptr])
    pow_ui = bind("ui_pow_ui", [ptr, ctypes.c_ulong, ctypes.c_ulong])
    sub_ui = bind("sub_ui", [ptr, ptr, ctypes.c_ulong])
    gcd = bind("gcd", [ptr, ptr, ptr])
    sizeinbase = bind("sizeinbase", [ptr, ctypes.c_int], ctypes.c_size_t)
    get_str = bind("get_str", [ctypes.c_char_p, ctypes.c_int, ptr], ctypes.c_void_p)
    a, b, g = MPZ(), MPZ(), MPZ()
    for x in (a, b, g):
        init(ctypes.byref(x))
    try:
        pow_ui(ctypes.byref(a), 2, N); sub_ui(ctypes.byref(a), ctypes.byref(a), 1)
        pow_ui(ctypes.byref(b), 3, N); sub_ui(ctypes.byref(b), ctypes.byref(b), 1)
        gcd(ctypes.byref(g), ctypes.byref(a), ctypes.byref(b))
        buf = ctypes.create_string_buffer(sizeinbase(ctypes.byref(g), 10) + 2)
        get_str(buf, 10, ctypes.byref(g))
        value = int(buf.value.decode())
    finally:
        for x in (a, b, g):
            clear(ctypes.byref(x))
    version = ctypes.c_char_p.in_dll(gmp, "__gmp_version").value.decode()
    return value, f"GMP {version} via ctypes ({name})"


def main(args):
    if args.output is not None and args.output.exists():
        raise SystemExit(f"{args.output} exists; choose another name")
    t0 = time.perf_counter()
    factors = json.loads(args.factors.read_bytes())
    N = factors["N"]
    if N != lcm(*range(1, 18)):
        raise ValueError("unexpected N")
    product = prod(p ** v for p, v in factors["factors"])
    if args.gmp is None:
        try:
            G, backend = gcd_gmpy2(N)
        except ImportError:
            G, backend = gcd_ctypes(N, None)
    else:
        G, backend = gcd_ctypes(N, args.gmp)
    decimal_text = args.decimal.read_text().strip()
    decimal_file = int(decimal_text)
    pairs = factors["factors"]
    dropped = prod(p ** v for p, v in pairs[:-1])                      # drop the largest listed prime
    raised = product * pairs[0][0]                                      # raise the valuation of the smallest prime
    flipped = int(decimal_text[:-1] + str((int(decimal_text[-1]) + 1) % 10))   # change the last digit
    checks = {
        "gcd_equals_factor_product": G == product,
        "gcd_equals_G_N_decimal_txt": G == decimal_file,
        "control_dropped_prime_rejected": G != dropped,
        "control_raised_valuation_rejected": G != raised,
        "control_changed_digit_rejected": G != flipped,
    }
    digits = len(str(G))
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "N": N, "checks": checks, "distinct_primes_in_list": len(factors["factors"]),
        "largest_prime": max(p for p, _ in factors["factors"]),
        "gcd_bits": G.bit_length(), "gcd_exact_decimal_digits": digits,
        "gcd_decimal_sha256": hashlib.sha256(str(G).encode()).hexdigest(),
        "factor_file_sha256": hashlib.sha256(args.factors.read_bytes()).hexdigest(),
        "backend": backend, "platform": platform.platform(),
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
    }
    text = json.dumps(result, indent=2)
    if args.output is None:
        print(text)
    else:
        with args.output.open("x", encoding="utf-8", newline="\n") as out:
            out.write(text + "\n")
        print(result["status"], "bits", result["gcd_bits"], "digits", digits)
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--factors", type=Path, default=ROOT / "data" / "pool_factors.json")
    ap.add_argument("--decimal", type=Path, default=ROOT / "data" / "G_N_decimal.txt")
    ap.add_argument("--gmp", default=None, help="path or name of the GMP shared library (skips gmpy2)")
    ap.add_argument("--output", type=Path, default=None, help="write JSON here (must not exist); default: print")
    main(ap.parse_args())
