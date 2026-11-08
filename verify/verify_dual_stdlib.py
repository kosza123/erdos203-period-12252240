"""Exact checker of the corrected article. Python standard library only.

Public copy of the checker written on 2026-11-02 (ChatGPT/Codex, within the author's AI-assisted workflow) for the
article in paper/. Changes in this copy: this header (the original first line called the checker "independent",
which is too strong), input defaults point to data/ and paper/ of this repository, --output is optional
(the result is printed as JSON only when it is omitted), the output file always uses LF line endings, and input
paths outside the repository are allowed.
The mathematical checks are unchanged.

usage:  python3 verify/verify_dual_stdlib.py [--output RESULT.json]
"""
import argparse
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parent.parent


def rel(path):
    try:
        return str(path.resolve().relative_to(ROOT).as_posix())
    except ValueError:
        return str(path)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sieve(limit):
    flags = bytearray(b'\1') * (limit + 1)
    flags[:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p:limit+1:p] = b'\0' * ((limit-p*p)//p+1)
    return [p for p, flag in enumerate(flags) if flag]


def order_from_N(base, prime, N, divisors):
    require(pow(base, N, prime) == 1, 'Prime outside declared pool')
    order = N
    for q in divisors:
        while order % q == 0 and pow(base, order // q, prime) == 1:
            order //= q
    require(pow(base, order, prime) == 1, 'Bad order witness')
    require(all(pow(base, order//q, prime) != 1 for q in divisors if order % q == 0),
            'Order is not minimal')
    return order


def check_pool(factors, data):
    N, modulus = data['N'], data['modulus']
    require(N == lcm(*range(1, 18)) and modulus == 60, 'Wrong target')
    pairs = factors['factors']
    require(factors['N'] == N and len(pairs) == 192, 'Wrong factor count')
    require(all(type(p) is int and type(v) is int and p > 3 and v > 0 for p, v in pairs), 'Malformed factors')
    primes = [p for p, _ in pairs]
    require(len(set(primes)) == len(primes), 'Duplicate factors')
    trial_primes = sieve(isqrt(max(primes)))
    for p in primes:
        require(all(p == q or p % q for q in trial_primes if q*q <= p), 'Composite factor')
    rows = []
    for p in primes:
        orders = [order_from_N(a, p, N, [2, 3, 5, 7, 11, 13, 17]) for a in (2, 3, 6)]
        e, h = lcm(*orders[:2]), orders[2]
        require(e % h == 0 and N % e == 0, 'Inconsistent subgroup size')
        rows.append({'p': p, 'ord2': orders[0], 'ord3': orders[1], 'e': e, 'h': h,
                     'i': e//h, 'd': gcd(e//h, modulus)})
    value = prod(p**v for p, v in pairs)
    require(max(primes) == 8168161 and value.bit_length() == 2432 and len(str(value)) == 732,
            'Article size claim mismatch')
    require({p: v for p, v in pairs if v != 1} == {5: 2, 7: 2, 11: 2, 13: 2, 17: 2},
            'Article multiplicity claim mismatch')
    return rows, value


def check_dual(data, groups):
    M, bound = data['modulus'], Fraction(data['bound'])
    support = [t for t in range(M) if t % 3 and t % 4]
    entries = data['dual']
    q = {(d, c): Fraction(value) for d, c, value in entries}
    require(len(q) == len(entries), 'Duplicate dual entries')
    require(all(d in groups and 0 <= c < d and v >= 0 for (d,c),v in q.items()), 'Bad dual entry')
    budgets = {}
    for d, R in groups.items():
        total = sum((value for (dd,c),value in q.items() if dd == d), Fraction())
        require(total <= d*R, 'Dual budget exceeded')
        budgets[d] = str(d*R-total)
    slack = {}
    for t in support:
        total = sum((q.get((d,t%d), Fraction()) for d in groups), Fraction())
        require(total >= bound, 'Dual lower bound violated')
        slack[t] = str(total-bound)
    return {'nonzero_entries': len(q), 'budget_inequalities': len(budgets),
            'support_inequalities': len(slack), 'budget_slack': budgets, 'support_slack': slack}


def check_primal(data, rows):
    M, N = data['modulus'], data['N']
    weights = {int(t): Fraction(value) for t, value in data['weights'].items()}
    require(sum(weights.values(), Fraction()) == 1, 'Weights do not sum to one')
    require(all(0 <= t < M and t%3 and t%4 and v >= 0 for t,v in weights.items()), 'Invalid support')
    groups = defaultdict(Fraction)
    for row in rows:
        if row['p'] not in (5,7):
            groups[row['d']] += Fraction(1,row['e'])
    claimed = {d: (Fraction(mu),Fraction(R),term) for d,mu,R,term in data['groups']}
    require(set(claimed) == set(groups), 'Wrong nonzero groups')
    calculated = []
    bound = Fraction()
    for d, R in sorted(groups.items()):
        mu = max(sum((v for t,v in weights.items() if t%d == c), Fraction()) for c in range(d))
        term = d*mu*R
        require((mu,R,N*term) == claimed[d], 'Article Table 1 mismatch')
        bound += term
        calculated.append([d,str(mu),str(R),int(N*term)])
    require(bound == Fraction(data['bound']) == Fraction(3961943,4084080) < 1, 'Wrong obstruction')
    return dict(groups), calculated


def check_projection(rows):
    # Enumerate actual residues for one least-cost representative of every nonzero d group.
    selected = {p: next(row for row in rows if row['p'] == p) for p in (5,7)}
    for d in {row['d'] for row in rows}:
        row = min((r for r in rows if r['d'] == d), key=lambda r:r['e'])
        selected[row['p']] = row
    examined = []
    wrong_d_detected = False
    for row in selected.values():
        p,h,i,d = (row[k] for k in ('p','h','i','d'))
        L = lcm(i,60)
        seen = defaultdict(lambda: defaultdict(int))
        power2 = 1
        for tau in range(L):
            residue = power2
            for u in range(h):
                seen[residue][tau%60] += 1
                residue = residue*6 % p
            power2 = power2*2 % p
        require(len(seen) == row['e'], 'Wrong image group size')
        for histogram in seen.values():
            a = min(histogram) % d
            expected = {t:1 for t in range(60) if t%d == a}
            require(dict(histogram) == expected, 'Direct projection enumeration disagrees')
            if d > 1 and dict(histogram) != {t:1 for t in range(60)}:
                wrong_d_detected = True
        examined.append({'prime':p,'e':row['e'],'d':d,'direct_modular_pairs':L*h})
    require(wrong_d_detected, 'Deliberately incorrect projection was not rejected')
    for a5 in range(4):
        for a7 in range(3):
            shifts = [s for s in range(12) if (a5+s)%4 == 0 and (a7+s)%3 == 0]
            require(len(shifts)==1, 'CRT normalization failed')
    return {'representatives':examined,'normalization_pairs':12,'wrong_projection_negative_control':'REJECTED'}


def negative_controls(data, groups):
    bad = deepcopy(data)
    bad['dual'] = [row for row in bad['dual'] if row[0] != 1]
    try:
        check_dual(bad, groups)
    except ValueError:
        return {'deleted_dual_mass':'REJECTED'}
    raise ValueError('Invalid dual accepted')


def main(args):
    wall, cpu = time.perf_counter(), time.process_time()
    data = json.loads(args.data.read_bytes())
    factors = json.loads(args.factors.read_bytes())
    require(digest(args.pdf) == data['source_pdf_sha256'], 'PDF pin mismatch')
    rows, product = check_pool(factors, data)
    groups, table = check_primal(data, rows)
    dual = check_dual(data, groups)
    geometry = check_projection(rows)
    negative = negative_controls(data, groups)
    divisors = [d for d in range(1,data['N']+1) if data['N'] % d == 0]
    require(len(divisors) == 480, 'Wrong divisor count')
    require(all(sum((Fraction(row['d'],row['e']) * max(sum((Fraction(v) for t,v in data['weights'].items() if int(t)%row['d']==c),Fraction()) for c in range(row['d'])) for row in rows if row['p'] not in (5,7) and D%row['e']==0),Fraction()) <= Fraction(data['bound']) for D in divisors), 'Divisor monotonicity failure')
    result = {'status':'PASS_EXACT_ARITHMETIC_AND_FINITE_PROJECTION_CHECKS', 'source_pdf_sha256':digest(args.pdf),
              'inputs_sha256':{rel(p):digest(p) for p in (args.data,args.factors)},
              'verifier_sha256':digest(Path(__file__)), 'N':data['N'],'distinct_proven_primes':len(rows),
              'max_prime':max(r['p'] for r in rows),'factor_product_bits':product.bit_length(),
              'factor_product_exact_decimal_digits':len(str(product)),'bound':data['bound'],
              'gap':str(1-Fraction(data['bound'])),'divisors_checked':len(divisors),'table':table,
              'dual':dual,'projection':geometry,'negative_controls':negative,
              'gcd_completeness':'REQUIRES_SEPARATE_DIRECT_GCD_CHECK',
              'elapsed_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,
              'model_api_calls':0}
    if args.output is None:
        print(json.dumps(result,indent=2))
    else:
        with args.output.open('x',encoding='utf-8',newline='\n') as output:
            json.dump(result,output,indent=2);output.write('\n')
        print(result['status'],result['bound'],result['gap'],result['elapsed_seconds'])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=ROOT/'data'/'optional_lp_dual.json')
    parser.add_argument('--factors',type=Path,default=ROOT/'data'/'pool_factors.json')
    parser.add_argument('--pdf',type=Path,default=ROOT/'paper'/'Erdos203_Article_Draft_2026-11-01.pdf')
    parser.add_argument('--output',type=Path,default=None,help='write JSON here (must not exist); default: print')
    main(parser.parse_args())
