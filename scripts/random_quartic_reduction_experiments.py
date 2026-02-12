# random_quartic_reduction_experiments.py
#
# Experimental script for sampling random smooth plane quartics over QQ and
# computing their p-adic stable reduction types.
#
# The script generates random homogeneous quartic forms with integral
# coefficients, filters for smooth curves over QQ, computes the geometric
# stable reduction at a chosen p-adic valuation, and groups the results by
# reduction type (including the cases "hyperelliptic" and "fail").
#
# It is intended for exploratory computations and data generation, not as a
# library module or automated test.

from collections import defaultdict
from random import Random

from sage.all import QQ, PolynomialRing, Curve
from semistable_model.curves.stable_reduction_of_quartics import stable_reduction_of_quartic


# ------------------------------------------------------------
# Random homogeneous quartics over QQ
# ------------------------------------------------------------


def random_int(rng, bound):
    """
    Sample an integer with typical size ~ bound.
    Uses a symmetric uniform distribution in [-bound, bound].
    """
    if bound <= 0:
        return 0
    return rng.randint(-bound, bound)


def random_quartic_form_QQ(
    *,
    seed=None,
    n_terms=8,
    coeff_bound=10
):
    r"""
    Return a random homogeneous quartic F(x,y,z) over QQ.

    Parameters:
      - n_terms: number of monomials used (<= 15 for quartics)
      - coeff_bound: coefficients are sampled from [-coeff_bound, coeff_bound]
      - ensure_nonzero: ensure F != 0
    """
    rng = Random(seed)
    R = PolynomialRing(QQ, names=("x", "y", "z"))
    mons = R.monomials_of_degree(4) # 15 monomials
    rng.shuffle(mons)
    mons = mons[: min(n_terms, len(mons))]

    F = R.zero()
    for m in mons:
        c = random_int(rng, coeff_bound)
        # allow zeros but try to avoid too many trivial terms
        if c == 0:
            continue
        F += QQ(c) * m

    if F == 0:
        # force one term
        F = QQ(1) * mons[0]

    return F


# ------------------------------------------------------------
# Experiment loop
# ------------------------------------------------------------

def experiment_random_quartics(
    *,
    n_samples=10,
    seed=0,
    n_terms=8,
    coeff_bound=10,
    v_K=None,
    max_tries_factor=50,
    verbose=True,
):
    r"""
    Generate random smooth quartics over QQ, compute stable reduction types,
    and return a dict: bucket -> list of StableReductionResult.

    Buckets include reduction types, plus "hyperelliptic" and "fail".

    Parameters:
      - n_samples: number of *smooth* quartics to process
      - seed: RNG seed
      - n_terms, coeff_bound: random generator knobs
      - v_K: a p-adic valuation on QQ (required)
      - max_tries_factor: max attempts = max_tries_factor * n_samples
      - verbose: print progress
    """
    if v_K is None:
        raise ValueError("Please pass v_K, e.g. v_K = QQ.valuation(2).")

    rng = Random(seed)

    buckets = defaultdict(list)

    tries = 0
    found = 0
    max_tries = max_tries_factor * n_samples

    while found < n_samples and tries < max_tries:
        tries += 1

        # vary seed per attempt to avoid repetition
        F = random_quartic_form_QQ(
            seed=rng.randint(0, 10**18),
            n_terms=n_terms,
            coeff_bound=coeff_bound
        )

        if not Curve(F).is_smooth():
            continue

        found += 1
        if verbose:
            print(f"[{found}/{n_samples}] smooth quartic found; running stable reduction...")
            print(f"F = {F}")
        res = stable_reduction_of_quartic(F, v_K)

        # bucket key
        if res.status == "ok":
            key = res.reduction_type
        else:
            key = res.status  # "hyperelliptic" or "fail"

        buckets[key].append(res)

        if verbose:
            print(f"  -> status={res.status}, type={getattr(res, 'reduction_type', None)}")

    if verbose:
        print("")
        print(f"Generated {found} smooth quartics in {tries} attempts.")
        print("Bucket sizes:")
        for k in sorted(buckets.keys()):
            print(f"  {k}: {len(buckets[k])}")

    return buckets


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
""""
if __name__ == "__main__":
    # Example: stable reduction at p=2 over QQ
    v2 = QQ.valuation(2)

    buckets = experiment_random_quartics(
        n_samples=10,
        v_K=v2
    )

    # Example: inspect one result from a bucket
    for k, lst in buckets.items():
        if lst:
            print("")
            print("Example from bucket:", k)
            print(lst[0])
            break
"""