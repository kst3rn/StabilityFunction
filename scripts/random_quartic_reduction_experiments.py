# random_quartic_reduction_experiments.py
#
# Experimental script for sampling random smooth plane quartics over QQ and
# computing their p-adic stable reduction types.
#
# The script generates random homogeneous quartic forms with integral
# coefficients, filters for smooth curves over QQ, computes the geometric
# stable reduction at a chosen p-adic valuation, and stores successful
# examples (status == "ok") in an append-only JSONL database.
#
# Intended for exploratory computations and data generation (not a test).
#
# Usage
# -----
# Typical invocation (example for p = 5):
#
#   sage -python random_quartic_reduction_experiments.py \
#       --n-samples 500 \
#       --prime 5 \
#       --seed 12345 \
#       --out /path/to/quartic_reduction_data/data/v1/p5.jsonl \
#       --stats /path/to/quartic_reduction_data/data/v1/p5_stats.json
#
# This:
#   - generates random smooth plane quartics over QQ,
#   - computes their geometric stable reduction at p,
#   - appends successful cases (status == "ok") to a JSONL file,
#   - writes a JSON stats summary (including status/type counts).
#
# The raw data should be stored in a separate data repository
# (e.g. “quartic_reduction_data”), not in the main code repository.
#
# The JSONL file contains one record per successful example.
# The stats file contains aggregate counts (ok/hyperelliptic/fail)
# and reduction-type frequencies.
#
# Run with --help for the full list of options.
#


import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
from random import Random

from sage.all import QQ, PolynomialRing, Curve
from semistable_model.curves.stable_reduction_of_quartics import stable_reduction_of_quartic


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def timed_call(func, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    end = time.perf_counter()
    return result, end - start


def append_jsonl(path, record):
    """
    Append a JSON-serializable dict as one line to a .jsonl file.
    Creates parent directories if needed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True) + "\n")
        f.flush()


def write_run_stats(path, stats):
    """
    Write a small JSON stats file for the run (overwritten).
    Creates parent directories if needed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(stats, f, sort_keys=True, indent=2)


# ------------------------------------------------------------
# Random homogeneous quartics over QQ
# ------------------------------------------------------------

def random_int(rng, bound):
    """Sample an integer uniformly from [-bound, bound]."""
    if bound <= 0:
        return 0
    return rng.randint(-bound, bound)


def random_quartic_form_QQ(*, rng, n_terms=8, coeff_bound=10):
    r"""
    Return a random homogeneous quartic F(x,y,z) over QQ.

    Parameters:
      - n_terms: number of monomials used (<= 15 for quartics)
      - coeff_bound: coefficients are sampled from [-coeff_bound, coeff_bound]
    """
    R = PolynomialRing(QQ, names=("x", "y", "z"))
    mons = list(R.monomials_of_degree(4))  # 15 monomials
    rng.shuffle(mons)
    mons = mons[: min(n_terms, len(mons))]

    F = R.zero()
    for m in mons:
        c = random_int(rng, coeff_bound)
        if c != 0:
            F += QQ(c) * m

    if F == 0:
        F = QQ(1) * mons[0]

    return F


# ------------------------------------------------------------
# Experiment loop
# ------------------------------------------------------------

def run_experiment(
    *,
    n_samples,
    prime,
    out_path=None,
    stats_path=None,
    checkpoint_every=25,
    seed=0,
    n_terms=8,
    coeff_bound=10,
    max_tries_factor=50,
    verbose=True,
    store_only_ok=True,
):
    r"""
    Generate random smooth quartics over QQ, compute stable reduction, and optionally
    store successful examples (status == "ok") in a JSONL file.

    Convenience features:
      - periodic checkpoint prints every ``checkpoint_every`` samples
      - optional JSON stats file written at the end (and updated on checkpoints)

    Returns:
      (buckets, stats) where buckets maps type/status to list of results (in memory)
      and stats is a small dict with counts/timings.
    """
    v_K = QQ.valuation(prime)
    rng = Random(seed)

    buckets = defaultdict(list)

    tries = 0
    found = 0
    stored = 0
    max_tries = max_tries_factor * n_samples

    t_total = 0.0
    t_start = time.perf_counter()

    def compute_counts():
        """
        Return (status_counts, type_counts, bucket_sizes).
        - status_counts: counts for statuses ok/hyperelliptic/fail
        - type_counts: counts for reduction types among ok cases
        - bucket_sizes: legacy mixed dict (types + statuses) for backwards compatibility
        """
        bucket_sizes = {k: len(v) for k, v in buckets.items()}

        status_counts = {
            "ok": 0,
            "hyperelliptic": bucket_sizes.get("hyperelliptic", 0),
            "fail": bucket_sizes.get("fail", 0),
        }

        type_counts = {}
        for k, n in bucket_sizes.items():
            if k in ("hyperelliptic", "fail"):
                continue
            # everything else are reduction types (since key=SR.reduction_type for ok)
            type_counts[k] = n
            status_counts["ok"] += n

        return status_counts, type_counts, bucket_sizes

    def checkpoint(final=False):
        elapsed = time.perf_counter() - t_start
        status_counts, type_counts, bucket_sizes = compute_counts()
        stats = {
            "prime": prime,
            "seed": seed,
            "n_samples_target": n_samples,
            "found": found,
            "tries": tries,
            "stored": stored,
            "time_total_sec": t_total,
            "time_elapsed_sec": elapsed,
            "n_terms": n_terms,
            "coeff_bound": coeff_bound,
            "max_tries_factor": max_tries_factor,
            "out_path": out_path,
            "final": bool(final),
            "status_counts": status_counts,
            "type_counts": type_counts,
            "bucket_sizes": bucket_sizes,  # keep for backward compatibility
        }
        if verbose:
            msg = f"[checkpoint] found={found}/{n_samples}, tries={tries}, stored={stored}, elapsed={elapsed:.1f}s"
            print(msg)
        if stats_path is not None:
            write_run_stats(stats_path, stats)

    while found < n_samples and tries < max_tries:
        tries += 1

        F = random_quartic_form_QQ(rng=rng, n_terms=n_terms, coeff_bound=coeff_bound)

        # accept only smooth quartics over QQ
        try:
            if not Curve(F).is_smooth():
                continue
        except Exception:
            continue

        found += 1
        if verbose:
            print(f"[{found}/{n_samples}] smooth quartic found")
            print(f"  F = {F}")

        SR, t = timed_call(stable_reduction_of_quartic, F, v_K)
        t_total += t

        if verbose:
            print(f"  time = {t:.3f} sec")
            print(f"  -> status={SR.status}, type={getattr(SR, 'reduction_type', None)}")

        # bucket key for in-memory summary
        key = SR.reduction_type if SR.status == "ok" else SR.status
        buckets[key].append(SR)

        # write to database
        if out_path is not None:
            if (not store_only_ok) or SR.status == "ok":
                rec = SR.to_json_record()
                if rec is not None:
                    rec["time_sec"] = t
                    append_jsonl(out_path, rec)
                    stored += 1

        # periodic checkpoint
        if checkpoint_every and (found % checkpoint_every == 0):
            checkpoint(final=False)

    # final checkpoint + summary
    checkpoint(final=True)

    if verbose:
        print("")
        print(f"Generated {found} smooth quartics in {tries} attempts.")
        print(f"Total time (stable reduction calls): {t_total:.1f} sec")
        if out_path is not None:
            print(f"Stored {stored} records in: {out_path}")
        if stats_path is not None:
            print(f"Wrote run stats to: {stats_path}")
        print("Bucket sizes:")
        for k in sorted(buckets.keys()):
            print(f"  {k}: {len(buckets[k])}")

    status_counts, type_counts, bucket_sizes = compute_counts()
    stats = {
        "prime": prime,
        "seed": seed,
        "n_samples_target": n_samples,
        "found": found,
        "tries": tries,
        "stored": stored,
        "time_total_sec": t_total,
        "time_elapsed_sec": time.perf_counter() - t_start,
        "status_counts": status_counts,
        "type_counts": type_counts,
        "bucket_sizes": bucket_sizes,
        "out_path": out_path,
        "stats_path": stats_path,
    }
    return buckets, stats


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Sample random smooth plane quartics over QQ and compute p-adic stable reduction types."
    )
    parser.add_argument("--n-samples", type=int, default=50,
                        help="Number of smooth quartics to process (default: 50)")
    parser.add_argument("--prime", type=int, default=2,
                        help="Prime p for the p-adic valuation on QQ (default: 2)")
    parser.add_argument("--seed", type=int, default=0,
                        help="RNG seed (default: 0)")
    parser.add_argument("--n-terms", type=int, default=8,
                        help="Number of monomials used in random quartic (default: 8)")
    parser.add_argument("--coeff-bound", type=int, default=10,
                        help="Coefficient bound for random quartic (default: 10)")
    parser.add_argument("--max-tries-factor", type=int, default=50,
                        help="Max attempts = factor * n_samples (default: 50)")
    parser.add_argument("--out", type=str, default=None,
                        help="Path to JSONL output file (append-only). If omitted, no file is written.")
    parser.add_argument("--stats", type=str, default=None,
                        help="Path to JSON stats file (overwritten). If omitted, no stats file is written.")
    parser.add_argument("--checkpoint-every", type=int, default=25,
                        help="Checkpoint interval in samples (default: 25; 0 disables).")
    parser.add_argument("--quiet", action="store_true",
                        help="Reduce output")
    parser.add_argument("--store-all", action="store_true",
                        help='Store also "hyperelliptic"/"fail" (not recommended).')

    args = parser.parse_args()

    run_experiment(
        n_samples=args.n_samples,
        prime=args.prime,
        out_path=args.out,
        stats_path=args.stats,
        checkpoint_every=args.checkpoint_every,
        seed=args.seed,
        n_terms=args.n_terms,
        coeff_bound=args.coeff_bound,
        max_tries_factor=args.max_tries_factor,
        verbose=(not args.quiet),
        store_only_ok=(not args.store_all),
    )


if __name__ == "__main__":
    main()
