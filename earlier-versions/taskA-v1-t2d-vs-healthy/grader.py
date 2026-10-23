"""Deterministic grader for mgx_t2d_drug_confounding.

Usage: python grader.py answer_key.json submission.json
Pass requires every field to be correct. Prints per-field results and PASS/FAIL.
"""
import json
import sys


def norm(x):
    return str(x).strip().lower().replace(" ", "_").replace("-", "_")


def grade(key, sub):
    results = {}
    results["primary_driver"] = norm(sub.get("primary_driver")) in [norm(a) for a in key["primary_driver"]["accept"]]
    try:
        v = float(sub.get("driver_fraction_T2D"))
        k = key["driver_fraction_T2D"]
        results["driver_fraction_T2D"] = abs(v - k["value"]) <= k["tolerance"]
    except (TypeError, ValueError):
        results["driver_fraction_T2D"] = False
    for block in ["T2D_intrinsic", "driver_effect_within_T2D"]:
        given = sub.get(block) or {}
        for species, rule in key[block].items():
            results[f"{block}.{species}"] = norm(given.get(species)) in [norm(a) for a in rule["accept"]]
    return results


if __name__ == "__main__":
    key = json.load(open(sys.argv[1]))
    sub = json.load(open(sys.argv[2]))
    res = grade(key, sub)
    for k, ok in res.items():
        print(f"{'ok  ' if ok else 'FAIL'} {k}")
    passed = all(res.values())
    print("PASS" if passed else "FAIL", f"({sum(res.values())}/{len(res)} fields)")
    sys.exit(0 if passed else 1)
