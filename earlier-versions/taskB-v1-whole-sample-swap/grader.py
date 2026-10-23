"""Deterministic grader for mgx_crc_sample_provenance.

Usage: python grader.py answer_key.json submission.json
Pass requires every field to be correct.
"""
import json
import sys


def norm(x):
    return str(x).strip().lower().replace(" ", "_").replace("-", "_")


def grade(key, sub):
    r = {}
    try:
        r["n_flagged_samples"] = int(sub.get("n_flagged_samples")) == key["n_flagged_samples"]["value"]
    except (TypeError, ValueError):
        r["n_flagged_samples"] = False
    given = sub.get("flagged_sample_ids")
    ids = {str(x).strip().upper() for x in given} if isinstance(given, list) else set()
    r["flagged_sample_ids"] = ids == set(key["flagged_sample_ids"]["accept"])
    r["flagged_sample_origin"] = norm(sub.get("flagged_sample_origin")) in [norm(a) for a in key["flagged_sample_origin"]["accept"]]
    block = sub.get("crc_vs_control") or {}
    for species, rule in key["crc_vs_control"].items():
        r[f"crc_vs_control.{species}"] = norm(block.get(species)) in [norm(a) for a in rule["accept"]]
    return r


if __name__ == "__main__":
    key = json.load(open(sys.argv[1]))
    sub = json.load(open(sys.argv[2]))
    res = grade(key, sub)
    for k, ok in res.items():
        print(f"{'ok  ' if ok else 'FAIL'} {k}")
    print("PASS" if all(res.values()) else "FAIL", f"({sum(res.values())}/{len(res)} fields)")
    sys.exit(0 if all(res.values()) else 1)
