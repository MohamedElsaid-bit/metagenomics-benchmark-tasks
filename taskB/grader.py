"""Deterministic grader for mgx_crc_partial_contamination.

Usage: python grader.py answer_key.json submission.json
Pass requires every field to be correct.

Flag-set rule: the submission must flag every contaminated sample (no misses) and may
flag at most `max_extra` additional ones. Missing a contaminated sample biases the case
control comparison; dropping a few clean samples only costs power, and was verified not
to change any graded species call (see scripts/t3_flagset_sensitivity.py).
"""
import json
import sys


def norm(x):
    return str(x).strip().lower().replace(" ", "_").replace("-", "_")


def grade(key, sub):
    r = {}
    given = sub.get("flagged_sample_ids")
    ids = {str(x).strip().upper() for x in given} if isinstance(given, list) else set()
    required = set(key["flagged_sample_ids"]["required"])
    max_extra = key["flagged_sample_ids"]["max_extra"]
    missed, extra = required - ids, ids - required
    r["flagged_sample_ids.no_missed"] = not missed
    r[f"flagged_sample_ids.extra_at_most_{max_extra}"] = len(extra) <= max_extra
    try:
        r["n_flagged_samples_consistent"] = int(sub.get("n_flagged_samples")) == len(ids)
    except (TypeError, ValueError):
        r["n_flagged_samples_consistent"] = False
    r["flagged_sample_origin"] = norm(sub.get("flagged_sample_origin")) in [norm(a) for a in key["flagged_sample_origin"]["accept"]]
    block = sub.get("crc_vs_control") or {}
    for species, rule in key["crc_vs_control"].items():
        r[f"crc_vs_control.{species}"] = norm(block.get(species)) in [norm(a) for a in rule["accept"]]
    return r, missed, extra


if __name__ == "__main__":
    key = json.load(open(sys.argv[1]))
    sub = json.load(open(sys.argv[2]))
    res, missed, extra = grade(key, sub)
    for k, ok in res.items():
        print(f"{'ok  ' if ok else 'FAIL'} {k}")
    if missed:
        print(f"     missed: {sorted(missed)}")
    if extra:
        print(f"     extra:  {sorted(extra)}")
    print("PASS" if all(res.values()) else "FAIL", f"({sum(res.values())}/{len(res)} fields)")
    sys.exit(0 if all(res.values()) else 1)
