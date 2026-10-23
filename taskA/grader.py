"""Deterministic grader for mgx_drug_coprescription_confounding."""
import json
import sys


def norm(x):
    return str(x).strip().lower().replace(" ", "_").replace("-", "_")


def grade(key, sub):
    r = {}
    block = sub.get("drug_labels") or {}
    for drug, rule in key["drug_labels"].items():
        r[f"drug_labels.{drug}"] = norm(block.get(drug)) in [norm(a) for a in rule["accept"]]
    for field in ["strongest_independent_drug", "drug_explaining_sulfonylurea_signal"]:
        r[field] = norm(sub.get(field)) in [norm(a) for a in key[field]["accept"]]
    try:
        v = float(sub.get("statin_metformin_overlap"))
        k = key["statin_metformin_overlap"]
        r["statin_metformin_overlap"] = abs(v - k["value"]) <= k["tolerance"]
    except (TypeError, ValueError):
        r["statin_metformin_overlap"] = False
    return r


if __name__ == "__main__":
    key = json.load(open(sys.argv[1]))
    sub = json.load(open(sys.argv[2]))
    res = grade(key, sub)
    for k, ok in res.items():
        print(f"{'ok  ' if ok else 'FAIL'} {k}")
    print("PASS" if all(res.values()) else "FAIL", f"({sum(res.values())}/{len(res)} fields)")
    sys.exit(0 if all(res.values()) else 1)
