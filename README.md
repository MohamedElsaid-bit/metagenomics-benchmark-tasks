# Latch task authoring exercise

Mohamed Elsaid

Start with **[SUBMISSION.md](SUBMISSION.md)**. It covers both tasks, what the testing showed, and the defects testing found.

## The two submitted tasks

| Folder | Task | Spec |
|---|---|---|
| `taskA/` | `mgx_drug_coprescription_confounding` | [taskA/TASK.md](taskA/TASK.md) |
| `taskB/` | `mgx_crc_partial_contamination` | [taskB/TASK.md](taskB/TASK.md) |

Each task folder contains:

```
TASK.md          prompt, dataset provenance, expected answer, robustness, skill, observed LLM behavior
data/            the exact files handed to the agent, and nothing else
answer_key.json  machine readable expected answer
grader.py        deterministic grader
llm_runs/        every test run: the agent's JSON answer and its own analysis notes
```

To grade a submission:

```
python grader.py answer_key.json <submission.json>
```

It prints one line per field and exits non zero on failure. Example, re-grading a stored run:

```
cd taskB && python grader.py answer_key.json llm_runs/r3_sonnet/answer.json
```

Requires Python 3 with pandas, scipy and statsmodels for the build and analysis scripts; the graders themselves need only the standard library.

## Supporting material

| Folder | Contents |
|---|---|
| `scripts/` | every build and analysis script, so all numbers in the task files are reproducible |
| `earlier-versions/taskB-v1-whole-sample-swap/` | first version of Task B, using whole mislabelled plaque samples instead of partial mixtures, 6 test runs, kept for the difficulty comparison |
| `earlier-versions/taskA-v1-t2d-vs-healthy/` | first version of Task A, framed as T2D versus healthy, 3 test runs |
| `results/` | robustness output tables |

Files ending `_PRIVATE` (`taskB/truth_PRIVATE.json`, the `id_map_PRIVATE.tsv` files) hold the ground truth and the mapping back to original sample IDs. They are author material for verification and are never given to the agent.

## Data

All data comes from curatedMetagenomicData (Bioconductor), which redistributes uniformly processed MetaPhlAn3 profiles under open access. No controlled access resource is used. Source studies and citations are listed in each task spec and in SUBMISSION.md.

## A note on publishing this openly

This repository is public so it can be reviewed without any access friction, which means the answer keys and the ground truth for the contaminated samples are public too. Both tasks should therefore be treated as worked examples rather than as live benchmark items.

If either task is wanted for actual use, it can be regenerated as an unseen variant. `scripts/t3_build_and_key.py` takes a fixed seed (`SEED`), the number of contaminated samples in each arm and the contamination fraction range; changing the seed produces a different set of affected samples, different anonymized IDs and a matching answer key, emitted in the same run so the two cannot drift apart. Task A is built from an unmodified public cohort, so a fresh variant there means changing the graded drug set or the analysis population rather than regenerating data.
