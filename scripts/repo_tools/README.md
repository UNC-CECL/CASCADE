# `repo_tools/` — tools that act on the repository itself

Not site content and not part of any analysis: these read the tree and report
on it. They are **run**, never imported.

| Tool | Does |
|---|---|
| `hat_layout_check.py` | Audits the tree against the seven rules in `ORGANIZATION.md`. |

```
python scripts/repo_tools/hat_layout_check.py            every rule
python scripts/repo_tools/hat_layout_check.py --rule 5   just one
python scripts/repo_tools/hat_layout_check.py --full     every offender
```

**It is advisory and always exits zero.** Chosen deliberately (Hannah,
2026-09-13): a check that blocks you mid-experiment gets disabled, and a
disabled check reports nothing. This one is meant to be run when you want a
picture, and it produces a worklist rather than an obstacle.

It finds the repo root by searching upward for `pyproject.toml`, so it does not
care where under the tree it is invoked from or where it is moved to.
