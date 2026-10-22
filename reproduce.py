"""Re-run every exact check of the article and compare with the stored outputs.

Usage:  python3 -B reproduce.py [--output FILE]

Everything runs in a temporary copy of the packet; nothing in the packet is
modified. Steps (about 30 seconds in total on a laptop):

  1. exact/check_exact.py: elimination of S from the vector field by a resultant,
     k = 1..12, in integer polynomial arithmetic (Python standard library only).
     Every per-k line (including a SHA256 of each eliminated polynomial) must equal
     results/exact-checks.jsonl; the final line must report PASS.
  2. independent/check_proof_algebra.py (SymPy): the separately written check of
     the formulas, k = 2..9, symbolic k, irreducibility, the identities of Lemma 5,
     and the single trajectory for k = 2, 3, 5. Output must equal
     results/check_proof_algebra.out and end with ALL PASS.
  3. independent/check_article_algebra.py (SymPy): a second separately written
     recomputation from the article text, k = 1..12. Output must equal
     results/check_article_algebra.out and end with ALL OK.
  4. independent/check_supplement_s9.py (SymPy): the input-output equations printed
     in Supplement S9 for k = 2, 3, 4 plus P must be zero. Output must equal
     results/check_supplement_s9.out.

All comparisons are exact (text equality); the computations are exact symbolic or
rational arithmetic. They are finite checks of the article's formulas, not a proof
of the general theorem.
"""
import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = [sys.executable, "-B"]


def run(script, cwd):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1")
    t = time.time()
    out = subprocess.run(PY + [script], cwd=cwd, env=env, capture_output=True, text=True)
    assert out.returncode == 0, (script, out.stderr[-3000:])
    return out.stdout, round(time.time() - t, 2)


def same_text(new, stored_path, name):
    stored = Path(stored_path).read_text()
    assert new == stored, (name, "output differs from the stored file", stored_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output")
    a = ap.parse_args()
    t0 = time.time()
    report = dict(status="RUNNING", steps={})
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "packet"
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git", "__pycache__"))

        # 1. exact elimination, k = 1..12
        out, sec = run("exact/check_exact.py", work)
        new = out.splitlines()
        stored = (ROOT / "results/exact-checks.jsonl").read_text().splitlines()
        assert len(new) == len(stored) == 13, len(new)
        assert new[:12] == stored[:12], "per-k lines differ"
        last_new, last_stored = json.loads(new[12]), json.loads(stored[12])
        assert last_new["status"] == last_stored["status"] == "PASS"
        # The stored line carries the SHA256 of the research original of the script; the public
        # copy differs only in its docstring (see source-derivatives.json).
        derivs = json.loads((ROOT / "source-derivatives.json").read_text())["records"]
        orig = {r["public"]["path"]: r["original_sha256"] for r in derivs}
        public_sha = hashlib.sha256((ROOT / "exact/check_exact.py").read_bytes()).hexdigest()
        assert last_new["own_code_sha256"] == public_sha
        assert last_stored["own_code_sha256"] == orig["exact/check_exact.py"]
        assert {k: v for k, v in last_new.items() if k != "own_code_sha256"} == \
               {k: v for k, v in last_stored.items() if k != "own_code_sha256"}
        rows = [json.loads(x) for x in new[:12]]
        assert all(r["resultant_equals_closed"] and r["vector_field_identity"] and r["primitive_checks"]
                   and r["recovery_coefficient_identities"] for r in rows)
        assert [r["same_alt_parameters"] for r in rows] == [True] + [False] * 11
        assert [r["source_s9_equal_up_to_sign"] for r in rows][1:4] == [True, True, True]
        report["steps"]["1-exact-elimination-k1-12"] = dict(status="PASS", seconds=sec)

        # 2-4. separately written SymPy checks
        for key, script, stored_out, final in (
                ("2-proof-algebra", "independent/check_proof_algebra.py", "results/check_proof_algebra.out", "ALL PASS"),
                ("3-article-algebra", "independent/check_article_algebra.py", "results/check_article_algebra.out", "ALL OK"),
                ("4-supplement-s9", "independent/check_supplement_s9.py", "results/check_supplement_s9.out", None)):
            out, sec = run(script, work)
            same_text(out, ROOT / stored_out, script)
            lines = out.strip().splitlines()
            if final:
                assert lines[-1] == final and not any(x.startswith("FAIL") for x in lines), script
            else:
                assert len(lines) == 3 and all(x.endswith("= 0") for x in lines), script
            report["steps"][key] = dict(status="PASS", seconds=sec)

    import sympy
    report.update(status="PASS", seconds=round(time.time() - t0, 1), python=platform.python_version(),
                  sympy=sympy.__version__)
    if a.output:
        Path(a.output).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(status="PASS", seconds=report["seconds"], steps=list(report["steps"]))))


if __name__ == "__main__":
    main()
