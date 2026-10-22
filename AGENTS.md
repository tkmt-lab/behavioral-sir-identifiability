# Entry point for a computational reader

Start with README.md, claims.json and VERIFICATION-RECORD.md. The main result is Theorem 6 of the article, paper-behavioral-sir-v1.1-en.tex.

**What the theorem says.** For every fixed, known integer k ≥ 2, the behavioral SIR model (1), observed through I(t) with known total population N0, is globally structurally identifiable. This holds for generic parameters and initial values, in the sense of Definition 1, and it carries over to the reals as an open dense set.

**What it rests on.** The proof is in prose and is not machine-checked. It depends on published results: Ovchinnikov–Pogudin–Thompson, Theorem 4.7, Corollary 5.7 and Proposition 3.1, numbered as in arXiv:2007.14787v2; and Hong–Ovchinnikov–Pogudin–Yap, Lemma 3.2, numbered as in arXiv:1801.08112v5.

**How not to describe it.** Do not:
- extend the result to non-integer k or to k treated as unknown;
- extend it to every (non-generic) parameter value or to estimation from noisy or discrete data;
- describe it as formally verified or peer reviewed;
- present the finite computations as a proof for all k.

**Running the checks.** Neither command needs network access.
- `python3 -B check_packet.py` checks the file inventory, the claim scopes and privacy.
- `python3 -B reproduce.py` takes about 30 seconds. It re-runs:
  - the exact elimination for k = 1..12 (`exact/check_exact.py`, standard library only);
  - two separately written SymPy checks (`independent/check_proof_algebra.py`, `independent/check_article_algebra.py`);
  - the comparison with the printed Supplement S9 equations (`independent/check_supplement_s9.py`).

  Each output is compared exactly with the stored file in `results/`.

**Evidence types.** claims.json keeps them apart:
- `prose-proof`: the article's proof.
- `ai-review`: two AI reader sessions, one on the proof notes and one on the article; not human peer review.
- `finite-exact-computation`: exact symbolic or rational checks for finitely many k or one trajectory.

**Equation numbers.** `independent/check_proof_algebra.py` labels equations by the numbering of the proof notes, which are not distributed. Its docstring maps them to the article's equations. `independent/check_article_algebra.py` uses the article's numbering.

**The k = 1 case.** The coefficients are invariant under (β, γ) → (αγ, β/α), which agrees with Theorem 3.8 of Pant, Saucedo and Pogudin. The recovery formulas (7) do not apply to k = 1.
