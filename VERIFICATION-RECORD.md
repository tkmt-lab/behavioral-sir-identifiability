# Verification record: behavioral-sir-identifiability 1.0.0

This record carries the process and provenance details that the article (version 1.1) leaves out. The article states the mathematics and its present scope. This file states who or what checked which part, with which tools, and when. All dates are 22 October 2026, times JST (UTC+9). The packet was assembled and checked as a local release candidate; publication is a separate step taken by the author.

## Status of the main result

Theorem 6 covers every fixed, known integer k ≥ 2 and is proved in prose. The proof uses three results of Ovchinnikov, Pogudin and Thompson, numbered as in arXiv:2007.14787v2:
- Theorem 4.7: identifiability equals input–output identifiability when there are no rational first integrals.
- Corollary 5.7: the coefficients of a monic characteristic presentation generate the input–output identifiable functions.
- Proposition 3.1: the field criterion for identifiability.

It also uses Lemma 3.2 of Hong, Ovchinnikov, Pogudin and Yap, numbered as in arXiv:1801.08112v5. The hypotheses of these results are verified in the article (Proposition 3 and Lemma 5).

There is no Lean or other machine-checked formalization and no human expert review. The finite computations below check formulas; they do not replace the proof.

## Computations and their independence

1. **Exact elimination, k = 1..12** (`exact/check_exact.py`; stored output `results/exact-checks.jsonl`).
   - Tools: Python 3.14.7, standard library only (sparse integer polynomials, `fractions.Fraction`), macOS arm64.
   - Run on 22 October 2026, 10:26 JST, in about one second under a memory cap, single-threaded. Re-run in a Python 3.11 environment for this packet with identical output.
   - What it does: eliminates S by a 2×2 resultant and compares every monomial with the closed coefficient formulas. It checks the vector-field identity as an exact polynomial identity, the recovery coefficients, the values used for primitivity, the scaled derivation, and the k = 1 two-point control. It also compares a separate transcription of the Supplement S9 equations for k = 2, 3, 4.
   - Independence: written by the same investigator who derived the formulas, after reading an earlier exploratory script. It is a re-implementation, not a blinded derivation.
   - An earlier exploratory script (not distributed) was also re-run for k = 1..12 and agreed on all coefficients in a direct comparison.
2. **Proof check by a separate reader** (`independent/check_proof_algebra.py`; output `results/check_proof_algebra.out`).
   - Who and when: an AI reader session that did not take part in writing the proof, finished at 10:55 JST.
   - Inputs: it read the statement, the proof notes and the original sources, which it fetched itself and matched against recorded SHA256 values. It did not read the elimination code, its results or the study report.
   - Tools: SymPy 1.13.3; about 6 s.
   - What it checks:
     - the input–output equation on the model for k = 2, 3, 5 and for symbolic k;
     - the coefficient table and its inverse for k = 2..9;
     - irreducibility for k = 2, 3, 5;
     - the identities of Lemma 5;
     - that no polynomial first integral of degree ≤ 6 exists for k = 2, 3;
     - one trajectory as an exact power series to order 40 for k = 2, 3, 5.
   - The stored output was regenerated for this packet with Python 3.11.14 and SymPy 1.13.3. The reader reported every item as PASS.
3. **Article check by a second separate reader** (`independent/check_article_algebra.py` and `independent/check_supplement_s9.py`; outputs in `results/`).
   - Who and when: a further AI reader session, acting as referee on the draft article (version 1), finished at 13:58 JST.
   - Inputs: it worked from the article text, fetched the sources itself, and read neither the elimination code nor its results.
   - Tools: the reader ran SymPy 1.14.0. The packet re-runs with SymPy 1.13.3 give byte-identical output.
   - What it checks:
     - the derivation for symbolic k;
     - the coefficient table for k = 1..12 and the inverse map;
     - irreducibility and the Lemma 2 values for k = 2..6;
     - the k = 1 invariance;
     - the scaled derivation;
     - one trajectory with its own series code (order 26, k = 2, 3).
   - Separately, the Supplement S9 equations for k = 2, 3, 4, transcribed by this reader, equal −P term by term.

## AI readings (not human peer review)

- **Proof notes.** The reader of item 2 returned "correct". It made five optional suggestions, all applied before the article was written; none changed the result:
  - state the reduction from three states to two;
  - name the criterion used for the states and the condition βI(0)(1+αI(0)) ≠ 0;
  - give a source for the identification of the input–output ideal;
  - correct one word in the coefficient description and note the zero coefficient;
  - mention the large-k sketch of Supplement S9.
- **Article, version 1.** The referee of item 3 returned "minor revision" with no mathematical errors. It checked every cited theorem, definition and equation number against the sources. Version 1.1 addresses all of its findings:
  - the abstract's description of the denominators;
  - wording around Eq. (4), Theorem 6 and Section 3;
  - a missing citation of Supplement S9;
  - one omitted step showing that the input–output ideal is the kernel, now with a citation of Hong–Ovchinnikov–Pogudin–Yap;
  - naming the tools of the checks;
  - the journal version of Ovchinnikov–Pogudin–Thompson;
  - credit given only through citations in the article;
  - this record.
- **Bibliography.** The journal data for Ovchinnikov–Pogudin–Thompson (AAECC 34(2):165–182, 2023, published online 2021) and for Hong–Ovchinnikov–Pogudin–Yap (CPAM 73(9):1831–1879, 2020) were checked at Crossref and arXiv by the referee. Theorem numbers in the article follow the arXiv versions. The journal numbering was not checked because it is behind a paywall.
- **Literature.** The novelty search was limited: the original preprint versions, the authors' code and keyword searches. It found no other proof for general k. It was not exhaustive.

## Article versions

- 1: first draft.
- 1.1 (this release): response to the referee (above). The mathematics, every equation and every number are unchanged from version 1.

## Present packet checks

- `check_packet.py` checks:
  - the closed file list and SHA256SUMS;
  - the claims.json evidence fingerprints and fixed scope;
  - the article's labels, references and citations;
  - the version and repository fields and the PDF validation record;
  - the stored outputs;
  - the absence of private paths and addresses (`publication_privacy.py`).
- `reproduce.py` re-runs all four computations and requires text-identical output. About 30 s in a Python 3.11 environment.
- Every page of the PDF was rendered at 80 dpi and inspected (`pdf-validation.json`). There are no overfull boxes and no undefined references.
- Before release the packet was also scanned with a non-public list of internal identifiers (not distributed).
- No Lean run and no new literature search were made for the packet. Bounded privacy scans are not a guarantee against every unknown identifier.

## Files and labels

- Public copies whose text differs from the research originals are listed with both hashes in `source-derivatives.json`. The differences are in docstrings and comments only: internal labels and file names removed, and an equation-number key added. The stored outputs are unchanged.
- The last line of `results/exact-checks.jsonl` records the SHA256 of the research original of `exact/check_exact.py`. `reproduce.py` checks that it matches the original hash recorded in `source-derivatives.json`.
- Not distributed: the proof notes, the AI reviews, the earlier exploratory scripts and research-process records.
