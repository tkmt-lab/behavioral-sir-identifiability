# Global identifiability of a behavioral SIR model for integer k ≥ 2

Release 1.0.0 of https://github.com/tkmt-lab/behavioral-sir-identifiability (tag `v1.0.0`). The article is [paper-behavioral-sir-v1.1-en.pdf](paper-behavioral-sir-v1.1-en.pdf) (version 1.1, 7 pages).

## In plain words

In this epidemic model, people reduce their contacts when many others are infected. The question is whether a perfect record of the number of infected people over time, I(t), together with the total population, is enough to determine the three model parameters and the hidden numbers of susceptible and recovered people. Pant, Saucedo and Pogudin showed that the answer is no when the saturation exponent is k = 1. Software checks for k = 2, 3, 4, 5, 10 and 50 suggested that the answer is yes, and they stated a conjecture for k ≥ 2. The article proves that the answer is yes for every whole number k ≥ 2.

## Exact claim

Model (k a fixed, known integer ≥ 2):

    S' = −β S I / (1 + α I)^k,   I' = β S I / (1 + α I)^k − γ I,   R' = γ I,

with output y = I and known N0 = S(0) + I(0) + R(0). **Theorem 6:** for generic parameters and initial values, α, β, γ and S(0), I(0), R(0) are globally structurally identifiable. Explicitly, the minimal input–output equation normalized so that the coefficient of y·y'' is 1 has

    [y² y''] = kα,   [y² y'] = β + kαγ,   [y³ y'] = k(k−1)α²γ,

so α, γ and β follow by rational formulas whose denominators are k and k(k−1)α². Then S = (y' + γy)(1 + αy)^k / (βy) and R = N0 − S − I.

## Scope

- **Proved in prose.** The proof has two parts. The first is a closed form of the input–output equation and a check that it is a monic characteristic presentation. The second is a lemma (Lemma 5) that the two-state model has no nonconstant rational first integral. It uses published theorems of Ovchinnikov, Pogudin and Thompson (Theorem 4.7, Corollary 5.7, Proposition 3.1) and Lemma 3.2 of Hong, Ovchinnikov, Pogudin and Yap.
- **Checked by computation, for finitely many cases.** The checks below cover k = 1..12, k = 2..9 and one trajectory. They confirm the formulas but do not replace the proof.
- **Read by AI.** Two AI reader sessions checked the proof and the article, each recomputing parts of it. This is not human peer review.
- **Not claimed:**
  - a Lean or other machine-checked proof;
  - human expert review;
  - non-integer k, or k unknown;
  - a description of the exceptional (non-generic) set;
  - estimation from discrete or noisy data;
  - priority beyond a limited literature search.

## Check in a minute

You need Python 3.9 or later and the packages in requirements.txt (`python3 -m pip install -r requirements.txt`). From the repository root, run:

```sh
python3 -B reproduce.py
```

It works in a temporary copy and takes about 30 seconds. It re-runs four checks, and each output must match the stored file exactly:

1. The exact elimination for k = 1..12, which uses only the standard library.
2. A SymPy check of the formulas, of Lemma 5's identities and of a single trajectory.
3. A second SymPy recomputation from the article text.
4. The comparison with the equations printed in Supplement S9 of the original preprint.

The last line of the output is

```text
{"status": "PASS", "seconds": ..., "steps": ["1-exact-elimination-k1-12", "2-proof-algebra", "3-article-algebra", "4-supplement-s9"]}
```

To see the recovery formulas at work for k = 3:

```sh
python3 -B -c "
import sympy as sp
a,b,g,y,v,w = sp.symbols('alpha beta gamma y v w'); k = 3; u = 1 + a*y
P = sp.Poly(sp.expand(y*u**k*w - (1-(k-1)*a*y)*u**(k-1)*v**2 + (k*a*g*y**2*u**(k-1) + b*y**2)*v + b*g*y**3), y, v, w)
c = dict(P.terms()); al = c[(2,0,1)]/k; ga = c[(3,1,0)]/(k*(k-1)*al**2); be = c[(2,1,0)] - k*al*ga
print(sp.simplify(al), sp.simplify(be), sp.simplify(ga))
"
```

Expected output: `alpha beta gamma`.

## Check the evidence

- Read [AGENTS.md](AGENTS.md), [claims.json](claims.json) and [VERIFICATION-RECORD.md](VERIFICATION-RECORD.md).
- `python3 -B check_packet.py` checks the file list and hashes, the claim scopes, the article's references and the absence of private paths.
- [exact/](exact) holds the elimination code. [independent/](independent) holds the two separately written SymPy checks and the Supplement S9 comparison. [results/](results) holds their stored outputs.

## Research and verification

The author is Tatsuhiko Takimoto, Independent researcher. AI systems assisted with the derivation and checking of the proof, the exact computer-algebra checks, independent readings of the proof and of the manuscript, literature organization and drafting. The author is responsible for the content and its stated scope. Process details are in [VERIFICATION-RECORD.md](VERIFICATION-RECORD.md).

## Built on prior work

- Binod Pant, Omar Saucedo and Gleb Pogudin posed the question and the conjecture. They proved the k = 1 result, ran the software computations for individual k and sketched the large-k argument. Their code is public. See [medRxiv 10.64898/2026.01.16.26344284](https://www.medrxiv.org/content/10.64898/2026.01.16.26344284v2) and its supplement.
- Alexey Ovchinnikov, Gleb Pogudin and Peter Thompson related input–output equations, characteristic presentations and identifiability ([arXiv:2007.14787](https://arxiv.org/abs/2007.14787v2); AAECC 34 (2023)). These results carry the proof.
- Hoon Hong, Alexey Ovchinnikov, Gleb Pogudin and Chee Yap laid the differential-algebra foundation of global identifiability used here ([arXiv:1801.08112](https://arxiv.org/abs/1801.08112v5); CPAM 73 (2020)).

I am grateful to them, and I hope this work adds something new to the field they opened. These acknowledgements imply no endorsement.

Code: Apache-2.0. Article, original text and generated data: CC BY 4.0. Contact: GitHub Issues at the repository.
