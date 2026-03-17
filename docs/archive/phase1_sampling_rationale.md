# Phase 1 PPI Sampling Rationale

## Purpose

This note documents the scientific rationale for using `20,000` global docking
trajectories per seed and `5` independent random seeds per receptor state for
Phase 1 PyRosetta docking.

Under this policy, the planned production budget is:

- `20,000` models per seed
- `5` seeds per receptor state
- `100,000` total models per receptor state
- `300,000` total models across the 3 receptor states

## Decision Summary

We recommend `20k x 5 seeds` rather than a single `100k` run for each receptor
state.

The total sampling budget is the same, but the multi-seed design is preferred
because it supports both:

1. sufficient global-search coverage for blind protein-protein docking, and
2. explicit convergence and robustness assessment across independent stochastic
   restarts.

In practical terms, `20k` models for one seed is acceptable as a first-pass
global docking budget, but it is on the lower side of the Rosetta-recommended
range for global docking. `100k` models per receptor state is a more defensible
production-level target for manuscript claims, especially when the docking is
fully blind and the interface is not prelocalized.

## Scientific Rationale

### 1. Global docking in Rosetta typically requires large decoy counts

Rosetta's protein-protein docking tutorial states that global docking requires
many runs to converge because of the large search space, typically
`10,000-100,000` trajectories. The RosettaDock application documentation gives
the same recommendation and further notes that for global docking one should
generate at least `10,000` decoys and ideally `100,000`.

This places `20,000` models per seed inside the accepted operating range, but
closer to the low-to-middle part of that range than to the high-confidence
production end.

### 2. A single 20k run is suitable for screening, but weak for strong negative claims

For blind docking, a `20k` run can be sufficient to detect strong, repeatedly
sampled interface basins if they are accessible under the protocol.
However, failure to recover a stable interface in a single `20k` run should not
be over-interpreted as evidence that no plausible interface exists. In other
words, `20k` is useful for positive signal discovery, but not ideal for
excluding alternatives.

### 3. Five independent 20k seeds are preferable to one 100k run

The preference for `5 x 20k` over `1 x 100k` is an inference from stochastic
search practice and from the scientific goal of measuring reproducibility, not
from a Rosetta rule that explicitly mandates multiple seeds.

The logic is:

- both designs spend the same total budget per receptor state (`100,000`
  models),
- but the multi-seed design provides five independent stochastic restarts,
- so it allows us to ask whether an interface patch is reproduced across
  independent runs rather than produced only once in one long trajectory set.

For manuscript interpretation, repeated recovery of the same receptor-side patch
across independent seeds is stronger evidence than a comparable number of
decoys generated from a single seed alone. This is especially useful when the
pipeline is later summarized at the residue or patch level.

### 4. The system size is still compatible with global docking

The current Phase 1 setup uses a receptor of about `309` residues and an
extended beta-meander partner of about `52` residues, for a combined docking
system of roughly `361` residues. Rosetta's tutorial notes that global docking
is most appropriate for relatively small complexes; this system remains within a
practical range for such a strategy.

### 5. Sampling budget does not remove the need to think about conformational flexibility

Rosetta global docking still assumes a fixed-backbone representation during the
core search. The tutorial explicitly notes that if backbone differences between
unbound and bound states are substantial, ensemble docking or related
flexibility-aware strategies become important.

Therefore, increasing from `20k` to `100k` per receptor state improves sampling
coverage and reproducibility, but it does not by itself solve errors caused by
backbone mismatch or partner flexibility. If seed-to-seed agreement remains poor
even at `100k` models per state, the limiting factor may be conformational
modeling rather than raw decoy count.

## Practical Interpretation For This Project

For this project, the following interpretation policy is reasonable:

- `20k x 1 seed`:
  suitable for smoke tests, pipeline validation, and first-pass interface
  discovery.
- `20k x 3 seeds` (`60k/state`):
  suitable as an intermediate production checkpoint when compute time is limited.
- `20k x 5 seeds` (`100k/state`):
  recommended for final production and manuscript-supported interpretation.

This makes `20k x 5 seeds` the best balance between:

- staying inside Rosetta's recommended global docking range,
- preserving independent-run reproducibility information, and
- avoiding over-reliance on one stochastic trajectory set.

## Manuscript-Ready Wording

### Methods paragraph

> For Phase 1 protein-protein docking, we performed blind global docking with
> PyRosetta/RosettaDock using 20,000 trajectories per random seed and five
> independent seeds for each receptor state, yielding 100,000 trajectories per
> state. This sampling level was chosen because Rosetta global docking typically
> requires on the order of 10,000-100,000 decoys for convergence, with 100,000
> decoys representing a commonly recommended upper production target for global
> searches. We distributed this budget across independent seeds rather than a
> single run in order to evaluate reproducibility of recovered interface patches
> across stochastic restarts.

### Results / interpretation paragraph

> Interface patches that were recovered across multiple independent seeds were
> interpreted as more robust than patches observed in only a single seed, even
> when the total number of sampled decoys was comparable. This multi-seed design
> allowed us to distinguish reproducible docking features from seed-specific
> stochastic outcomes.

### Limitations paragraph

> Increasing the total number of rigid-body docking trajectories improves search
> coverage but does not fully address limitations associated with backbone
> mismatch or partner flexibility. Accordingly, lack of convergence at a given
> site was interpreted cautiously, particularly for a flexible beta-meander-like
> partner, and was not treated as definitive evidence of absence of binding.

## Recommended Reviewer-Facing Claim

The safest strong claim is not "20,000 trajectories were sufficient", but
rather:

> We used a total sampling budget of 100,000 trajectories per receptor state,
> distributed across five independent seeds, to balance broad global-search
> coverage with explicit reproducibility assessment across stochastic restarts.

That formulation is both scientifically accurate and easier to defend than a
claim based only on a single-seed trajectory count.

## References

Primary sources used for this rationale:

1. Rosetta protein-protein docking tutorial:
   https://docs.rosettacommons.org/demos/latest/tutorials/Protein-Protein-Docking/Protein-Protein-Docking
2. RosettaDock application documentation:
   https://docs.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol
3. Marze NA, Roy Burman SS, Sheffler W, Gray JJ. Efficient flexible backbone
   protein-protein docking for challenging targets. Bioinformatics. 2018.
   RosettaDock 4.0 benchmark paper:
   https://pmc.ncbi.nlm.nih.gov/articles/PMC6184633/
4. Alam N, Goldstein O, Xia B, Porter KA, Kozakov D, Schueler-Furman O.
   High-resolution global peptide-protein docking using fragments-based
   PIPER-FlexPepDock. PLoS Comput Biol. 2017.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC5760072/
