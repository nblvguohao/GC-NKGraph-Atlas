# GC-NKGraph-Atlas — Positioning & Manuscript Framing (revised)

> Replaces Section 0 (Project Identity) and Section 21 (Manuscript Framing).
> The change in this revision: the project is no longer framed as "a gastric
> cancer NK graph pipeline with a gated topology add-on." It is framed as a
> **computational operationalization of a specific, published immune-evasion
> mechanism**, validated by recovering that mechanism in its origin system
> (liver) and extended to a new, unclaimed system (gastric cancer).

---

## 1. Core positioning (one paragraph)

A landmark wet-lab finding (Zheng et al., Nat Immunol 2023, DOI
10.1038/s41590-023-01462-9) showed that tumors evade NK cytotoxicity by altering
the surface topology of NK cells: dysregulated tumor serine metabolism lowers
sphingomyelin in intratumoral NK membranes, which collapses membrane protrusions,
prevents lytic immune-synapse formation, and abolishes killing. That mechanism
was established with cutting-edge imaging and single-cell mass spectrometry on a
small number of liver-cancer patients. **This project asks whether that
mechanism can be reconstructed at scale from public transcriptomes, and whether
it generalizes to gastric cancer** — using a single-cell-informed, NK-aware,
heterogeneous graph-learning framework that (a) recovers the serine → SM →
protrusion-machinery → cytotoxicity axis in liver as a positive control,
(b) tests the same axis in gastric cancer as the novel result, and
(c) prioritizes tumor-intrinsic candidate targets in that axis with wet-lab-ready
validation assays.

The central deliverable is **not** a prognosis signature and **not** a public
atlas. It is a reproducible protocol that turns a mechanistic discovery into a
scalable target-discovery engine.

---

## 2. Two-arm design (this is the spine of the paper)

```text
Arm A — POSITIVE CONTROL (liver / HCC)
    The credibility anchor. Show the framework recovers the published
    serine->SM->topology->cytotoxicity axis in the system where it was proven.
    Cohorts: TCGA-LIHC (+ optional HCC scRNA with intratumoral-vs-peritumoral NK).
    Success = pre-registered axis hypotheses recovered in the expected direction.

Arm B — NOVEL EXTENSION (gastric cancer)
    The new contribution. Test whether the SAME axis operates in gastric cancer,
    a digestive-tract cancer NOT yet on the mechanism's published extension list
    (which already includes lung, colon, ovarian — see Clin Transl Med 2023,
    DOI 10.1002/ctm2.1395). Output gastric-specific tumor-intrinsic targets.
```

Why this structure matters: recovering a target lab's own published mechanism
from independent public data is the single most credible thing a computational
project can demonstrate to that lab. It also resolves the previous tension
between the "HCC-NK" heritage and the "gastric cancer" focus — liver is now the
principled control, gastric is the principled extension.

---

## 3. Recommended title

Primary:

**GC-NKGraph-Atlas: reconstructing the serine–sphingomyelin–membrane-topology
axis of NK-cell immune evasion from tumor transcriptomes, from liver to gastric
cancer**

Alternatives:

- **A single-cell-informed tumor–NK heterogeneous graph reconstructs a physical
  immune-evasion mechanism and prioritizes gastric-cancer targets**
- **Operationalizing NK surface-topology evasion: a knowledge-guided graph
  framework recovers a metabolic-topology axis in liver and extends it to gastric
  cancer**

---

## 4. Core claims (conservative, mechanism-grounded)

State only what is supported by generated evidence. Each claim maps to a
figure/table in `manuscript/notes/main_claims.md`.

1. We define a reproducible, cell-type-attributed transcriptional proxy for the
   serine → sphingomyelin → membrane-protrusion-machinery → cytotoxicity axis of
   NK immune evasion.
2. In liver cancer (the mechanism's origin system), the framework **recovers**
   the axis: pre-registered relationships between tumor serine capacity, NK SM
   balance, NK protrusion machinery, and cytotoxicity replicate in the expected
   direction from independent public cohorts.
3. We construct a single-cell-informed tumor–NK heterogeneous graph in which a
   mechanism-grounded `metabolic_crosstalk` edge (tumor serine program → NK
   topology state) is justified by the biology, not by generic priors — giving
   the heterogeneous-graph architecture a principled reason to exist.
4. The graph model improves or stabilizes NK-state prediction versus strong
   tabular and graph baselines under internal and external validation.
5. The same axis is tested in gastric cancer; where it holds, we prioritize
   tumor-intrinsic candidate targets supported by multi-evidence patterns
   (tumor specificity, axis membership, NK-dysfunction correlation, spatial-niche
   support, ligand–receptor context, graph attention, druggability), each with a
   recommended wet-lab validation assay.
6. A model-based in-silico "SM-restoration" readout stratifies samples by
   predicted benefit from an SM-restoration + Tim3-blockade combination logic —
   presented as a hypothesis for experimental testing, not a validated predictor.

---

## 5. What is novel vs what is a control (state this explicitly)

```text
Positive control (not novel, but essential):
  - recovery of the published serine-SM-topology-cytotoxicity axis in liver.

Novel:
  - scalable transcriptional reconstruction of that axis from bulk + scRNA;
  - the mechanism-grounded heterogeneous-graph formulation (metabolic_crosstalk edge);
  - the gastric-cancer extension and gastric-specific target list;
  - the mechanism-card abstraction that makes the whole pipeline transferable.
```

---

## 6. Avoid overclaiming (hard limits)

- Do **not** claim the model predicts physical membrane topology, microvilli
  density, or SM metabolite content from transcriptome. The permitted wording is
  "transcriptional program permissive-of / associated-with the topology
  phenotype."
- Do **not** present transcription as a substitute for the metabolite-level
  crosstalk; serine/SM flux requires metabolomics / the origin lab's single-cell
  MS. Transcription captures machinery/capacity only.
- Do **not** hard-code the direction of the tumor-serine → NK-SM crosstalk;
  report the sign calibrated on the liver positive control.
- Do **not** claim targets are experimentally validated.
- Do **not** claim graph attention is causal evidence.
- Do **not** present survival association as mechanism.
- Do **not** make any tumor-vs-NK axis claim on cell-type-unresolved data;
  such results are labeled MIXED_UNRESOLVED.

---

## 7. Strategic note — why this is the right first deliverable

(Internal note; not for the manuscript.)

This project is designed to function as a demonstration piece for a specific
target group whose flagship work is the surface-topology mechanism above. The
logic:

- **It speaks their language and builds on their finding**, rather than proposing
  an unrelated method. Recovering their own mechanism from independent data is
  the highest-trust signal a computational collaborator can send.
- **It offers leverage they lack.** Their mechanism was proven with specialized
  imaging/MS on few patients; a framework that reconstructs it at cohort scale
  and generates new, testable, cancer-specific targets multiplies the reach of
  every mechanistic paper they publish.
- **It targets an unclaimed extension.** The mechanism is already reported to
  extend to lung/colon/ovarian; gastric cancer (same digestive tract) is a
  natural next step not yet on that list.
- **It is transferable (the mechanism-card).** The offer is not one paper but a
  reusable engine: "each wet-lab mechanism you publish, I can turn into a
  scalable computational target-discovery run."

Two cautions kept in view: (1) the honesty boundaries in Section 6 are
load-bearing — a mechanistic lab will distrust any conflation of correlation with
mechanism or proxy with ground truth, so the disciplined qualifiers are an asset,
not a hedge; (2) a deliverable is a means of opening a conversation, not a
substitute for one — the document earns a meeting; the collaboration is built in
person.
