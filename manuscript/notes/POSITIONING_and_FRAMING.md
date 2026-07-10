# GC-NKGraph-Atlas — Positioning & Manuscript Framing (revised)

> Replaces Section 0 (Project Identity) and Section 21 (Manuscript Framing).
> The change in this revision: the project is no longer framed as "a gastric
> cancer NK graph pipeline with a gated topology add-on" or as full reconstruction
> of a physical mechanism. It is framed as a **computational operationalization
> and transcriptional-reach map** for a specific, published immune-evasion
> mechanism, using liver as a partial-recovery control and gastric cancer as the
> extension setting.

---

## 1. Core positioning (one paragraph)

A landmark wet-lab finding (Zheng et al., Nat Immunol 2023, DOI
10.1038/s41590-023-01462-9) showed that tumors evade NK cytotoxicity by altering
the surface topology of NK cells: dysregulated tumor serine metabolism lowers
sphingomyelin in intratumoral NK membranes, which collapses membrane protrusions,
prevents lytic immune-synapse formation, and abolishes killing. That mechanism
was established with cutting-edge imaging and single-cell mass spectrometry on a
small number of liver-cancer patients. **This project asks which layers of that
mechanism can be surveyed at scale from public transcriptomes, which layers
require cell-type resolution, and which remain outside transcriptional reach**.
The framework uses a single-cell-informed, NK-aware heterogeneous graph to (a)
recover the protrusion-machinery → cytotoxicity effector layer in liver as a
partial-recovery control, (b) test the same recoverable layer in gastric cancer,
and (c) prioritize putative tumor-intrinsic candidate targets with wet-lab-ready
validation assays.

The central deliverable is **not** a prognosis signature, **not** a public atlas,
and **not** a claim that transcriptomes predict membrane topology. It is a
reproducible protocol that turns a mechanistic discovery into a bounded,
testable, scalable target-discovery workflow.

---

## 2. Two-arm design (this is the spine of the paper)

```text
Arm A — POSITIVE CONTROL (liver / HCC)
    The credibility anchor. Test how much of the published
    serine->SM->topology->cytotoxicity axis is visible in transcriptomes.
    Cohorts: TCGA-LIHC (+ optional HCC scRNA with intratumoral-vs-peritumoral NK).
    Result = partial recovery: effector layer yes, weak cell-resolved metabolic
    coupling, physical topology no.

Arm B — NOVEL EXTENSION (gastric cancer)
    The new contribution. Test whether the SAME axis operates in gastric cancer,
    a digestive-tract cancer NOT yet on the mechanism's published extension list
    (which already includes lung, colon, ovarian — see Clin Transl Med 2023,
    DOI 10.1002/ctm2.1395). Output gastric-specific tumor-intrinsic targets.
```

Why this structure matters: showing exactly which parts of a target lab's
published mechanism are, and are not, visible from independent public
transcriptomes is a higher-trust computational contribution than claiming full
recovery. It also resolves the previous tension between the "HCC-NK" heritage and
the "gastric cancer" focus: liver is the principled control, gastric is the
principled extension.

---

## 3. Recommended title

Primary:

**Mapping the transcriptional reach of the serine–sphingomyelin–membrane-topology
axis of NK-cell immune evasion: a single-cell-informed heterogeneous graph
framework from liver to gastric cancer**

Alternatives:

- **A single-cell-informed tumor–NK heterogeneous graph maps the transcriptomic
  reach of a physical immune-evasion mechanism**
- **Operationalizing NK surface-topology evasion: a knowledge-guided graph
  framework identifies recoverable and non-recoverable layers from liver to
  gastric cancer**

---

## 4. Core claims (conservative, mechanism-grounded)

State only what is supported by generated evidence. Each claim maps to a
figure/table in `manuscript/notes/main_claims.md`.

1. We define a reproducible, cell-type-attributed transcriptional proxy for the
   serine → sphingomyelin → membrane-protrusion-machinery → cytotoxicity axis of
   NK immune evasion.
2. In liver cancer (the mechanism's origin system), the framework **partially
   recovers** the axis: the protrusion-machinery → cytotoxicity effector layer is
   robust, the SM-balance → protrusion relationship is weak and detectable only
   after cell-type resolution, and the physical topology phenotype is not
   recovered from machinery transcription.
3. We construct a single-cell-informed tumor–NK heterogeneous graph in which a
   mechanism-grounded `metabolic_crosstalk` edge (tumor serine program → NK
   topology state) is justified by the biology, not by generic priors — giving
   the heterogeneous-graph architecture a principled reason to exist.
4. The graph model is statistically on par with the strongest tree baselines on
   NK-state prediction, outperforms weaker linear/kernel/shallow baselines, and
   provides mechanism-structured gene embeddings for axis analysis.
5. The recoverable effector layer is tested in gastric cancer; where it holds,
   we prioritize putative tumor-intrinsic candidate targets supported by
   multi-evidence patterns
   (tumor specificity, axis membership, NK-dysfunction correlation, spatial-niche
   support, ligand–receptor context, graph attention, druggability), each with a
   recommended wet-lab validation assay.
6. A model-based in-silico "SM-restoration" readout is retained as a hypothesis
   for experimental testing, not as a validated predictor or primary manuscript
   claim.

---

## 5. What is novel vs what is a control (state this explicitly)

```text
Positive control (not novel, but essential):
  - partial recovery and boundary mapping of the published
    serine-SM-topology-cytotoxicity axis in liver.

Novel:
  - scalable transcriptional reach mapping of that axis from bulk + scRNA;
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

## 7. Strategic note for authors — why this is the right first deliverable

This project is designed to function as a demonstration piece for a specific
target group whose flagship work is the surface-topology mechanism above. The
logic:

- **It speaks their language and builds on their finding**, rather than proposing
  an unrelated method. Recovering their own mechanism from independent data is
  the highest-trust signal a computational collaborator can send.
- **It offers leverage they lack.** Their mechanism was proven with specialized
  imaging/MS on few patients; a framework that maps its transcriptional reach at
  cohort scale and generates new, testable, cancer-specific hypotheses multiplies
  the reach of every mechanistic paper they publish.
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
