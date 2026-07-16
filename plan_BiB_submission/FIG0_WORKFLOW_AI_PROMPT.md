# AI image-generation prompt — Figure 0 (workflow / architecture overview)

## Why the previous single-prompt approach failed

The first version of this document asked one image-generation call to produce a
complete five-stage pipeline in a single image: five distinct icon clusters, a
strict left-to-right reading order, a fork into two branches that converge again,
per-stage colors, and no embedded text. That is far more compositional control
than current text-to-image models (Midjourney, DALL·E, Ideogram, Stable Diffusion)
can reliably hold in one generation — they are good at rendering *one* coherent
scene or icon, not at laying out six independent, precisely-labeled sub-scenes in
a fixed sequence with correct connective structure. In practice this produces an
image that looks atmospheric but is missing stages, merges two stages into one
icon, ignores the fork, or invents extra unrelated elements.

**The fix is the same one real BioRender/Nature graphical abstracts actually use:
generate each stage as its own small, simple icon image (models are reliable at
this), then assemble the six pieces yourself on a template canvas** in
PowerPoint/Figma/Inkscape/Illustrator, where you have exact control over position,
size, arrows, and labels. This is not a workaround — it is literally how BioRender
itself works (a library of individually-drawn icons a user arranges by hand), and
it is why the reference-paper schematics look precise: no one asked an image model
to generate the whole panel at once.

---

## Step 1 — generate six separate icon images

Generate each of these independently (one image-gen call per row). Each prompt
asks for a single self-contained icon/motif, square or near-square canvas,
transparent or white background, no text — this is well within what current
models render reliably.

Shared style suffix to append to every one of the six prompts below:

```
flat 2D vector icon, BioRender-style scientific illustration, single self-contained
motif centered on a plain white background, clean thin outlines (1.5-2pt), flat
solid fill only, no gradients, no drop shadows, no 3D bevels, no photorealism, no
embedded text or letters, no watermark, minimalist and print-clean at small size
```

| # | Icon prompt (the specific motif — combine with the shared suffix above) | Accent color to request |
|---|---|---|
| 1 | A small stack of 3 database-cylinder icons in a row, each with a faint horizontal bar-chart glyph on its face, next to one flat rectangular document/index-card icon with a small gear-and-checkmark badge in its corner | soft rose/pink `#F0DCE8` |
| 2 | A simple flat microscope icon on the left, next to a loose cluster of small circular dots in 3 muted colors forming a rough scatter blob (like a UMAP plot), with one small sub-cluster of dots subtly circled or highlighted | mustard/gold `#E8A93A` |
| 3 | An abstract network diagram icon: 6-8 small circles of 2 colors connected by thin straight lines, with exactly one connecting line rendered as a thicker dashed line in a highlighted color, the whole thing inside a soft rounded-rectangle outline | steel blue `#1F77A8` |
| 4 | A minimalist layered neural-network icon: 3 vertical columns of small circles connected by thin lines in a funnel shape, narrowing left-to-right into a single highlighted output circle | teal/emerald `#1B8A6B` |
| 5a | A single flat liver-organ icon, simple outline style, with a small circular checkmark badge overlapping its bottom-right corner | soft sky blue outline `#56B4E9` |
| 5b | A single flat stomach-organ icon, simple outline style, with a small circular magnifying-glass badge overlapping its bottom-right corner | green outline `#009E73` |
| 6 | A small cluster of 4 hexagonal "molecule" icons of slightly different sizes, loosely grouped, with one small syringe/pipette icon placed beside the cluster | burnt orange/vermillion `#D55E00` |

Generate 2–3 variants of each and keep the cleanest, simplest-silhouette result —
you want something that still reads clearly at 1–2cm on a printed page.

Style-keyword suffixes by tool:
- Midjourney: `--style raw --stylize 100 --ar 1:1 --no text,words,letters,3d,shadow,gradient,cartoon,photo,background`
- DALL·E / ChatGPT image / Ideogram: append `flat 2D vector icon, isolated on white, no embedded text, minimalist line icon, scientific illustration style`
- Stable Diffusion: negative prompt `text, watermark, signature, 3d render, photorealistic, gradient, drop shadow, cartoon character, multiple scenes, collage`

## Step 2 — assemble the six icons on a template canvas

Open a blank canvas in PowerPoint / Figma / Inkscape / Illustrator, page size
matching `fig0_workflow.pdf` (currently ~14 × 6.2in landscape — check
`src/figures/make_workflow_figure.py` for the exact current dimensions if you
want to match precisely). Place the six generated icons using this layout,
identical in structure to the current placeholder figure so no other files need
to change:

```
                [mechanism-card icon, #1]
                          |
                          v
[icon #1]  ->  [icon #2]  ->  [icon #3]  ->  [icon #4]  ->  [icon #6]
 DATA           scRNA          GRAPH          MODEL           TARGETS
                                                              /      \
                                                     [icon #5a]    [icon #5b]
                                                      Arm A          Arm B
                                                      liver          gastric
```

Concretely:
1. Place icon **#1** top-left as a small standalone card feeding down into the
   main row (this is the mechanism-card prior).
2. Place icons **#1** (data), **#2** (scRNA), **#3** (graph), **#4** (GNN), and
   **#6** (targets) left-to-right in one horizontal row, evenly spaced, each with
   a plain straight or gently-curved arrow connecting it to the next.
3. Below icon **#6**, branch two arrows down to icons **#5a** (liver) and **#5b**
   (gastric), placed side by side.
4. Add the text labels yourself (do not rely on the image model for these) using
   the exact wording in the table below — this guarantees correct spelling and
   lets you match the manuscript's typography.

## Step 3 — add labels (exact wording)

| Position | Bold label | Small caption underneath |
|---|---|---|
| mechanism-card icon | **mechanism-card** | machine-readable prior |
| Stage 1 | **Data + mechanism card** | *Phases 1–2* — TCGA-STAD, TCGA-LIHC, GEO gastric cohorts + machine-readable mechanism-card prior |
| Stage 2 | **Single-cell integration** | *Phases 3–7* — scRNA integration, NK-atlas annotation, 7-module SST-axis proxy scoring, trajectory |
| Stage 3 | **Heterogeneous graph** | *Phase 8* — PPI / ligand-receptor / TF edges + mechanism-grounded `metabolic_crosstalk` edge + SST-axis edges |
| Stage 4 | **Graph neural network** | *Phases 9–10* — baseline comparison + GNN gene-embedding model → NK-state classifier |
| Stage 5 (targets icon) | **Candidate targets** | *Phases 11–14R* — SST-axis scoring, de-circularized prioritization, 37 candidates, assay recommendations |
| Branch a | **Arm A — liver** | positive control (HCC) |
| Branch b | **Arm B — gastric** | extension cohort |

Title above the whole banner: **GC-NKGraph-Atlas: from mechanism card to
candidate targets**

## Step 4 — recolor to the manuscript's exact palette

Recolor each icon's fill/outline to match the Okabe–Ito palette already used in
`make_figures.py`, so the architecture figure and the results figures read as one
visual family — don't leave whatever approximate hue the image generator produced:
`sky #56B4E9`, `orange #E69F00`, `blue #0072B2`, `green #009E73`, `vermillion
#D55E00`, neutral grey `#8C8C8C`. Most vector editors let you select-by-color and
swap the fill in one step.

## Step 5 — export

Export the assembled canvas as vector (SVG/EPS) if your editor supports it (Figma
and Inkscape both do), then export/convert to `fig0_workflow.pdf` at the same page
size as the current placeholder so it drops into `main.tex` without changing the
`\includegraphics` width settings. If your editor only exports PNG, export at
≥600dpi at the target print size so line edges stay crisp in the two-column PDF.

## If you'd rather try a single-shot generation first (faster, less reliable)

Some newer tools (Ideogram, ChatGPT image with iterative editing) handle simple
layouts better than pure Midjourney/SD. If you want to try one shot before falling
back to the per-icon assembly above, use a **much simpler** ask than the original
version — three stages, not five, and no fork:

```
Flat 2D vector scientific infographic, BioRender style, wide horizontal banner,
white background, no gradients, no 3D, no photorealism, no embedded text.
Exactly three icon clusters in a horizontal row connected by two simple arrows:
(1) a small stack of database-cylinder icons, (2) an abstract network/graph icon
of small connected circles, (3) a minimalist layered neural-network icon. Flat
solid colors only: soft pink, steel blue, teal. Clean thin outlines, generous
white space, no other elements.
```

Expect to still need manual touch-up (adding the two remaining stages, the fork,
and all labels) — even this reduced version is a starting sketch, not a finished
figure. For a figure this structured, the six-icon assembly path in Steps 1–5
above is the reliable route to something you can actually submit.
