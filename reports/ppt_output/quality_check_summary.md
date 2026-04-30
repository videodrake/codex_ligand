# PPT quality check summary

Checked: 2026-04-30

## Output checked

- `reports/ppt_output/EGFR_MYO1D_receptor_modeling_clustering_deck_cluster12_selected.pptx`
- `reports/ppt_output/EGFR_MYO1D_progress_until_slide9.pptx`
- `reports/ppt_output/previews/slide01.png` through `slide09.png`

## Package inspection

- Slide count: 9
- Media count: 12
- Zero-byte media: 0
- Invalid PNG media: 0
- Final PPTX size: 2,199,124 bytes

## Visual inspection notes

- Slides 7-9 now reflect the selected Cluster 1/2 representative-structure plan.
- Slide 7 uses figures from `3.rotate_EGFR/10/cluster_results/70-100ns_cut0.173`.
- Slide 7 reports C1 as 121/301 frames, 40.2%, representative time 79.3 ns.
- Slide 7 reports C2 as 68/301 frames, 22.6%, representative time 93.5 ns.
- Slide 8 now separates completed work from MYO1D docking preparation.
- Slide 9 now moves from receptor-decision language to C1/C2-based MYO1D docking next steps.

## Tool note

The bundled `check_presentation_quality.js` script could not complete on this Windows environment because its internal `unzip -Z1` call failed. A PowerShell/.NET ZIP package inspection was used instead.

