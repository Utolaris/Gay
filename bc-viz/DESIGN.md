# bc-viz design

## Concept
Instrument-panel viewer of MaleCNS activity under direct AN09B017b/c stimulation.
Not a marketing page: a research instrument. The living point-cloud brain is the hero.

## Palette
- --void #070A10 background
- --ink #E8EEF5 primary text
- --haze #6B7A8C secondary
- --line #1C2633 rules
- --stim #5EEBFF b/c stimulation (cyan)
- --drive #7CFFB2 excitatory drive
- --brake #FF6B9D inhibitory / mAL
- --p1 #FFD166 P1 readout

## Type
- body: system-ui, 'Segoe UI', sans-serif
- mono: ui-monospace, 'SF Mono', Consolas, monospace — metrics, IDs

## Layout
Full-bleed WebGL canvas. Left HUD: condition, time, legend. Bottom strip: population spike sparklines. No cards, hairline rules only.

## Signature
Activity wavefronts radiating from AN09B017b/c through the point cloud, with edge pulses only on the top b/c→target projection.

## Interaction
- orbit / zoom
- play / pause / scrub time
- toggle edge pulses, labels, dim inactive
- switch condition: bc_tonic vs WT (no drive)

## Tech
Bun + Vite + three.js. Offline static build. Prefer reduced-motion: auto-pause.
