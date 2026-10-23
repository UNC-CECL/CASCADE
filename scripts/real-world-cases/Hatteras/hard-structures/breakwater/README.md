# breakwater — a detached offshore breakwater (not started)

Notes from 2026-10-12, when the revetment was scoped first. Nothing is built.

## What a breakwater does

1. **Shelters the lee.** Wave height behind it drops by the transmission
   coefficient K_t (about 0 emergent, 0.3–0.7 low-crested or submerged).
2. **Grows a salient or tombolo.** Waves diffract into the shadow and drift
   slows there. Rule of thumb on length L_s over distance offshore X: below
   ~0.5–0.7 a small salient, ~1 a full salient, above ~1–1.5 a tombolo
   (Dally & Pope; Suh & Dalrymple; Hsu & Silvester).
3. **Acts as a leaky groin.** The salient traps drift; the coast downdrift erodes.
4. **Cuts storm erosion and runup in the lee.**

## A simple version in the model

- **BRIE:** scale wave height by K_t in the sheltered domains. Transport goes as
  roughly H^2.4–2.5, so drift drops there; the salient and downdrift erosion
  would follow from the drift gradient.
- **Barrier3D:** lower storm Rhigh/Rlow in those domains (the per-domain
  `StormSeries` setter, as in the revetment's deferred storm note).
- A real segment (100–300 m) is smaller than a 500 m domain, so, as the groin
  showed, only its net effect on drift can be represented.
