# 4D Polytope — Manual

Reference for the N-panel sidebar tab **"4D Polytope"**. Everything below is listed in the order it appears in the panel.

---

## Form

Sets *what kind* of 4D object gets built. Depending on the choice, the next 1–2 rows in the panel change.

### Form = Regular
- **Polytope** — one of the six regular convex 4-polytopes: 5-cell, Tesseract (8-cell), 16-cell, 24-cell, 120-cell, 600-cell.
- Only with this Form can **Dual Compound** and **Hopf Rings** further down be used at all (see below).

### Form = Uniform (semi-regular)
- **Symmetry** — which Coxeter reflection group is used: A4 (5-cell family), B4 (Tesseract/16-cell family), F4 (24-cell family), H4 (120-/600-cell family), D4 (branching diagram, Jenn3D's "Y-family" — gives the demitesseract and the omnitruncated demitesseract, among others).
- **Ringed Nodes** — a 4-digit bit code (e.g. "1000", "0110") for which of the 4 nodes of the Coxeter diagram are "ringed". Each of the 15 combinations gives a different uniform polytope (truncations, rectifications, omnitruncations of the base shape). Some single ringings correspond directly to the Regular polytopes (e.g. B4 + "1000" = Tesseract).

### Form = Hyperprism
- **Base** — the 3D solid that gets duplicated/extruded along the 4th axis. All 5 Platonic solids plus all 13 Archimedean solids are available (including the 2 chiral snub solids, built from their known closed-form coordinates rather than Wythoff ringing, which can't reach them). Note: these seeds only work with Struts/Ball and Stick/Wireframe, not Leonardo (more on that below).
- **Fourth-Axis Extent** — how far the second copy is shifted along the 4th axis (the prism's height in the W direction).

### Form = Duoprism
- **Sides P** / **Sides Q** — the number of sides of the two polygons whose Cartesian product is built (e.g. P=4, Q=4 gives exactly the Tesseract). Duoprisms sit on a Clifford torus of the 3-sphere and look especially good with curved edges (nested rings of circles).

---

## Edges

- **Edges (Curved/Straight)**
  - *Curved (stereographic)*: vertices sit on the 3-sphere, edges are great-circle arcs, projected exactly stereographically.
  - *Straight (perspective)*: a direct 4D perspective projection; at a small distance the result approaches a Schlegel diagram.
- **Projection Distance** *(Straight only)* — eye distance along W. Values near 1.0 produce a Schlegel diagram.

---

## Style (render mode) — and its parameters

Sets *how* the projected edge/face framework is actually built.

- **Struts** — solid square-section tubes along the edges, no vertex markers. Fastest.
- **Ball and Stick** — like Struts, plus a marker at every vertex.
- **Wireframe** — plain polylines with no volume, the cheapest option for very large/animated objects.
- **Leonardo (da Vinci)** — flat, framed panels for each 2D face of the polytope (an open "window" look). Only works on polytopes with a single, uniform face size (not on the Archimedean hyperprism seeds, not on Uniform/Duoprism).

Depending on the Style, different follow-up parameters appear:

**With Struts / Ball and Stick:**
| Parameter | Effect |
|---|---|
| Arc Segments | Number of segments per edge (more = smoother arcs under Curved style, but more geometry). |
| Strut Radius | Thickness of the edge tubes. |
| Strut Sides *(Ball and Stick only)* | Number of sides in the cross-section of each round strut. |
| Taper With Projection | Scales strut thickness by the local projection factor (areas near the projection pole render thicker — more realistic, optionally switched off). |
| Sphere Size *(Ball and Stick only)* | Size of the vertex markers relative to the strut radius. |
| Vertex Shape *(Ball and Stick only)* | **Sphere** (round, ~48 faces) or **Cube** (6 faces, roughly 8x cheaper — recommended for large objects like the 120-cell and for animation). |
| Scale | Overall size of the object. |

**With Wireframe:** just Arc Segments and Scale.

**With Leonardo:** Border (frame width relative to face size), Panel Thickness, Taper With Projection, Scale.

---

## Rotate XW / YW / ZW / XY

Four independent 4D rotation planes, applied **before** projecting down to 3D, in degrees. These are the only values worth animating (keyframes work, see below). Worth knowing:

- The rotation is applied **after** any Half Cutaway — a 180° turn can rotate a vertex that was previously "cut away" right back into view.
- With curved projection there's also an automatic, invisible correction rotation that keeps any vertex from landing exactly on the projection pole (which would blow up a single strut to a huge size). That happens on its own in the background and isn't an adjustable value.

---

## Cutaway, Compound & Rings

- **Half (Cutaway)** — keeps only the elements on one side of the equatorial hyperplane (W ≤ 0, equator included), **before** projecting and rotating. This keeps the result stable no matter how you rotate afterwards (the rotation then spins the already-cut remainder further, see the note above under Rotate).
- **Dual Compound** *(Form = Regular only)* — also builds the dual polytope (vertices at the cell centers of the original), in a second material slot.
- **Hopf Rings** *(Form = Regular, Polytope = 120-cell only)* — renders N of the 12 rings of 10 dodecahedral cells that partition the 120-cell along its Hopf fibers, as solid, shrunken cells with one material per ring. (Grayed out otherwise — hover for a tooltip explaining why.)
  - **Ring Cell Scale** *(only when Hopf Rings > 0)* — shrink factor of each ring cell toward its 4D center (gaps make the ring structure easier to read).
  - **Rings Only** *(only when Hopf Rings > 0)* — hides the normal edge framework and shows only the solid ring cells.

---

## Materials (automatic, no UI control)

- **Struts/Wireframe/Leonardo, without Dual/Rings:** no extra material needed.
- **With Dual Compound and/or Hopf Rings:** "Polytope4D Primal", "Polytope4D Dual", "Polytope4D Ring 1…12" — one slot each.
- **Ball and Stick:** always two simple placeholder materials, "Polytope4D Sticks" and "Polytope4D Spheres" (regardless of Dual/Rings) — ready for your own styling later.

---

## Bottom section

- **Auto Update** — when on, the mesh rebuilds immediately on every parameter change, including keyframe playback/scrubbing. On large objects (120-cell, many Arc Segments) it's worth switching this off during setup and using **Rebuild** instead, to avoid stutter.
- **Rebuild** — rebuilds the mesh once, manually. Needed when Auto Update is off.
- **Last build …s — N vertices, M edges** — build time and the actual vertex/edge count of the last build (the main body only, not the dual or rings).
- **Reset to Defaults** — resets every parameter to its starting value and rebuilds.
- **Detach from Engine** — freezes the current mesh and hides this panel for the object (the mesh itself is left untouched, only the live link is removed). After this, the object is a perfectly ordinary, static mesh.

---

## Animation (short version)

Rotate XW/YW/ZW/XY are ordinary Blender properties — right-click → **Insert Keyframe** works directly. Since keyframes don't trigger the sliders' own update callback (a general Blender quirk), a `frame_change_post` handler runs in the background and rebuilds every engine object with Auto Update on whenever the frame changes (both during playback *and* manual scrubbing). For smooth scrubbing on large objects: switch off Auto Update while working on keyframes, use **Vertex Shape = Cube** instead of Sphere, and/or lower Arc Segments.
