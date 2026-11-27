# Tracking for planning/writeup_plan.tex

## Delivery dates (e9 work spec `work_spec_e9_2026-11-26.md`, sent and accepted 26 Nov; figure IDs are the new ones)

| Item | For | Owner | Due | Output |
|---|---|---|---|---|
| N1 BR and bubble coherence (V, V+, V−, V+ − V−) | F6 | b0 | 26 Nov | Elements `asme2027_figs/<item>/` |
| N5 coherent share of μ_eff, 193M | 4.1 (ii) text | a3 | 27 Nov | same |
| N3 KH centre-plane instantaneous view | F8 | b0 | 27 Nov | same |
| 193M ε and per-batch fields on the W1 point set | F14 SE, E13 check | fc | 27 Nov 18:00 | Elements |
| N3 mean cross-planes x/D 2, 4, 10; N4 Q/λ2 volume | F8, F4 | b0 | 28 Nov | same |
| N2 wall C_f space–time and snapshots, x/D 0–14.9, 193M + 20M LES | F10 | a3 | 28 Nov | same |
| N6 slow mode and film meander, 193M | 4.1 (i) text | a3 | 29 Nov | same |
| N2 extension to x/D 25 (Elements read speed is the risk) | F10 | a3 | 30 Nov | same |
| N7 wall-layer coolant balance on the 193M b010 (owner: strong physical evidence if it holds; decide after) | 4.1 or 4.3 | fc | 28 Nov 18:00 | `asme2027_figs/wall_balance/`, FINDINGS WL4+ |
| W1 item 3 scale-resolving totals | 4.2 bound sentence | fc | 1 Dec | |
| W2 directional diffusivities, 193M and 100M | F19 | fc | 2 Dec | |
| STAR base, QCR, QCR + Pr_t 0.4 gradient exports | 4.3 a⊥, F4, F8 RANS side | STAR session | running since 26 Nov 11:10 | Elements `rans50_exports/` |

### Owner rules of 26 Nov to apply in the 28 Nov pass

- **193M is b010 only:** remove 193M b006 values from the plan text (F19 "b006 → b010", near-wall Pr_t "b006
  against b010", F5 "b006 is within", 4 Results common bullet) and from `fig_prt_data.py`; convergence statements
  use per-batch b010 data.
- **Domain ends at x/D 20** (x/D 20–25 is a buffer; nothing reported there). Plan lines to change (line numbers
  of 26 Nov): 396 region definition "far film (to 25)"; 669 "laminar to x/D = 25"; 724, 732 wall-state plan view
  to x/D 20; 764 "0.043 at 25"; 803 η map 0–25; 931, 958, 1009 far film 10–25; 982, 1144, 1158, 1252, 1257
  "lateral η̄ x/D 20–25" → x/D 15–20 (STAR to recompute for every run); 1065 centre-plane figure to x/D 20.
  Numbers to recompute, not relabel: μ_t split far film (estimator on 10–20; was effectively 10–20.4), volume
  Pr_t far film (fc), lateral η̄ of every a posteriori run and of the 193M (STAR), figure axes (fig_wall_state,
  fig_eta_wall, fig_eta_x, fig_centre_factors, fig_bl_params already to 20, fig_spots).
- **Volume Pr_t in the paper and abstract** (applied 26 Nov; far-film value pending the 10–20 recompute).
- **Every delivery as soon as possible** (e9 told).

### Received 26 Nov, for the 28 Nov revision of 4.2 (not yet in the plan)

- **Volume classic Pr_t,R of the 193M (fc, W1 item 4, 1eb725f, FINDINGS C26):**
  `4dns_filmcooled/analysis/eps_prt_c193b/results/prt_classic.txt`. Volume least squares, z/D ≥ 0.3, point-wise
  regions. b010 (b006): jet 0.300 (0.307), near film 0.225 (0.229), far film x/D 10–25 0.471 (0.481). No
  jackknife SE; b006 → b010 shift 0.004–0.010. Replaces the centre-plane 0.44 / 0.63 / 0.80 (plan F15 and step
  2, abstract "0.44 in the jet and 0.6–0.8 in the film"). Consistent with the x/D 10 cross-plane value 0.33.
  All three are below the free-shear 0.5; step 2 wording ("near the free-shear value in the jet") and the 4.3
  reading of the Pr_t sweep (saturation at 0.3–0.4) to be revisited. Far-film region differs: 10–25 here, 10–20
  on the centre plane.
- **N6 slow mode and film meander, 193M b010 (a3, 5a26246):** `/mnt/Elements/DNS/analysis/asme2027_figs/n6_slow_mode/`
  (`summary.txt`, `c193b_slow.npz`). Lobe difference V+ − V− peaks at St 0.037–0.046; lobes in antiphase
  (coherence 0.50 at St 0.039, −169°); total V has no peak there. Footprint shift Y peaks at St 0.02–0.05.
  Coherence of the lobe difference with Y: 0.28 at St 0.039 for x/D 3–8 (95 % level 0.206 with 14 overlapping
  segments, 0.312 with 9 non-overlapping), 0.16 for x/D 8–15; V–Y at most 0.19. Marginal. The 50M "79 % of the
  slow-mode energy in the film" cannot be repeated without a 193M SPOD; substitute: wall η has 0.43 of its variance
  in St 0.02–0.06 over x/D 3–15. Frequencies below St 0.03 not resolved. For the 4.1 (i) "effect on the
  mainstream" bullet.
- **N2 spots (a3, 99eff14, FINDINGS SP1–SP5, cc72e97):** wall C_f data x/D −3 to 14.9 (near-wall volume ends at
  14.94). **Owner decision (via a3, 26 Nov):** celerities and merger statistics come from the du/dt method to
  x/D 20 (193M: 0.46 leading, 0.40 trailing, N 18, computed to x/D 20.5; recompute to 20); wall C_f is used only
  for the space–time and plan-view figure. Note: `fig_spots.py` now draws the C_f version; the du/dt version is
  in git (`e0913c1:figures/code/fig_spots.py`). Candidate `fig_spots_z010.png` (C_f estimated from u at z/D 0.1,
  to x/D 20, calibrated against the wall C_f on the overlap) is being built.
- **N5 coherent share of μ_eff, 193M b010 (a3, d9d5f22, FINDINGS CL1–CL5):** group medians, Lav σ 0.026 | σ 0.036:
  hole 0.030 | 0.045, jet 0.055 | 0.077, film x/D 3–20 0.043 | 0.059 (far film 10–20: 0.043 | 0.061), approach
  layer 0.002 | 0.005. Heat (Γ_κ) film 0.048 | 0.070. k share: jet 0.071 | 0.095, film 0.052 | 0.073. Replaces the
  50M "0.06–0.12". Paper value: coherent share 0.03–0.08 of μ_eff in the hole, jet and film (state the filter
  width). No strain cross-check for the 193M.
  **Owner decision 26 Nov (via a3, FINDINGS LV1–LV3a, f6bcddb):** quote the σ 0.036 values, with the Lav-σ values
  alongside: hole 0.045 (0.030), jet 0.077 (0.055), film x/D 3–20 0.059 (0.043), far film 10–20 0.061 (0.043).
  Rung comparison at σ 0.036: 193M within 0.013 of the 50M in the jet and hole; wake x/D 7.5–10 0.037 against 0.075.
- **Spots figure (now combined `fig_spots.png`):** top z/D 0.1 C_f estimate to x/D 20, bottom wall C_f snapshots
  to 14.9. Caption numbers: z+ ≈ 8 laminar / 12.5 turbulent; estimate/true 0.94 laminar, 0.64 turbulent (193M);
  z/D 0.1 leads the wall by 1.2–2.5 D/U. du/dt celerities at x/D ≤ 20: 193M 0.46 / 0.40 (N 18), LES 0.41 / 0.33
  (N 6); mergers N 2 per case (describe, do not quantify). Plan entry F10 to be rewritten.
- **N8 span-summed Pr_t sums (fc, received 26 Nov):** `/mnt/Elements/DNS/analysis/asme2027_figs/n8_lsq/n8_lsq_c193b_b010.npz`.
  Candidate `fig_prt_classic_span.png` drawn (map + profiles against z+, log axis; u_τ 3.84 m/s near film, 4.31 m/s
  far film). Owner: replaces the centre-plane F13 once the per-batch band is in (per-batch sums due 27 Nov
  10:00–12:00). Span values: near film 0.80 at z/D 0.1, 0.39 at 0.2, 0.18 at 0.3; far film 0.98, 0.89, 0.83.
  Far film: fitted heat diffusivity negative at z/D 0.02–0.08 (z+ 2.7–10.7), counter-gradient wall-normal heat
  flux; the claim rests on few bins until the band exists. Step 2 is to mention the far-film layer near 0.9
  (owner, 26 Nov). Owner: the paper and abstract use the volume values.
- **193M ε point file (fc, received 26 Nov 18:24):** `/mnt/Elements/DNS/analysis/asme2027_w1/eps_points_c193b_w1_b010.npz`
  (b010 only; 3,141,253 points; 10 delete-one time-batch fields). Far film to be cut at x/D 20. In use for the
  time-batch SE of the μ_t split (`planning/maths_audit/estimator/estimator_timejk.py`) and the dissipation-length
  check.
- **STAR x/D 20 window:** all tables and `fig_eta_x.py` use x/D 19.5–20.0 (193M centreline 0.387, mid-pitch 0.024–0.025;
  RANS 0.497 / 0.006).
- **η map (`fig_eta_wall.png`):** x/D −2 to 20; smooth hole outline with border and label.

### 4.3 plan (written 26 Nov) and run status

Plan names → STAR runs: "Pr_t 0.4 near the hole" = RC-S10 local; "Pr_t 0.4 in the film only" = P1, RC-S14, job
39592916 (24 h from 26 Nov evening); "3DNS heat flux" = run C, RC-S13; "3DNS heat flux in the film only" = P2,
RC-S15, job 40253285; "3DNS stresses" = run A (smoke failed 26 Nov, guard retry 41747008 pending); "3DNS
stresses and heat flux" = run B, job 39747096 running; "γ = 1" = RC-S9 run F 50M, job 40253286, results about
27 Nov 19:00. Numbers: `4dns_filmcooled/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns_substC/compare/substC_eta_windows.txt`.
STAR pre-run predictions: PLAN lines 96–99, 116–122, 674–675; sub-additivity A + C ≥ B (PLAN 127–133).
Open for the owner:
- Runs A, B, C and P2 impose stresses and heat flux from the **50M** 3DNS (C50b b006), not the 193M. The η
  references are the 193M. **Decided 26 Nov (owner):** stated plainly as the 50M in 4.3, no further comment; for
  the post-review revision and the thesis the tables are rebuilt from the 193M b010 and the runs rerun (e9 told).
- Run A: owner confirmed 26 Nov that it stays; 4.3 keeps runs A, B, C and the A + C ≥ B test.
- No run changes the length scale alone (stated in 4.3 (vi)).
- F20 (`fig_aposteriori.py`) not drawn; the lateral profile exists for run C as STAR's `fig_substC_eta_lateral.png`.
- QCR decomposition: Sonnet agent running on the STAR gradient exports (`planning/maths_audit/qcr/`).
- Decided 26 Nov: bubble figure shows total V; page-16 note meant H and C_f; 50M bubble numbers dropped.

**Dates revised 26 Nov after the owner's ASAP instruction (supersede the table above):** ε file today ~17:00;
N5 today 18:00–19:00; N8 cumulative b010 today ~18:30; N3 planes and N4 tonight ~22:00 (latest 27 Nov 10:00);
N8 per-batch and N7 27 Nov ~08:00; W1 item 3 27 Nov ~12:00; W2 27 Nov ~14:00; N2 x/D 0–14.9 after a
spot-celerity check, extension to x/D 20 end of 27 Nov (latest 28 Nov midday); G4 time-scale ratio on the 193M
28 Nov ~10:00. The x/D 20 rule and the b010-only rule have been passed to fc, b0, a3 and the STAR session.

Figure freeze 4 Dec. The owner's plan freeze was lifted on 26 Nov; the earlier freeze note is void.
Open check: the page-16 note of the 26a PDF ("this should probably be first ... a bit small and hard to read")
sits on the heading of the old F11 (profiles), directly under the old F10 (H, C_f). It was applied to the H/C_f
figure (now F7: first in the mainstream chain, enlarged, text-and-arrow labels); owner to confirm.

## Update 2026-11-26a (owner feedback on the annotated plan)

Backup of the plan before this update:
`/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/writeup_plan_backup_pre26a.tex`.
Old figure IDs in the sections below this one are the pre-update IDs; translate with this table.

### Old ID to new ID

| old | new | figure / change |
|---|---|---|
| F1 | removed | transition map; justified in text through Kar et al. (2026) and Clark et al. (2024) |
| F2 | F1 | config and regions |
| F3 | F2 | mesh classes |
| F4 | F3 | resolution |
| new | F4 | `fig_overview_3d.png`, 3D view RANS 50M against 193M (not drawn) |
| F5 | F5 | hole axis |
| F6 | removed | diffuser flow; numbers kept in one sentence in 4.1 (i) |
| F7 | F6 | bubble motion, reworked |
| F10 | F7 | `H` and `C_f` along the plate, now first in the mainstream chain |
| new | F8 | `fig_flow_features.png` (not drawn) |
| F8 | F9 | wall state |
| F9 | F10 | spots, reworked |
| F11 | removed | profiles at x/D 10; one sentence with its numbers in 4.1 (ii) |
| F12 | F11 | eta against x/D |
| F13 | F12 | eta on the wall |
| F14 | F13 | classic Pr_t, panels (a) and (e) kept, renamed (a) and (b) |
| F15 | F14 | mu_t split |
| F16 | F15 | regional Pr_t |
| F17 | F16 | cross-planes |
| F18 | F17 | centre-plane factors |
| F19 | F18 | budget lines |
| F20 | F19 | spanwise diffusion |

In 4.2 only figure numbers changed, plus: "in every run" to "in every Pr_t and QCR run" (mid-pitch eta
0.006-0.012); the reference "F14 (e)" is now "F13 (b)" because of the new panel letters.
Run C (RC-S13) is not added to that sentence; it belongs to 4.3.

### Requests for e9 (awaiting owner approval)

1. Bubble figure F6: BR and bubble signal on one time axis, and the coherence between them, on the long continuous
   stretch of the 193M record. Candidate signals: total V, each lobe, lobe difference V+ - V-; the strongest is
   used. Report the stretch length, number of segments n and the 95 % level.
2. Effect of the bubble on the mainstream, on the 193M b010 (4.1 (i), last bubble bullet): share of the slow-mode
   energy in the film; coherence of the lobe switching with the film meander and the sideways shift of the cooled
   footprint. The 50M values that stood in the plan are kept here: slow mode at St 0.04 has 79 % of its energy in
   the film; coherence <= 0.10; bubble motion modulates the coolant flow by about 1 %.
3. Periodic (coherent) share of mu_eff on the 193M (method of Lav et al., 2019); the plan states 0.06-0.12 in the
   hole, jet and film, to be recomputed (4.1 (ii), F8 WHAT; the same number is in 4.2).
4. Spots figure F10: time-resolved wall C_f (z = 0, or the first cell centre, plane to be stated) on the plate for
   the 193M b010 window and for the 20M LES, along the mid-pitch line and in plan view, at the record rate of the
   existing z/D = 0.1 sampler (193M: 3992 records over 400 D/U). Spot-edge celerities to be re-estimated on C_f.
5. F8: mean streamwise vorticity and in-plane streamlines on x/D 2, 4, 10 for the 193M (RANS 50M exists); one
   instantaneous 193M field (Q or T) showing the Kelvin-Helmholtz roll-up, with the time chosen.
6. F4: 193M mean Q (or lambda_2) field on the full domain for the isosurface, wall eta, hole outline; same for
   RANS 50M; camera and cut-plane frames (centre plane, x/D 4, 10, 20) identical.

Scan for other comparisons that use the 50M or another rung as the reference (numbers not deleted; for the owner to
decide):

- 4.1 (i), first bubble bullet: "10.9 % of the diffuser on the 50M" and the 50M volume 0.761 beside the 193M
  value; the 50M lobe correlation -0.59 in F6 WHAT.
- 3.1 (ii): "RANS against 50M 3DNS on the same mesh isolates modelling; 50M, 100M and 193M isolate resolution":
  a design statement; the 50M 3DNS is a comparison rung there, not the reference of any result.
- 4.2 step 1 and F14/F19: the RANS 50M is compared with the 193M (correct). The "scale-resolving rungs lie within
  +-2 SE of zero" sentence is also against the 193M.
- 4.2 sweep bullet: RC-S9 gains (5M +1.3 %, 50M smoke test +1.26 %) are against the unmodified RANS of the same
  mesh, not against the 193M.
- 4.1 (iii): "the scale-resolving rungs agree within 0.008" and the 20M, 50M, 100M values are against each other,
  with the 193M as reference.
- Removed from the plan with the F7 rework (kept here): coherence of V with BR on the 193M 0.78 at St 0.273 and
  0.39 at St 0.078 against 0.21 (n = 14 segments of 25.6 D/U, 50 % overlap, inside 13 gap-free stretches); 50M
  0.82, 0.35, level 0.24; long-segment estimate (256-record segments, n = 4, level 0.63) V peak St 0.293,
  coherence 0.90 at St 0.293 and 0.85 at St 0.273, no peak near St 0.04. The 1995 snapshots at 0.2 D/U
  (99.97 % of b010) described the old estimate.

### Pending

- Reference for the coherence significance level (1 - 0.05^(1/(n-1)); F6) needed. No citation invented.
- Number needed from the cited papers: the error in laterally averaged effectiveness of best-practice RANS for
  shaped holes (Introduction, P1, third bullet). The bullet is written without a number.
- Intro P5a: "RANS under-predicts the lateral spreading they set" reuses the claim of P4 (Acharya et al., 2001;
  Harrison and Bogard, 2008a); check that the cited papers support it for mean vortex features.

### Restyle notes for the figure agent (not applied; figures and scripts untouched)

- F7 (`fig_bl_params`): enlarge (plan width 0.95); Blasius and turbulent correlations, and any other line that is
  not a case, are labelled by text and arrows in the plot, not in the legend; the legend lists the cases only.
- F13 (`fig_prt_classic`): keep only the centre-plane map (old panel a) and the x-averaged Pr_t against z/D (old
  panel e); drop old panels b, c, d; relabel (a), (b); the x/D = 10 medians 0.72 and 0.65 go to the text.
- F10 (`fig_spots`): space-time diagram of C_f (not du/dt); 20M LES beside the 193M; hole position marked in
  every panel; shared x axis and marks between top and bottom panels; plane or height stated in the caption.
  The snapshot quantity in panel (c) is assumed to be C_f (owner note covered the space-time diagram only).
- F6 (`fig_bubble_motion`): two panels (BR and bubble signal in time; one coherence line with the 95 % level as one
  clear line and the exceedance frequencies marked).
- New scripts: `fig_overview_3d.py` (F4), `fig_flow_features.py` (F8).

### Other

- "Strongest" bubble signal (F6) is read as the signal with the highest coherence with BR; the owner note did not
  give the criterion.
- Intro: the in-hole literature bullet (Thole 1998, Gritsch 1998b, Haydt 2017/2018, Oliver 2017) was removed from
  P5. Thole 1998 moved to P5a; Haydt 2017/2018 to P8; Gritsch 1998b and Oliver 2017 stay cited in 4.1 (i) only.
- Intro paragraph P5 is now P5a (mean features) and P5b (instantaneous features); P6 to P8 keep their numbers.
- The aims were reordered; the old aims 3 and 4 (locate and explain; identify which matter) are merged into the
  new aim 2.
- `\setcounter{figure}{1}` before Section 3 is now `{0}` (F1 is the config figure); LaTeX figure numbers match the
  F IDs.

**Plan frozen for owner review from 2026-11-26 evening; the owner lifts the freeze at about 15:00 on 2026-11-27.
No edits to writeup_plan.tex until the owner says so.** Incoming data during the freeze (fc 193M ε and per-batch fields, 27 Nov 18:00; STAR ω/limiter and QCR
exports) are checked and recorded here, not applied to the plan.

## Received during the freeze (not yet in the plan)

- **Run C, RC-S13 (STAR session, 26 Nov):** SST stresses with the 50M 3DNS heat flux imposed in all regions.
  FINDINGS "Step 4 run C" SC1–SC6; `4dns_filmcooled/starccm/rans_50M_Lm25b_BR15_FST_ZPG_p0dns_substC/compare/substC_eta_windows.txt`.
  - Converged (SC1). Recovered x/D 7–20 η̄ deficit against the 193M b010: 76.3 % (QCR 72.1 %, Pr_t 0.4 44.1 %).
  - η̄ windows 2.5–7 / 7–15 / 15–20: 0.292 / 0.220 / 0.191 (193M 0.304 / 0.237 / 0.208; SST 0.272 / 0.168 / 0.138).
  - Centreline η at x/D 20: 0.335 (193M 0.383, SST 0.494). Mid-pitch η at x/D 20: 0.0365 (193M 0.0256; every
    other RANS run 0.006–0.012). Lateral η̄ x/D 20–25: 0.182 (193M 0.194).
  - Run C is the only RANS run with coolant at mid-pitch; it over-spreads (centreline low, mid-pitch high).
    STAR's reading (SC6, medium): the window mean overstates the match. Runs A (stresses) and B (both) pending.
  - For 4.3 and the plan's "mid-pitch η stays 0.006–0.012 in every run" sentence (true for the Pr_t/QCR runs only).
- **e9 rundowns (72f6d5b):** `4dns_filmcooled/audit/supporting/rundowns/R26_WALL_LAYER_BALANCE.pdf`,
  `ALL_RUNDOWNS.pdf`. Wall-layer coolant balance (FINDINGS WL1–WL3): in the 3DNS the turbulent flux of coolant
  through z = 0.5 D exceeds the mean upwash from x/D ≈ 4; the RANS mean upwash is 2× (x/D 5) to 11× (x/D 10) the
  3DNS and its turbulent flux 0.4× at x/D 2. Balance accept test not passed (residual up to 16 %). Not read in
  full; candidate evidence for 4.1/4.3.
- **STAR gradient exports:** base export running since 11:10; paths to follow.

Companion to the plan. The plan holds the paper's content only. This file holds status, dates, owner decisions,
pending work, provenance, names and open questions that were removed from the plan on 2026-11-26 (owner
instruction: plain plan; everything else in a separate document for the agent and Claude). Items are grouped by
plan location. Category keys: S status, D date/owner decision, P pending/due, V provenance/source, N agent/session/
person/commit/finding/job id, F file path of analysis data, O open question, H history/superseded.

Backup of the plan before the rewrite:
`/tmp/claude-1001/-home-sk999-Documents-ASME-2027/f4c92ceb-6df8-4254-96c3-6552d15d0ecf/scratchpad/writeup_plan_backup_2026-11-26.tex`

## Front matter and lists

- [D] Title block carried the date "26 November 2026".
- [P] Intro text: "Not yet planned: 4.3 A Posteriori Tests (after the P1, P2 and RC-S9 runs) and 5 Conclusions (after 4.3)."
- [V] "Each agreed section is written into manuscript/main_claude.tex."
- [S] Status legend: exists = usable as it is; rework = content right, data or layout to be updated; new = not yet
  made; draft = made for this plan, to be finalised; interim = drawn from the data available now.
- [S] Figure list status column (dropped):
  F1 exists; F2 draft; F3 draft; F4 interim; F5 rework; F6 rework; F7 interim; F8 rework; F9 new; F10 new;
  F11 rework; F12 rework; F13 new; F14 new; F15 rework; F16 new; F17 new; F18 rework; F19 new; F20 rework.
  (Entries in the body carried finer tags: F7 "new", F14 "rework", F5-F20 "Complete." in captions; see below.)
- [S] Table list: T1 draft. Table list "source" column: "mesh files, hole_discharge" [F].
- [S] Equation list column "in the draft" (dropped): E1 new; E2 new; E3 numbered (eq:re); E4 numbered (eq:br);
  E5 inline, numbered in the paper; E6 numbered (eq:N); E7 numbered (eq:xi); E8 numbered (eq:eta), changed to
  eta_K,193; E9 methods eq:meff; E10 new; E11 methods eq:three; E12 methods eq:split; E13 new (2026-11-26);
  E14 methods appendix A.

## 1 Introduction

- [H] Purpose: "The RANS failure leads the argument" reworded to "The RANS errors lead the argument".
- [H] P4 "the Boussinesq relation fails in the hole ..." reworded to "does not hold" (neutral voice).
- [V] P4 equations "(new; standard forms, ...)": the word "new" removed.
- [F] P6: values read from the papers, `4dns_filmcooled/write_up/literature/papers/`.
- F1 entry:
  - [S] tag "exists", "reused".
  - [V] "Cut from the published PDF as a vector crop; the caption cites both papers."
  - [V] "The one paper figure not drawn with set_mpl_params; its colours are those of the original." (the colour
    remark stays in HOW).

## 2 Research Aims

- [P] "Wording of the aims is set with the Results plan; the conclusions mirror them." Plan keeps "The conclusions
  mirror the aims."

## 3.1 Configuration

- [F] "(from 01_configuration.tex; E3, E4 numbered there, E5 inline there)".
- [V] E5 purpose: "Numbered in the paper, since the Results refer to it."
- [S] T1 tag "draft"; caption "Draft of T1."
- [F] T1 caption: CGNS files in `/mnt/Elements/DNS/meshes/`; BR source `analysis/hole_discharge/README.md`.
- [S] F2 tag "draft", "full page width" kept as layout; caption "Draft of F2".
- [V] F2 provenance: "Okabe-Ito fills" kept in WHAT.
- [P] F2 "To finalise: hole outline from cl_outline.npy [F]; label sizes checked at the full page width of the
  ASME layout."

## 3.2 Mesh family

- [F] "(from 03_mesh_resolution.tex; the decisions in HANDOFF.md are kept)"; "Equations (from 03_mesh_resolution.tex,
  both numbered there)"; "In text, no table (replaces the class table of 03_mesh_resolution.tex)".
- [S] F3 tag "draft"; caption "Draft of F3".
- [P] F3 "To finalise: class numbers are small at print size."

## 3.3 Achieved resolution

- [F] "(from 03_mesh_resolution.tex, changed to the single reference eta_K)"; "(replaces the resolution table of
  03_mesh_resolution.tex)".
- [H] "The lower-bound caveat of the draft (filtered or modelled dissipation on the 20M LES and 50M 3DNS) no
  longer applies, and no LES statistics are needed for the 20M mesh."
- [H] "All numbers from the new computation (above); the draft values (batches b005/b006, own eta_K, placeholder 20M)
  are not used."
- [P] "Method to compute it (for the analysis agents)": `dns/python/resolution_regions.py` uses each run's own
  epsilon and "is to be extended" [F].
- [S] F4 tag "interim".
- [P] F4: "The interim version below uses each run's own eta_K and the histogram regions x/D 3-10 and 10-25. The
  final version uses eta_K,193 for every mesh, adds the 5M curve and uses the region bounds of F2." Caption said
  "Interim F4: each run with its own eta_K." The plan now describes the final version (eta_K,193, five meshes).
- [O] To state in 3.3 text: "To state:" wording on eta_K,193 uncertainty in the diffuser converted to a plain fact.

## 4 Results (common)

- [P] "Section 4.3 (a posteriori tests) is planned once the tests are done."
- [F] "All 3DNS values come from the stat3 files ...; RANS values from the STAR cell exports." Plan says
  "time-averaged statistics" and "cell exports".

## 4.1

- [V] (i) "Effect on the mainstream (50M analysis, not yet repeated on the 193M)": the pending repeat on the 193M.
- F5:
  - [S] tag "rework"; caption "Complete."
  - [F] Source: `analysis/hole_discharge/plots.py` F3; "redrawn to the figure rules, 5M RANS and 193M added".
  - [F] caption "stat3 b010".
- F6:
  - [S] tag "rework"; caption "Complete."
  - [F] Source: `hole_discharge` F4, F5.
- F7:
  - [S] tag "new" (body), "interim" (list); caption "Complete."
  - [N] Source: `hole_discharge` W3 (b0, commit 381ca65; FINDINGS S15, S16, S18). The S17 mrDMD numbers are not used in the paper.
  - [V] 1995 snapshots at 0.2 D/U, 99.97 % of the b010 window (kept in HOW).
- F8:
  - [S] tag "rework"; caption "Complete."
  - [F] Source: `layer_state` M1 (gamma maps); k maps new from stat3 and cell exports.
- F9:
  - [S] tag "new"; caption "Complete."
  - [F] Source: new, from `layer_state` (Otsu detector) applied to the 193M z/D=0.1 sampler.
- F10:
  - [S] tag "new"; caption "Complete."
  - [F] Source: stat3 means and RANS exports (approach-layer values in `setup_differences`, SU3, SU7; plate values in `layer_state`, LS3); RANS `bl_summary.csv`.
  - [U] Panel layout: the entry said "H (a) and C_f (b)", the caption said "H (a, b) and C_f (c, d)". The plan uses the caption layout (four panels).
- F11:
  - [S] tag "rework"; caption "Complete."
  - [F] Source: `layer_state` L5, L6; extended to the 20M LES and 193M.
- F12:
  - [S] tag "rework"; caption "Complete."
  - [F] Source: `validation_50M` V1; 100M and 193M from `validation_dns`; centreline and mid-pitch panels new.
  - [N] Experiment dataset id727 (DR 1.5, BR 1.5, Tu 0.5 %); plan says "experiment".
- F13: [S] tag "new"; caption "Complete."
- (iv) [D] "Owner decision (2026-11-25b): no effective diffusivity alpha_eff and no tensor form of Pr_t anywhere in the
  paper; one linear story." Plan keeps the rule without the attribution.
- F14:
  - [S] tag "rework"; caption "Complete."; "candidate fig_prt_classic.png".

## 4.2

- [D] Purpose: "(owner, 2026-11-25b)"; "(owner, 2026-11-26: neutral voice, no advocacy)". Wording changed to
  "two terms of opposite sign" and "higher than the reference value".
- [N] (i) "(audit, planning/maths_audit/)" [F] after the sentence on reporting L_al, L_eff and cos(theta).
- [P] (i) "The 193M volume check waits for fc's epsilon file (27 Nov)."
- [D] (ii) "Order: ... only (owner, 2026-11-26)".
- [P] (ii) Scale-resolving bound: "(W1 item 3, 1 Dec)".
- [F] Step 1: "`maths_audit/estimator`, mode C"; [P] "time-batch SE from fc, 27 Nov".
- [H] "Superseded (2026-11-26): the thresholded mean-of-ln values of W1 item 1 (jet +0.57, film -0.50, -0.70); the jet
  sign depended on that estimator. fc's W1 table stays in FINDINGS C25."
- [P] "Volume values of the classic Pr_t,R from W1 item 4 (27 Nov) replace the plane values." The plan states the plane
  values.
- [P] Step 3: "W2 (2 Dec) adds the 100M and the volume regions."
- [N] RANS sweep bullet: "(STAR session, 2026-11-25; FINDINGS LP6, reference 193M b010)".
- [P] "4.3 adds two runs approved on 2026-11-25: P1 ..., P2 ...; both 50M SST, start by 28 Nov." "RC-S9 (approved
  2026-11-20) ...; 50M production ends about 26 Nov."
- F15: unchanged; no items.
- F16:
  - [S] tag "new"; caption "Complete; volume values replace the plane values when W1 item 4 lands."
  - [P] Source: Pr_t analysis; volume values from W1 item 4 (27 Nov).
  - [U] F16 caption near-film 95 % interval is 0.56-0.74; the 4.2 text gives 0.57-0.73 for the same quantity.
    Both left as written.
- F17: [S] tag "new"; caption "Complete."; format line moved into HOW.
- F18: [S] tag "rework"; caption "Complete."; [F] Source: `meff_*` F2; "Caption states both medians."
- F19: [S] tag "new"; caption "Complete."; [F] Source: `dns/python/tke_budget_pf.py` and `budgets_c50b` B1, extended to the 193M.
- F20:
  - [S] tag "rework", "simpler"; caption "Complete."
  - [P] "x/D=4 not read (19.9 GB per batch); W2 (2 Dec) adds the 100M and volume regions."
  - [V] "text must say so" (mu_t term compares with the classic wall-normal viscosity) kept as a plain statement.

## 4.3

- [P] "Planned once the tests are done (P1 zonal Pr_t, P2 film-only substitution, RC-S9 gamma=1; with the LP6 sweep)."
- [O] "Open question for 4.3 (2026-11-26): QCR recovers 72 % of the x/D 7-20 eta-bar deficit against 36-48 % for a
  Pr_t setting alone, and no run changes the length scale; how 4.3 tests the step-1 finding is to be decided after
  the split audit."

## 5 Conclusions

- [P] "Planned after 4.3." Plan: one line on contents.

## Unsure (content or history; left in the plan)

1. F4 caption: the plan describes the final F4 (eta_K,193 for every mesh, five meshes). The drawn file is the
   interim version (own eta_K per run). Caption lists batches b006/b010 as before. Check when F4 is redrawn.
2. F10 panel layout (see above).
3. F16 near-film 95 % interval, 0.56-0.74 (caption) against 0.57-0.73 (4.2 text).
4. Batch labels b006/b010 stay in captions and text as data definitions (reference = b010).
5. "Not seen before: the rung-to-rung BR spread ..." (4.1 (i)) and "new here: ..." (4.2) are novelty claims, kept as content.
6. 4.2 sweep bullet still lists the P1, P2 and RC-S9 runs and their numbers (5M +1.3 %, 50M smoke +1.26 %) as
   content of 4.2; the approval and schedule text was removed.
7. 4.1 (i) "Effect on the mainstream (50M analysis)": kept as a result; the 193M repeat is listed above under pending.
8. F3 WHAT lists the class numbers 1, 9, 11, 2, 7, 8, 10 (from the 3.2 text) and F4/F10/F12/F13 WHAT lines reuse
   numbers already in the plan; no new numbers were added.

## Tally of items moved (approximate, counting each bullet above)

- S status tags and status columns: 4 list-level items (figure, table, equation, legend) plus 20 figure-entry
  tags/"Complete." captions and 3 draft/interim caption prefixes: 27.
- D dates and owner decisions: 6.
- P pending/due/to-finalise: 14.
- V provenance and source bullets: 7.
- N agent, session, commit, finding, id names: 4.
- F file paths of analysis data and Source bullets: 20.
- H history/superseded: 5.
- O open questions: 2.
- U unsure items listed: 8.
