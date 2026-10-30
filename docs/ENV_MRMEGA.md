# env-MR-MEGA with the study protocol as the environmental covariate

Decision (30 Oct 2026): env-MR-MEGA is on the path, not deferred. The environmental variable is the
**study protocol** — the treatment protocol each participant was enrolled on. This file states the model,
the unit of meta-regression, how the covariate is built from the handoff bundle, and what the fixture test
has to show. The implementation is Track B task CP-29 (bash + SLURM, in the `mrmega` conda env).

## The model

MR-MEGA (Mägi et al. 2017) meta-regresses each variant's per-study effect on the first `pc_axes` axes of
genetic variation derived from study allele frequencies: `beta_k = a + sum_j b_j * PC_jk + e_k`, weighted
by `1 / se_k^2`. env-MR-MEGA (Wang et al., Communications Biology 2024, PMID 39543362) adds
**study-level** environmental covariates to the same regression: `beta_k = a + sum_j b_j * PC_jk + sum_m
c_m * E_mk + e_k`. Three tests come out per variant: the association (all coefficients), ancestry
heterogeneity (the `b_j`), and environmental heterogeneity (the `c_m`). The covariates are one value per
study, never per individual, so the method runs on the same per-study summary statistics MR-MEGA uses.

## The unit, and the covariate

The unit (a "study") is an ancestry stratum from `ancestry/strata.tsv`, exactly as for MR-MEGA. The
protocol covariate is therefore a **stratum-level protocol mix**: for a protocol variable with levels
`p_1 … p_P`, the covariates are the proportions of each stratum's samples on `p_2 … p_P` (one level is
the reference, so the proportions are not collinear with the intercept). They are computed from
`phenotypes/phenotypes.tsv` column `protocol` restricted to the samples that entered that stratum's
GWAS for that trait (the complete-case set), so the covariate describes the people whose effect it is
attached to. The table is written per trait as `results/mrmega_env/<trait>.study_covariates.tsv`
(`study  pc1..pc_k  pct_<level>…  n`).

Sensitivity unit (reported, not primary): **stratum × protocol** as the study, where every cell reaches
the stratum floor (50); cells below the floor are not run. This is the "Option B" of the earlier design
notes and it is what makes the protocol a categorical study-level factor rather than a proportion.

## What upstream must ship

`phenotypes.tsv` carries a required coded column `protocol` (a short code per protocol, never free text;
`NA` allowed and treated as missing). The handoff checker refuses a phenotype file without it. Levels and
their meaning are listed in `docs/handoff/protocol_levels.tsv` once the REDIAL protocol set is known.

## What the fixture test shows (CP-29, on the synthetic set)

- `M-01` simulates a `protocol` column with three levels whose mix differs across the simulated strata.
- The covariate builder reproduces the simulated per-stratum proportions exactly (golden file).
- env-MR-MEGA runs on all three traits and writes the three p-values per variant; on the null rung of the
  simulation the environmental-heterogeneity p-values are uniform (QQ against the null, lambda within
  0.95–1.05); on a rung where one causal variant's effect is simulated to scale with the protocol mix, that
  variant is the top environmental-heterogeneity hit.
- The run record names the protocol levels, the reference level, and the per-stratum proportions used.
