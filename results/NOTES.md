# Notes on result files

## `table8_oneway_reproduced_*` and `table8_oneway_with_northmac_*`

These are the core outputs for the strict 1-source analysis on the retained
public-HO target panel:

- `Croatian`
- `Serbian_Serb`
- `Romanian`
- `Bulgarian`
- `Albanian`
- `Greek`
- `Cypriot`

`table8_oneway_reproduced_*` contains Olalde's six published substrate columns
on this strict AADR-label panel. `table8_oneway_with_northmac_*` appends
`NorthMacedonia_IA` as a seventh source column under the same setup.
`table8_oneway_summary.tsv` reports, for each retained target label, the best
p-value among the six published columns and the p-value for the added
`NorthMacedonia_IA` column.

## `table8_oneway_with_northmac_allsnps_false_*`

These files repeat the same strict 1-source screen under qpAdm's default
SNP-intersection mode, equivalent to `allsnps=FALSE`. The pass/fail pattern is
unchanged relative to the main `allsnps=TRUE` output: every published substrate
column fails every retained target row, while `NorthMacedonia_IA` passes all
seven retained rows. The p-values do change, with the closest
`NorthMacedonia_IA` pass being `Serbian_Serb` (`p = 0.05019`).

`table8_oneway_allsnps_true_false_comparison.tsv` joins the main summary with
the `allsnps=FALSE` summary for quick inspection.

## `T9_allsnps_sensitivity.tsv`

This file is treated as an upstream `allsnps=TRUE` cache by `scripts/03_qpadm_run.R`.
For 1-source rows, ignore `z_w_summary`: the fitted weight is fixed to 1 by
construction, so the reported z-scores are numerical artifacts of dividing by a
near-zero standard error. Only the `pvalue` field is meaningful for the current
1-source pipeline.
