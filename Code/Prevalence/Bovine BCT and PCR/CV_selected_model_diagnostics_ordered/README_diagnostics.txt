AAT selected-covariate + spatial + IID validation: two figures
Generated: 2026-10-24 22:06:59.233385
Draw IDs: 1
Detailed point-plot draw: 1
Missing analyses: none
Coverage_pooled.csv uses the ratio of all covered evaluations to all evaluations, not the unweighted mean of fold percentages.
The primary aggregation gives each eligible cell equal weight and each surveyed site equal weight within a cell.
Only selected_spatial_iid and equal_sites evaluations are included; all comparator models and sample-size-weighted results are excluded.
The script reads saved predictions and performs no fitting or new posterior/predictive simulation.
All intervals scored against observations are posterior predictive intervals (new-site IID plus binomial sampling).
expected_prevalence_lo95/hi95 and p_lo95/hi95 are deliberately NOT scored against observed counts/fractions.
Coverage is conditional on the included diagnostic-correction draws. Repeated draws are not independent replicates.
No binomial/Wilson confidence intervals are assigned to coverage: spatial and shared-fit dependence would invalidate a naive independent-trials interpretation.
Fold/draw points describe heterogeneity; they are not inferential confidence bounds.
The 100-km scripts retain 200-km withholding blocks. Original non-aggregated scripts typically use 100-km withholding blocks.
Thus the three primary targets do not automatically share the same CV design; consult Fold_geometry.csv and the saved assignments.
Grid-cell means validate sampled populations, not unbiased areal mean prevalence. Finer cells may exclude more sparse areas.
80% and 95% are the only common saved interval levels; extra levels cannot be recovered from endpoint summaries.
Coverage and errors do not by themselves establish accurate spatial discrimination or calibrate mapped expected-prevalence credible intervals.
Coverage pools the included draws; the predictive plot uses representative_draw, with metrics computed only from its displayed points.
Source intervals, observations and coverage flags are exported in Interval_coverage_audit.csv.
Predictive_plot_data.csv records the within-panel rank and outside95 flag for each displayed unit.
The ordered plot joins pointwise intervals for display; it is not a simultaneous band or a spatial/temporal trajectory.
Vector PDFs are suitable for layout/editing; PNGs are 600 dpi by default.
Input_inventory.csv records the exact source paths and available/included draws.
