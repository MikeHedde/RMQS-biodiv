# Theoretical model of taxonomic uncertainty

## Purpose

This is a stand-alone, generative simulation framework for studying how the
architecture of a regional taxonomic pool, community structure, identification
error and reporting workflow jointly create apparent biodiversity change.

It is deliberately separate from the empirical Collembola workflow.

The empirical workflow answers:

> What happens when documented error and workflow scenarios are applied to an observed dataset?

This theoretical module answers:

> Across a broad, explicit parameter space, what types and amplitudes of artefact should be expected from different taxonomic architectures and observation processes?

The same framework can be used for Collembola, earthworms, ants or another taxon.

## Core model

The model is:

```
true community N
      |
      | identification process K
      v
error-affected species matrix
      |
      | reporting / workflow process C
      v
observed matrix Z = N K C
```

The simulation has four independent parameter blocks.

1. **Regional taxonomic architecture**
   - number of regional species, genera and families;
   - unevenness of species among genera;
   - unevenness of genera among families;
   - fraction of difficult genera.

2. **True community**
   - number of sites and observed species;
   - occupancy distribution;
   - abundance distribution;
   - abundance–occupancy relationship;
   - optional environmental-gradient structure.

3. **Identification process**
   - abundance-weighted mean error rate;
   - candidate pool restricted to observed congeneric species or open to regional congeneric species;
   - rare-species bias;
   - difficult-genus bias.

4. **Reporting / workflow**
   - species-level reporting;
   - mixed species + genus RTUs;
   - removal of unresolved records;
   - rare species reported as genus;
   - difficult genera reported as genus;
   - full genus or family aggregation.

## What it is not

This is not a calibration model for Collembola. The simulation engine should first
be explored over a theoretical parameter space. Empirical taxa are then overlaid
on that space using descriptive properties such as regional richness, taxonomic
architecture, number of sampled sites and observed richness.

The model becomes externally predictive only after it successfully anticipates
the responses of taxa, regions or datasets not used to define the theoretical
experiment.

## Files

- `00_theory_config.R`: isolated dependencies and paths.
- `01_taxonomic_architecture.R`: regional species–genus–family hierarchy.
- `02_community_generator.R`: true site × species matrix.
- `03_observation_process.R`: identification error and workflow.
- `04_metrics.R`: alpha, gamma, occupancy and beta diversity.
- `05_simulation_engine.R`: scenario-pair and replication engine.
- `06_experiment_design.R`: random theoretical world design.
- `07_empirical_bridge.R`: generic import/overlay functions only.
- `08_plotting.R`: theoretical figures.
- `000_run_theory_demo.R`: a self-contained demonstration.

## First run

Open `000_run_theory_demo.R` and source it.

The demo creates 60 theoretical worlds and five stochastic replicates per world.
Increase these values only after checking runtime.

## Empirical comparison

The bridge is intentionally one-way.

1. Run the empirical workflow as usual.
2. Read its output table `scenario_summary_by_iter.csv`.
3. Use `empirical_blowes_from_summary()` to compute empirical points.
4. Use `empirical_taxon_profile()` to describe the taxon.
5. Select comparable theoretical worlds with `select_theoretical_worlds()`.
6. Overlay the empirical scenario points on the theoretical cloud.

No empirical result is used to tune the simulation's response variables.

## Recommended next scientific step

Before applying the module to any taxon, define a small, defensible set of
theoretical contrasts:

- low versus high species-per-genus concentration;
- small versus large regional confusion pool;
- low versus high occupancy heterogeneity;
- observed-pool versus regional-pool identification errors;
- genus versus family coarsening.

This makes the first paper about mechanisms rather than about a large
unstructured simulation grid.

## Second-stage mechanism screen

After the first theoretical simulation, run `001_explain_theory_demo.R`.

This produces an exploratory parameter screen for the two core mechanisms:

- gamma inflation and apparent differentiation under regional-pool
  congeneric error;
- increased mean occupancy and apparent homogenisation under family-level
  aggregation.

The screen uses rank associations across theoretical worlds and response-curve
figures. Its role is to identify the small number of mechanisms to hold fixed
or vary systematically in the next, designed factorial simulation; it is not a
calibration or a final causal model.

## Controlled one-factor experiment

The exploratory random-world screen identifies candidate mechanisms but mixes
many parameter changes. The next step is a controlled one-factor-at-a-time
(OFAT) experiment implemented in:

- `10_controlled_one_factor_experiment.R`
- `11_controlled_one_factor_plots.R`
- `002_run_controlled_one_factor.R`

This experiment varies one focal axis at a time, holding the other three at a
reference value:

1. regional / observed species ratio;
2. mean species per genus;
3. mean genera per family;
4. initial mean occupancy.

Mean genera per family is used instead of mean species per family because it
is independent of mean species per genus. Mean species per family is therefore
an emergent property of the two taxonomic-depth axes.

The controlled run crosses these axes with:

- five nominal identification-error rates (1, 3, 5, 10 and 20% by default);
- observed-pool, regional-pool, rare-weighted and difficult-genus-biased
  congeneric errors;
- mixed RTU reporting and dropping unresolved records at the same rate range;
- fixed workflows: reporting rare taxa at genus, reporting difficult genera at
  genus, full genus aggregation and full family aggregation.

Run `002_run_controlled_one_factor.R` after the basic module has been tested.
The first version uses five stochastic replicates per design point. Increase
`N_REPLICATES_PER_WORLD` to 20 or more for final figures.

This is intentionally an OFAT design. It estimates direct mechanism-specific
response curves, not interactions among the four axes. A subsequent, smaller
factorial interaction experiment should be run only after these direct effects
identify the parameters worth crossing.


## Cross-taxon empirical profiles and taxon-specific theory

The controlled OFAT experiment establishes general mechanisms. The next stage
is to give each empirical group its own **descriptive profile**, without using
empirical uncertainty outcomes as calibration targets.

New files:

- `12_taxon_profile_builder.R`: reads a standard site/species/abundance file and
  a standard regional taxonomy; calculates the four central axes and the
  effective congeneric confusion capacity.
- `13_taxon_specific_simulation.R`: applies the shared theoretical scenarios to
  synthetic communities built on each group’s *real regional hierarchy*.
- `14_taxon_profile_plots.R`: profile and cross-taxon comparison figures.
- `003_run_taxon_profiles.R`: profile builder and optional taxon-specific
  simulation launcher.
- `bridge_from_empirical_workflow/export_collembola_profile_inputs.R`: optional
  one-way exporter for the current Collembola workflow.
- `inputs/taxon_profile_inputs.R`: edit this configuration to add Earthworms
  and Ants once their standard inputs exist.

### Standard input format

Each group needs two CSV files:

```text
community_long.csv: site, species, abundance
regional_taxonomy.csv: species, genus, family[, difficult_genus]
```

`species` must use the same identifier in both files. The taxonomy should be the
regional pool appropriate to the question (e.g. mainland France), not only the
species observed in the monitoring network.

### Recommended use

1. In the empirical Collembola project, source the one-way exporter after steps
   01–02. Copy its two CSV files to this module's `inputs/` directory.
2. Add equivalent standard files for Earthworms and Ants.
3. Edit `inputs/taxon_profile_inputs.R` to activate the available groups.
4. Run `003_run_taxon_profiles.R`.

The generated cross-taxon figures compare **expected artefact susceptibility**
conditional on real taxonomy and sampling descriptors. They are not a fit to,
or validation against, the empirical scenario responses. Validation comes later
by comparing those simulated expectations with the independently computed
empirical scenario analyses for each group.
