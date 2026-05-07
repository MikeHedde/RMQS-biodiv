library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

gr <- grViz("
digraph dag {
  graph [rankdir = LR, bgcolor='white', pad=0.15, nodesep=0.35, ranksep=0.55]
  node  [shape=box, style='rounded,filled', fillcolor='white', color='gray25',
         fontname='Helvetica', fontsize=10]
  edge  [color='gray35', arrowsize=0.7]

  # --- Sampling design (constant across years)
  DESIGN [label='Site-selection design\\n(RMQS grid, land-use representativeness)']

  # --- Years as independent realizations
  subgraph cluster_years {
    label='Two independent realizations';
    style='rounded,dashed'; color='gray70';
    Y2024 [label='2024 realization\\n(core protocols, fewer reps)']
    Y2025 [label='2025 realization\\n(extended / modular protocols)']
  }

  DESIGN -> Y2024
  DESIGN -> Y2025

  # --- Drivers and conditions (vary by realization)
  subgraph cluster_drivers {
    label='Ecological drivers & sampling conditions';
    style='rounded,dashed'; color='gray70';
    ENV  [label='ENV\\n(land use, soil, climate, vegetation)']
    MET  [label='Weather / phenology\\n(during sampling window)']
  }
  Y2024 -> ENV
  Y2025 -> ENV
  Y2024 -> MET
  Y2025 -> MET

  # --- Latent compartments
  subgraph cluster_latent {
    label='Latent compartments (true state)';
    style='rounded,dashed'; color='gray70';
    LSOIL [label='Soil community\\n(latent)']
    LLIT  [label='Litter community\\n(latent)']
    LSURF [label='Surface-active community\\n(latent)']
    LVEG  [label='Vegetation-layer community\\n(latent)']
  }
  ENV -> LSOIL; ENV -> LLIT; ENV -> LSURF; ENV -> LVEG
  MET -> LSURF; MET -> LVEG  # activity-sensitive

  # --- Observation filter (explicit chain, compact)
  subgraph cluster_obsfilter {
    label='Observation filter (method-specific)';
    style='rounded,dashed'; color='gray70';
    DET [label='Detectability components\\np_access × p_capture × p_ID']
    EFF [label='Effort\\n(reps, duration, volume/area)']
  }

  # --- Observed datasets (protocol modules)
  subgraph cluster_obs {
    label='Observed data (Y): protocol modules';
    style='rounded,dashed'; color='gray70';
    HS  [label='Hand-sorting std\\n(soil block)']
    HQ  [label='Hand-sorting qual\\n(list completion)']
    CORE[label='Soil cores + Macfadyen']
    LIT [label='Litter + Macfadyen']
    PIT [label='Pitfall traps\\n(Barber / GPD / etc.)']
    DVAC[label='DVAC (open habitats)']
  }

  # Latent -> Protocol modules (dominant links)
  LSOIL -> HS
  LSOIL -> CORE
  LLIT  -> LIT
  LSURF -> PIT
  LVEG  -> DVAC
  LSOIL -> HQ

  # Filters -> Observations
  DET -> HS; DET -> HQ; DET -> CORE; DET -> LIT; DET -> PIT; DET -> DVAC
  EFF -> HS; EFF -> HQ; EFF -> CORE; EFF -> LIT; EFF -> PIT; EFF -> DVAC

  # --- Indicators and downstream inference
  subgraph cluster_out {
    label='Inference products';
    style='rounded,dashed'; color='gray70';
    COV  [label='Coverage / completeness']
    ALPH [label='Alpha diversity\\n(Hill q0/q1/q2)']
    BETA [label='Community structure\\n(ordination / dissimilarity)']
    SDM  [label='SDM / jSDM\\n(species–environment)']
    NET  [label='Co-occurrence /\\ninteraction networks']
  }

  HS -> COV; CORE -> COV; LIT -> COV; PIT -> COV; DVAC -> COV; HQ -> COV
  COV -> ALPH
  HS -> BETA; CORE -> BETA; LIT -> BETA; PIT -> BETA; DVAC -> BETA
  BETA -> NET
  ALPH -> SDM
  BETA -> SDM

  # --- Year-specific availability of modules (annotation only)
  edge [style=dotted, color='gray55']
  Y2024 -> HS; Y2024 -> CORE; Y2024 -> PIT
  Y2025 -> HS; Y2025 -> HQ; Y2025 -> CORE; Y2025 -> LIT; Y2025 -> PIT; Y2025 -> DVAC
  edge [style=solid, color='gray35']
}
")

gr
svg_txt <- export_svg(gr)
writeLines(svg_txt, "figures/MEE_protocols/1.DAG/DAG_RMQS_multimethod_updated.svg")
rsvg_pdf(charToRaw(svg_txt), file = "figures/MEE_protocols/1.DAG/DAG_RMQS_multimethod_updated.pdf")
