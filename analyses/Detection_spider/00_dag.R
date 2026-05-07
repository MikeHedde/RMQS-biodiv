library(rsvg)
library(DiagrammeRsvg)

dag_plot <- grViz("
digraph DAG {
  graph [rankdir=LR, splines=true, overlap=false, outputorder=edgesfirst]
  node  [shape=ellipse, fontname=Helvetica, fontsize=12]
  edge  [fontname=Helvetica, fontsize=10, penwidth=1, arrowsize=0.8]

  # -------------------------
  # Nodes
  # -------------------------
  M   [label='Protocol (M)']
  E   [label='Effort (E)']
  p   [label='Detectability (p̂)']
  Y   [label='Observed data (Y)']
  C   [label='Coverage (C)']
  D   [label='Diversity (D)\\n(Hill q=0,1,2 at C=0.95)']

  Env [label='Season / site covariates\\n(DOY, habitat, etc.)']
  Z   [label='Occurrence (Z)']
  N   [label='Local abundance (N)']
  A   [label='Activity / mobility (A)']
  T   [label='Traits (T)\\n(size, wings)']

  # -------------------------
  # Highlighted nodes (studied causal chain)
  # -------------------------
  M [fontsize=16, style=bold]
  E [fontsize=16, style=bold]
  p [fontsize=16, style=bold]
  C [fontsize=16, style=bold]
  D [fontsize=16, style=bold]

  # -------------------------
  # Full DAG edges (context)
  # -------------------------
  Env -> Z
  Env -> N
  Z   -> N
  Env -> A
  Z   -> A
  N   -> A

  T   -> p
  A   -> p
  Env -> p

  Z   -> Y
  N   -> Y
  p   -> Y
  M   -> Y   [label='device window']

  Y   -> C

  # -------------------------
  # Core studied pathway (thick + bold edges)
  # -------------------------
  M -> p [penwidth=3.2, style=bold]
  E -> p [penwidth=3.2, style=bold]
  p -> Y [penwidth=3.2, style=bold]
  Y -> C [penwidth=3.2, style=bold]
  C -> D [penwidth=3.2, style=bold]

  # -------------------------
  # Layout hints (optional but improves readability)
  # -------------------------
  {rank=same; Env; T; M; E}
  {rank=same; Z; N; A}
  {rank=same; p}
  {rank=same; Y}
  {rank=same; C}
  {rank=same; D}
}
")

svg_txt <- export_svg(dag_plot)
rsvg_png(
  charToRaw(svg_txt),
  file = "figures/detectability/carabidae/carabidae_dag.png",
  width = 2400
)
