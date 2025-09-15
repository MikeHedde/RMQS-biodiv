
# app.R — RMQS Araneae Explorer (Shiny)
# Public-ready: reads CSV exported from RMQS (update_araneae.csv), maps columns,
# and serves interactive filters + map + tables.
#
# How to run locally:
   install.packages(c("shiny","shinydashboard","DT","leaflet","dplyr",
                      "tidyr","stringr","readr","sf","plotly","purrr"))
   shiny::runApp("Shiny/app")
#
# Data source:
# - Local: "data/update_araneae.csv"
# - Or set env var RMQS_CSV_URL to a raw GitHub URL to read remotely.

packages <- c("shiny","shinydashboard","DT","leaflet","dplyr","tidyr","stringr","readr","sf","plotly","purrr")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

RMQS_LOCAL <- Sys.getenv("RMQS_LOCAL", unset = "data/raw-data/lab files/databases update/update_araneae.csv")
RMQS_URL   <- Sys.getenv("RMQS_CSV_URL", unset = "")

REQUIRED_COLS <- c("site_id","date","latitude","longitude","department","order","family","genus","species","abundance")

# Read and map the RMQS CSV (semicolon, Latin1 typical for French Windows exports)
# Lecture CSV (;, Latin1) et mapping colonnes -> schéma commun


filtered <- reactive({
df <- apply_filters(raw_data, input$f_order, input$f_family, input$f_species, input$f_dep)
if (isTRUE(input$only_occ)) df <- df |> dplyr::filter(is.na(abundance) | abundance > 0)
df
})


output$tbl <- DT::renderDT({
df <- filtered() |>
dplyr::mutate(date = as.character(date)) |>
dplyr::select(site_id, date, department, order, family, species_full, abundance, latitude, longitude)
DT::datatable(df, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE), filter = "top")
})


output$dl_filtered <- downloadHandler(
filename = function() sprintf("rmqs_araneae_filtre_%s.csv", Sys.Date()),
content = function(file) readr::write_csv(filtered(), file)
)


# Carte
mapdata <- reactive({ apply_filters(raw_data, input$m_order, input$m_family, input$m_species, input$m_dep) })
output$map <- leaflet::renderLeaflet({
df <- mapdata()
pal <- colorFactor("Set2", df$order)
leaflet(df) |>
addProviderTiles(providers$CartoDB.Positron) |>
addCircleMarkers(~longitude, ~latitude, color = ~pal(order), radius = 6, stroke = FALSE, fillOpacity = 0.85,
clusterOptions = if (isTRUE(input$cluster)) markerClusterOptions() else NULL,
label = ~lapply(paste0("<b>", species_full, "</b><br/>", order, " · ", family,
"<br/><i>Site</i>: ", site_id, "<br/><i>Dep</i>: ", department), htmltools::HTML)) |>
addLegend("bottomright", pal = pal, values = ~order, title = "Ordre")
})


# Graphes synthèse
output$plt_top_species <- plotly::renderPlotly({
top <- featured_species(raw_data, 8)
p <- ggplot2::ggplot(top, ggplot2::aes(x = reorder(species_full, n_sites), y = n_sites, fill = family)) +
ggplot2::geom_col() + ggplot2::coord_flip() +
ggplot2::labs(x = NULL, y = "Nombre de sites") + ggplot2::guides(fill = ggplot2::guide_legend(title = "Famille"))
plotly::ggplotly(p, tooltip = c("x","y","fill"))
})


output$plt_rich_dep <- plotly::renderPlotly({
df <- raw_data |>
dplyr::filter(is.na(abundance) | abundance > 0) |>
dplyr::distinct(department, site_id, order, family, species_full) |>
dplyr::group_by(department) |>
dplyr::summarise(richness = dplyr::n_distinct(species_full), .groups = "drop")
p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(department, richness), y = richness)) +
ggplot2::geom_col() + ggplot2::coord_flip() +
ggplot2::labs(x = "Département", y = "Richesse spécifique (occurrences uniques)")
plotly::ggplotly(p, tooltip = c("x","y"))
})


output$tbl_top <- DT::renderDT({
dt <- raw_data |>
dplyr::group_by(order, family, genus, species_full) |>
dplyr::summarise(n_sites = dplyr::n_distinct(site_id), n_records = dplyr::n(), .groups = "drop") |>
dplyr::arrange(dplyr::desc(n_sites))
DT::datatable(dt, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
})
}


shinyApp(ui, server)