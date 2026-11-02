args <- commandArgs(TRUE)
skill <- file.path(Sys.getenv("HOME"),".agents/skills/jins-routine-plotting")
for (f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Stacked_composition/plot.R"))
  source(file.path(skill,f))
jin_setup()
dt <- fread(file.path(args[1],"tables/overview_region_lineages.tsv"))
dt <- dt[source_context=="tumor_study_specimen_not_region_diagnosis" & axis=="lesion_location" & region %in% c("inside","outside")]
pal <- fread(file.path(Sys.getenv("MP10_PROJECT",getwd()),"neighborhood_niche_visualization_v1/data_tables/lineage_palette.tsv"))
# Explicit project mapping, retaining the 9 reference major lineages.
cols <- setNames(pal$color,pal$lineage)
colors <- c(Epithelial=unname(cols["Epithelial"]),Fibroblast=unname(cols["CAF"]),Myeloid=unname(cols["Myeloid"]),
            Endothelial=unname(cols["EC"]),Mural=unname(cols["Mural"]),B=unname(cols["B/Plasma"]),Plasma="#A38EBC",
            T_NK=unname(cols["T/NK"]),Mast=unname(cols["Mast"]))
p <- plot_stacked_composition(dt,"region","feature","fraction",condition_order=c("inside","outside"),category_order=names(colors),
    colors=colors,xlab="v3.4 candidate location",ylab="Mean reference fraction",legend_title="Lineage")
jin_export(p,file.path(args[1],"figures","summary","lineage_composition"),5.4,4.7,legend_plot=p,panel_width=4.1)

