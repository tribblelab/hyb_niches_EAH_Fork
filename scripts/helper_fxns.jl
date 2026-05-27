using CSV, DataFrames
using RCall

R"""
library(gatoRs)
library(ggplot2)
library(sf)
library(ggspatial)
library(gridExtra)
library(CoordinateCleaner)
library(readxl)
library(dplyr) #needs to be loaded in last, so it defaults to correct `filter` fxn
"""

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_taxa_data(traits_path) -> (taxa_traits, taxa_to_nativerange_dict)

Load taxa traits CSV and build a scientificName → botanicalCountries lookup.
"""
function load_taxa_data(traits_path::String)
    taxa_traits = DataFrame(CSV.File(traits_path))
    taxa_to_nativerange_dict = Dict(
        row.scientificName => row.botanicalCountries for row in eachrow(taxa_traits)
    )
    return taxa_traits, taxa_to_nativerange_dict
end

"""
    load_bot_regions(shapefile_path)

Load TDWG botanical country polygons into R (once, globally).
"""
function load_bot_regions(shapefile_path::String)
    @rput shapefile_path
    R"bot_regions = st_read(shapefile_path, quiet = TRUE)"
end

"""
    load_georef_df(taxa_filename, georef_dir) -> DataFrame or nothing

Look for a georeferenced points file (.csv or .xlsx) for a taxon, filter to
viable rows, and return a DataFrame with Float64 lat/lon columns — or `nothing`
if no usable file is found.
"""
function load_georef_df(taxa_filename::String, georef_dir::String)
    for ext in [".csv", ".xlsx"]
        candidates = filter(
            f -> startswith(f, taxa_filename) && endswith(f, ext),
            readdir(georef_dir)
        )
        isempty(candidates) && continue

        georef_path = joinpath(georef_dir, first(candidates))
        println("  Found georeferenced file: $georef_path")

        try
            df = if ext == ".csv"
                DataFrame(CSV.File(georef_path))
            else
                @rput georef_path
                R"georef_r <- readxl::read_excel(georef_path)"
                @rget georef_r
                georef_r
            end

            viable_col = findfirst(c -> uppercase(string(c)) == "VIABLE", names(df))
            if isnothing(viable_col)
                println("  WARNING: No 'Viable' column found in $georef_path — skipping georef points")
                return nothing
            end

            rename!(df, names(df)[viable_col] => :Viable)
            df[!, :Viable] = map(v -> uppercase(strip(string(v))), df.Viable)
            filter!(row -> row.Viable == "TRUE", df)

            df[!, :latitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.latitude)
            df[!, :longitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.longitude)
            filter!(row -> !ismissing(row.latitude) && !ismissing(row.longitude), df)

            println("  Georeferenced viable points: $(nrow(df))")
            return nrow(df) > 0 ? df : nothing

        catch e
            println("  WARNING: Could not read georef file $georef_path: $e")
            return nothing
        end
    end
    return nothing
end

"""
    resolve_clean_file(filename, clean_dir) -> String

Return the best available cleaned CSV path for a taxon in `clean_dir`:
prefers `*_georef_merged.csv`, falls back to any `*_cleaned.csv`.
Returns an empty string if nothing is found.
"""
function resolve_clean_file(filename::String, clean_dir::String)
    taxon_stem = replace(filename, r"-\d{4}_\d{2}.*$" => "")
    georef_candidates = filter(
        f -> startswith(f, taxon_stem) && endswith(f, "_georef_merged.csv"),
        readdir(clean_dir)
    )
    cleaned_candidates = filter(
        f -> startswith(f, taxon_stem) && occursin("cleaned", f) && endswith(f, ".csv"),
        readdir(clean_dir)
    )
    if !isempty(georef_candidates)
        return joinpath(clean_dir, first(georef_candidates))
    elseif !isempty(cleaned_candidates)
        return joinpath(clean_dir, first(cleaned_candidates))
    else
        return ""
    end
end

"""
    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)

Transfer all data needed by the R plotting block to R's environment.
"""
function push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)
    @rput taxa raw_file clean_file country_codes
    has_georef = !isnothing(georef_df) && nrow(georef_df) > 0
    @rput has_georef
    if has_georef
        georef_df_plot = georef_df[!, [:latitude, :longitude]]
        @rput georef_df_plot
    end
end

# Core R plotting block (rendered into a grid.arrange object in R)

const RPLOT_BLOCK = """
    raw_df  <- read.csv(raw_file)
    clean_df <- read.csv(clean_file)
    raw_df  <- raw_df[!is.na(raw_df\$latitude) & !is.na(raw_df\$longitude), ]

    if (nrow(raw_df) == 0 || nrow(clean_df) == 0) stop("empty data")

    raw_lon_min <- min(raw_df\$longitude, na.rm=TRUE) - 3
    raw_lon_max <- max(raw_df\$longitude, na.rm=TRUE) + 3
    raw_lat_min <- min(raw_df\$latitude,  na.rm=TRUE) - 3
    raw_lat_max <- max(raw_df\$latitude,  na.rm=TRUE) + 3

    all_clean_lon <- clean_df\$longitude
    all_clean_lat <- clean_df\$latitude
    if (has_georef) {
        georef_df_plot\$longitude <- as.numeric(georef_df_plot\$longitude)
        georef_df_plot\$latitude  <- as.numeric(georef_df_plot\$latitude)
        georef_df_plot <- georef_df_plot[!is.na(georef_df_plot\$latitude) &
                                         !is.na(georef_df_plot\$longitude), ]
        all_clean_lon <- c(all_clean_lon, georef_df_plot\$longitude)
        all_clean_lat <- c(all_clean_lat, georef_df_plot\$latitude)
    }
    clean_lon_min <- min(all_clean_lon, na.rm=TRUE) - 2
    clean_lon_max <- max(all_clean_lon, na.rm=TRUE) + 2
    clean_lat_min <- min(all_clean_lat, na.rm=TRUE) - 2
    clean_lat_max <- max(all_clean_lat, na.rm=TRUE) + 2

    world    <- annotation_borders(database="world", colour="gray80", fill="gray80")
    countries <- annotation_borders(database="world", colour="gray40", fill=NA, size=0.5)

    bot_underlay <- bot_regions[bot_regions\$LEVEL3_COD %in% country_codes, ]
    underlay_layer <- if (nrow(bot_underlay) > 0) {
        geom_sf(data=bot_underlay, fill="khaki1", color="goldenrod4",
                alpha=0.25, linewidth=0.3)
    } else { NULL }

    p1 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=raw_df, aes(x=longitude, y=latitude),
                   color="blue", size=1.5, alpha=0.6) +
        coord_sf(xlim=c(raw_lon_min, raw_lon_max),
                 ylim=c(raw_lat_min, raw_lat_max)) +
        labs(title=paste0("Raw (n=", nrow(raw_df), ")"),
             x="Longitude", y="Latitude") +
        theme_minimal() +
        theme(plot.title=element_text(size=10))

    clean_title <- if (has_georef) {
        paste0("Cleaned (n=", nrow(clean_df), " + ", nrow(georef_df_plot), " georef)")
    } else {
        paste0("Cleaned (n=", nrow(clean_df), ")")
    }

    p2 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=clean_df, aes(x=longitude, y=latitude),
                   color="darkgreen", size=1.5, alpha=0.6) +
        { if (has_georef) {
              geom_point(data=georef_df_plot, aes(x=longitude, y=latitude),
                         color="darkorange", shape=17, size=2, alpha=0.8)
          } else { NULL }
        } +
        coord_sf(xlim=c(clean_lon_min, clean_lon_max),
                 ylim=c(clean_lat_min, clean_lat_max)) +
        labs(title=clean_title, x="Longitude", y="Latitude") +
        annotation_scale(location="bl") +
        annotation_north_arrow(location="tl",
                               height=unit(0.8,"cm"), width=unit(0.8,"cm")) +
        theme_minimal() +
        theme(plot.title=element_text(size=10))

    combined <- grid.arrange(p1, p2, ncol=2, top=taxa)
    print(combined)
"""

# useful functions

"""
    plot_all_occurrence_maps(;
        traits_path    = "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
        raw_dir        = "data/occurrence_data/pt_occs_raw",
        clean_dir      = "data/occurrence_data/pt_occs_clean",
        georef_dir     = "data/occurrence_data/pt_occs_georeferenced",
        shapefile_path = "data/occurrence_data/bot_country_shapefiles/level3.shp",
        pdf_file       = "data/occurrence_data/occurrence_maps_before_after.pdf",
        date_suffix    = "2026_04_23")

Iterate over every taxon in `traits_path`, build a before/after occurrence map
for each, and write all pages to `pdf_file`.  Returns the number of taxa
successfully plotted.
"""
function plot_all_occurrence_maps(;
    traits_path="data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
    raw_dir="data/occurrence_data/pt_occs_raw",
    clean_dir="data/occurrence_data/pt_occs_clean",
    georef_dir="data/occurrence_data/pt_occs_georeferenced",
    shapefile_path="data/occurrence_data/bot_country_shapefiles/level3.shp",
    pdf_file="data/occurrence_data/occurrence_maps_before_after.pdf",
    date_suffix="2026_04_23"
)
    taxa_traits, taxa_to_nativerange_dict = load_taxa_data(traits_path)
    load_bot_regions(shapefile_path)

    @rput pdf_file
    R"pdf(pdf_file, width=14, height=7)"

    plotted_count = 0

    for (idx, taxa) in enumerate(taxa_traits.scientificName)
        println("Processing $idx/$(length(taxa_traits.scientificName)): $taxa")

        filename = replace(taxa, " " => "_")
        raw_file = joinpath(raw_dir, "$(filename)-$(date_suffix).csv")
        clean_file = resolve_clean_file(filename, clean_dir)

        if !isfile(raw_file) || isempty(clean_file)
            println("  Skipping - files not found")
            continue
        end

        raw_df = filter(
            row -> !ismissing(row.latitude) && !ismissing(row.longitude),
            DataFrame(CSV.File(raw_file))
        )
        clean_df = DataFrame(CSV.File(clean_file))

        if nrow(raw_df) == 0 || nrow(clean_df) == 0
            println("  Skipping - no valid coordinates")
            continue
        end

        georef_df = load_georef_df(filename, georef_dir)
        country_codes = get(taxa_to_nativerange_dict, taxa, String[])

        push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)
        reval(RPLOT_BLOCK)
        plotted_count += 1
    end

    R"dev.off()"
    println("\nDone. Plotted $plotted_count taxa → $pdf_file")
    return plotted_count
end


"""
    plot_species(taxa;
        traits_path    = "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
        raw_dir        = "data/occurrence_data/pt_occs_raw",
        clean_dir      = "data/occurrence_data/pt_occs_clean",
        georef_dir     = "data/occurrence_data/pt_occs_georeferenced",
        shapefile_path = "data/occurrence_data/bot_country_shapefiles/level3.shp",
        date_suffix    = "2026_04_23")

Display a before/after occurrence map for a single `taxa` string in the R
graphics viewer (no PDF output).

# Example
```julia
plot_species("Calceolaria alba")
```
"""
function plot_species(taxa::String;
    traits_path="data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
    raw_dir="data/occurrence_data/pt_occs_raw",
    clean_dir="data/occurrence_data/pt_occs_clean",
    georef_dir="data/occurrence_data/pt_occs_georeferenced",
    shapefile_path="data/occurrence_data/bot_country_shapefiles/level3.shp",
    date_suffix="2026_04_23"
)
    _, taxa_to_nativerange_dict = load_taxa_data(traits_path)
    load_bot_regions(shapefile_path)

    filename = replace(taxa, " " => "_")
    raw_file = joinpath(raw_dir, "$(filename)-$(date_suffix).csv")
    clean_file = resolve_clean_file(filename, clean_dir)

    isfile(raw_file) || error("Raw file not found: $raw_file")
    isempty(clean_file) && error("No cleaned or georef_merged file found for: $taxa")

    raw_df = filter(
        row -> !ismissing(row.latitude) && !ismissing(row.longitude),
        DataFrame(CSV.File(raw_file))
    )
    clean_df = DataFrame(CSV.File(clean_file))

    nrow(raw_df) > 0 || error("No valid coordinates in raw file for $taxa")
    nrow(clean_df) > 0 || error("No rows in cleaned file for $taxa")

    georef_df = load_georef_df(filename, georef_dir)
    country_codes = get(taxa_to_nativerange_dict, taxa, String[])

    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)
    reval(RPLOT_BLOCK)

    println("Plot displayed for: $taxa")
    return nothing
end

"""
    prepare_geolocate_files(;
        input_dir  = "data/occurrence_data/pt_occs_to_georeference",
        output_dir = "data/occurrence_data/pt_occs_to_georeference_geolocate_format")
 
Read every CSV in `input_dir`, reformat columns to GEOLocate batch input
format, and write the result to `output_dir` (same filename).
 
Expected input columns: locality, country, stateProvince, county, latitude,
longitude, ID, scientificName, basisOfRecord.
 
Output columns (GEOLocate format): "locality string", country, state, county,
latitude, longitude, "correction status", precision, "error polygon",
"multiple results", ID, name, basis.
"""
function prepare_geolocate_files(;
    input_dir="data/occurrence_data/pt_occs_to_georeference",
    output_dir="data/occurrence_data/pt_occs_to_georeference_geolocate_format"
)
    isdir(output_dir) || mkpath(output_dir)

    files = DataFrames.filter(f -> endswith(f, ".csv"), readdir(input_dir))

    if isempty(files)
        println("No CSV files found in $input_dir")
        return 0
    end

    converted_count = 0

    for (idx, file) in enumerate(files)
        println("Processing $idx/$(length(files)): $file")

        input_file = joinpath(input_dir, file)
        output_file = joinpath(output_dir, file)

        @rput input_file output_file

        R"""
        rawdf_GeoRef <- read.csv(input_file)

        if (nrow(rawdf_GeoRef) == 0) {
            cat("  Skipping - empty file\n")
            next
        }

        rawdf_GeoRef <- rawdf_GeoRef %>%
            dplyr::select("locality string" = locality,
                          country,
                          state = stateProvince,
                          county,
                          latitude,
                          longitude,
                          ID,
                          name = scientificName,
                          basis = basisOfRecord)

        rawdf_GeoRef$'correction status' <- ""
        rawdf_GeoRef$precision           <- ""
        rawdf_GeoRef$'error polygon'     <- ""
        rawdf_GeoRef$'multiple results'  <- ""

        rawdf_GeoRef2 <- rawdf_GeoRef[, c("locality string", "country",
                                           "state", "county", "latitude",
                                           "longitude", "correction status",
                                           "precision", "error polygon",
                                           "multiple results", "ID",
                                           "name", "basis")]

        write.csv(rawdf_GeoRef2, output_file, row.names = FALSE)
        """

        converted_count += 1
    end

    println("\nDone. Converted $converted_count files → $output_dir")
    return converted_count
end