plotting_packages <- c(
    "repr",
    "data.table",
    "tidyr",
    "dplyr",
    "ggplot2",
    "ggrepel",
    "ggbreak",
    "ggtext",
    "scales",
    "cowplot",
    "RColorBrewer",
    "latex2exp",
    "glue",
    "stringr"
)

differential_expression_packages <- c(
    "TCGAbiolinks",
    "SummarizedExperiment",
    "DESeq2",
    "fgsea",
    "GSVA",
    "GSEABase",
    "biomaRt"
)

# Identify packages absent from the user's currently accessible R
get_missing_packages <- function(package_names) {
    installed_package_names <- rownames(installed.packages())
    setdiff(package_names, installed_package_names)
}

get_available_packages <- function(repos) {
    tryCatch(
        rownames(available.packages(repos = repos)),
        error = function(e) character(0),
        warning = function(w) character(0)
    )
}

# Install any missing packages into a project-local library.
install_missing_packages <- function(location, package_names) {
    local_library <- normalizePath(paste0(location,"/.Rlibs"),
                                   winslash = "/",
                                   mustWork = FALSE)
    dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
    .libPaths(unique(c(local_library, .libPaths())))

    missing_packages <- get_missing_packages(package_names)

    if (length(missing_packages) == 0) {
        # message("All plotting packages are already available in .libPaths().")
        return(invisible(local_library))
    }

    if (is.null(getOption("repos")) || identical(getOption("repos")[["CRAN"]], "@CRAN@")) {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
    }

    cran_packages <- get_available_packages(getOption("repos"))
    packages_on_cran <- intersect(missing_packages, cran_packages)

    if (length(packages_on_cran) > 0) {
        message("Installing CRAN packages into local library: ", paste(packages_on_cran, collapse = ", "))
        install.packages(packages_on_cran, lib = local_library, dependencies = TRUE)
    }

    remaining_packages <- get_missing_packages(missing_packages)

    if (length(remaining_packages) == 0) {
        return(invisible(local_library))
    }

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", lib = local_library, dependencies = TRUE)
    }

    bioc_repos <- tryCatch(BiocManager::repositories(), error = function(e) NULL)
    bioc_packages <- if (is.null(bioc_repos)) character(0) else get_available_packages(bioc_repos)
    packages_on_bioc <- intersect(remaining_packages, bioc_packages)

    if (length(packages_on_bioc) > 0) {
        message("Installing Bioconductor packages into local library: ", paste(packages_on_bioc, collapse = ", "))
        BiocManager::install(packages_on_bioc, lib = local_library, ask = FALSE, update = FALSE)
    }

    still_missing <- get_missing_packages(missing_packages)
    if (length(still_missing) > 0) {
        warning("These packages were not found in CRAN or Bioconductor: ",
                paste(still_missing, collapse = ", "))
    }

    invisible(local_library)
}
