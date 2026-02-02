#' @title Create a Pathway-Based Covariance Matrix from KEGG orthologs (KOs)
#' @description Create a pathway-based variance-covariance matrix from a set of KEGG orthologs (KOs).
#' @details This function leverages the KEGG metabolic pathway database <https://www.genome.jp/kegg/pathway.html> to estimate a variance-covariance structure for sets of KOs identified by high-throughput sequencing. The output matrix can be interpreted through its (1) Diagonal Elements (Variance), which represent the functional versatility of the KO (i.e., a high variance on the diagonal means the KO is involved in many different pathways (a "hub" gene); a low variance means the KO is specialized to only one or two specific pathways); and through its (2) Off-Diagonal Elements (Covariance), which prepresent pathway partnerships. Positive high covariance indicates that the KOs are "pathway partners." They almost always appear together (e.g., subunits of the same complex or sequential enzymes in a linear pathway). Zeros indicated that the KOs share no biological processes (e.g., a DNA replication gene vs. a glycolysis gene).
#' @param ko_list A character vector of KO identifiers (e.g., c("K00844", "K01803"))
#' @param strict_metabolic Logical. If TRUE, limits analysis to only manually drawn reference metabolic pathways ('map' prefix).
#' @return A list containing the covariance matrix and the binary pathway profile matrix.
#' @importFrom utils read.delim
#' @importFrom stats cov
#' @export
#' @examples
#'
#' # Define KOs
#' # K00844, K00845: Hexokinase/Glucokinase (Glycolysis)
#' # K00164: Pyruvate Dehydrogenase (TCA Cycle)
#' # K00500: Fatty Acid Synthase (Fatty Acid Biosynthesis)
#' my_kos <- c("K00844", "K00845", "K00164", "K00500")
#'
#' # Run the function, without metabolic filtering
#'
#' library(keggcov)
#'
#' results <- get_pathway_covariance(my_kos)
#'
#' str(results)
#' rownames(results$binary_matrix)
#'
#' if (!is.null(results)) {
#'     # Print the Covariance Matrix (Rounded for readability)
#'     print(round(results$cov_matrix, 3))
#'
#'     # Visualization
#'     heatmap(results$cov_matrix,
#'         main = "Pathway Co-Membership Covariance",
#'         symm = TRUE,
#'         col = cm.colors(256)
#'     )
#' }
#'
#' # Run the function, with metabolic filtering
#' results <- get_pathway_covariance(my_kos, strict_metabolic = TRUE)
#'
#' str(results)
#'
#' if (!is.null(results)) {
#'     # Print the Covariance Matrix (Rounded for readability)
#'     print(round(results$cov_matrix, 3))
#'
#'     # Visualization
#'     heatmap(results$cov_matrix,
#'         main = "Metabolic Pathway Co-Membership Covariance",
#'         symm = TRUE,
#'         col = cm.colors(256)
#'     )
#' }
#'
get_pathway_covariance <- function(ko_list, strict_metabolic = FALSE) {
    # 1. Fetch KO-Pathway links directly from KEGG API
    #    The 'link' operation returns a tab-separated text file (TSV)
    url <- "https://rest.kegg.jp/link/pathway/ko"

    message(paste("Fetching data from:", url))

    # specific colClasses speeds up reading and ensures text format
    tryCatch(
        {
            link_data <- read.delim(url,
                header = FALSE, col.names = c("ko", "pathway"),
                stringsAsFactors = FALSE, colClasses = "character"
            )
        },
        error = function(e) {
            stop("Failed to fetch data from KEGG API. Check your internet connection.")
        }
    )

    # 2. Clean the data (Remove prefixes 'ko:' and 'path:')
    #    Using base R 'sub' for dependency-free string manipulation
    link_data$ko <- sub("ko:", "", link_data$ko)
    link_data$pathway <- sub("path:", "", link_data$pathway)

    message("Filtering for selected KOs...")

    # 3. Filter for only the KOs provided by the user
    filtered_df <- link_data[link_data$ko %in% ko_list, ]

    if (nrow(filtered_df) == 0) {
        # It is possible the KOs are valid but have no pathway maps (e.g. uncharacterized)
        warning("None of the provided KOs were found associated with any KEGG pathways.")
        return(NULL)
    }

    # 4. Construct the Binary Matrix
    unique_pathways <- unique(filtered_df$pathway)

    # Initialize matrix: Rows = Pathways, Cols = KOs
    pathway_matrix <- matrix(0,
        nrow = length(unique_pathways),
        ncol = length(ko_list)
    )

    rownames(pathway_matrix) <- unique_pathways
    colnames(pathway_matrix) <- ko_list

    # Fill the matrix
    # We iterate through the filtered associations to set 1s
    for (i in seq_len(nrow(filtered_df))) {
        p_id <- filtered_df$pathway[i]
        k_id <- filtered_df$ko[i]

        # Check if the column exists (it should, based on initialization)
        if (k_id %in% colnames(pathway_matrix)) {
            pathway_matrix[p_id, k_id] <- 1
        }
    }

    message(paste("Constructed matrix with", length(unique_pathways), "pathways."))

    # 5. (Optional) Strict Metabolic Filtering
    if (strict_metabolic) {
        message("Applying strict metabolic filtering...")
        # TODO: Implement metabolic filtering logic
        valid_metabolic_maps <- grepl("^map", rownames(pathway_matrix))
        pathway_matrix <- pathway_matrix[valid_metabolic_maps, ]

        message(paste("Restricted pathway matrix to", nrow(pathway_matrix), "metabolic pathways."))

        if (nrow(pathway_matrix) == 0) {
            warning("None of the provided KOs were found associated with any metabolic KEGG pathways. Try strict_metabolic = FALSE")
            return(NULL)
        }
    }

    # 6. Calculate Covariance
    #    We calculate covariance between Columns (KOs)
    cov_matrix <- cov(pathway_matrix)

    return(list(
        cov_matrix = cov_matrix,
        binary_matrix = pathway_matrix
    ))
}

#' @title Get Detailed Pathway Information for KOs
#' @description Retrieves detailed metabolic pathway information for a list of KOs.
#' @param ko_list A character vector of KO identifiers.
#' @return A data frame with columns: KO, Pathway_ID, Pathway_Description.
#' @export
#' @examples
#'
#' library(keggcov)
#' my_kos <- c("K00844", "K00164")
#' df <- get_detailed_pathway_info(my_kos)
#' head(df)
#'
get_detailed_pathway_info <- function(ko_list) {
    # 1. Fetch KO-Pathway links
    url_links <- "https://rest.kegg.jp/link/pathway/ko"
    message(paste("Fetching links from:", url_links))

    tryCatch(
        {
            links <- read.delim(url_links,
                header = FALSE, col.names = c("ko", "pathway"),
                stringsAsFactors = FALSE, colClasses = "character"
            )
        },
        error = function(e) {
            stop("Failed to fetch links from KEGG API.")
        }
    )

    links$ko <- sub("ko:", "", links$ko)
    links$pathway <- sub("path:", "", links$pathway)
    # Normalize pathway ID to numeric string (remove 'ko' and everything before numeric)
    links$pathway_id <- gsub("[^0-9]", "", links$pathway)

    # Filter for user KOs
    links <- links[links$ko %in% ko_list, ]

    if (nrow(links) == 0) {
        warning("No pathways found for the provided KOs.")
        return(NULL)
    }

    # 2. Fetch Pathway Descriptions
    url_names <- "https://rest.kegg.jp/list/pathway"
    message(paste("Fetching pathway names from:", url_names))

    tryCatch(
        {
            names_df <- read.delim(url_names,
                header = FALSE, col.names = c("pathway", "description"),
                stringsAsFactors = FALSE, colClasses = "character"
            )
        },
        error = function(e) {
            stop("Failed to fetch pathway names from KEGG API.")
        }
    )

    names_df$pathway <- sub("path:", "", names_df$pathway)
    # Normalize pathway ID to numeric string (remove 'map' and everything before numeric)
    names_df$pathway_id <- gsub("[^0-9]", "", names_df$pathway)

    # 3. Merge by normalized ID
    result <- merge(links, names_df, by = "pathway_id", all.x = TRUE)

    # Reorder columns
    result <- result[, c("ko", "pathway.x", "description")]
    colnames(result) <- c("KO", "Pathway_ID", "Pathway_Description")

    return(result)
}
