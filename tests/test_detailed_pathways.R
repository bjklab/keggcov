# Test script for get_detailed_pathway_info

library(keggcov)

# Define KOs
# K00844: Hexokinase
# K00164: Pyruvate Dehydrogenase
my_kos <- c("K00844", "K00164")

message("Testing get_detailed_pathway_info with K00844 and K00164...")

# Run the function
result <- get_detailed_pathway_info(my_kos)

# Check results
if (!is.null(result)) {
    message("Success! Result dataframe:")
    print(head(result))

    # Check dimensions
    message(paste("Rows:", nrow(result)))
    message(paste("Cols:", ncol(result)))

    # Check columns
    expected_cols <- c("KO", "Pathway_ID", "Pathway_Description")
    if (all(expected_cols %in% colnames(result))) {
        message("Column names are correct.")
    } else {
        stop("Column names are incorrect!")
    }

    # Check for NA descriptions
    if (any(is.na(result$Pathway_Description))) {
        warning("Some Pathway Descriptions are NA. Partial match failure?")
        print(head(result[is.na(result$Pathway_Description), ]))
    } else {
        message("All Pathway Descriptions are populated.")
    }
} else {
    stop("Result is NULL!")
}

message("Test complete.")
