install.packages("Boruta")      # Install Boruta for feature selection
install.packages("data.table")  # Install data.table for efficient data loading
library(Boruta)                 # Load Boruta
library(data.table)             # Load data.table

# Load the training data
trainSet <- as.data.frame(fread("/Users/swetarai/GSE_merging_data/train.csv"))
View(trainSet)
# View the structure of the data
str(trainSet)

# Separate predictors (x) and response (y)
# Assuming 'group' is the response variable and the rest are predictors
x <- trainSet[, -which(names(trainSet) %in% c("Disease_state", "Gender","batch","V1","tissue","platform","sample"))]  # Exclude 'group' and 'sample' columns
y <- trainSet$Disease_state  # Response variable

# Convert response to factor (if it's not already)
y <- as.factor(y)

# Check the levels of the response variable
levels(y)

set.seed(123)  # Set seed for reproducibility

# Run Boruta
boruta_output <- Boruta(x, y, 
                        doTrace = 2)  # doTrace = 2 prints detailed output

# Print the Boruta output
print(boruta_output)

# Get the list of confirmed features
confirmed_features <- getSelectedAttributes(boruta_output, withTentative = FALSE)
print(confirmed_features)

# Get the list of tentative features
tentative_features <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(tentative_features)

# Get a summary of the Boruta results
boruta_summary <- attStats(boruta_output)
print(boruta_summary)

# Plot feature importance
plot(boruta_output, 
     cex.axis = 0.7,  # Adjust axis label size
     las = 2,         # Rotate axis labels
     xlab = "",       # Remove x-axis label
     main = "Feature Importance")  # Add a title


# Get the list of confirmed features
confirmed_features <- getSelectedAttributes(boruta_output, withTentative = FALSE)

# Print the confirmed features
print(confirmed_features)


# Save the confirmed features to a .txt file
write.txt(confirmed_features, 
            file = "/Users/swetarai/GSE_merging_data/boruta_selected_genes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)




