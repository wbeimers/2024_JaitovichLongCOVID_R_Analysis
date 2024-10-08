---------------------------------------------------------------
  # .......................................................  ######
#"Simple Linear Model"####
# .......................................................  ######
-------------------------------------------------------------------------
  
  transposed_norm_data <- t(expression_log2_normalized_runorder)
  
  transposed_norm_data1 <- as.data.frame(transposed_norm_data)
  colnames(transposed_norm_data1) <- Features$feature
  transposed_norm_data1$Filename <- rownames(transposed_norm_data1)
  rownames(transposed_norm_data1) <- NULL
  
  
  
  # Merge the transposed data frame with Final_metadata based on the common sample identifier (Filename)
  Linear_model_data <- merge(Final_metadata, transposed_norm_data1, by = "Filename")
  Linear_model_data <- read.csv("Linear_model_data_noQC.csv", header = TRUE)
  
  #Removing any weird symbols in my column names, especially the name of the features in order to make the for loop work
  # Set the locale to "C" for character handling
  Sys.setlocale(locale = "C")
  
  # Assuming Linear_model_data is your data frame
  
  # Remove the question mark symbol from all column names
  colnames(Linear_model_data) <- gsub("[^A-Za-z0-9_]", "", colnames(Linear_model_data))
  
  # Create an empty list to store the linear models
  linear_models <- list()
  
  # Loop through each intensity column and fit a linear model
  for (i in 28:ncol(Linear_model_data)) {  # Assuming that  intensity columns start from column 28
    intensity_col <- colnames(Linear_model_data)[i]  # Get the name of the current intensity column
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex ")  # Create the formula (first variable after ~ matters most)
    linear_models[[intensity_col]] <- lm(formula, data = Linear_model_data)  # Fit the linear model
  }
  
  
  # Summarize the model
  # Create an empty data frame to store summary information
  summary_df <- data.frame(Intensity_Column = character(),
                           Coefficient = numeric(),
                           Std.Error = numeric(),
                           t.Value = numeric(),
                           P.Value = numeric(),
                           stringsAsFactors = FALSE)
  
  # Loop through each intensity column and store summary information
  for (i in 1:7812) {
    intensity_col <- paste("V", i, sep = "")
    
    # Check if the linear model exists in the list
    if (intensity_col %in% names(linear_models)) {
      summary_data <- summary(linear_models[[intensity_col]])
      
      # Extract relevant information
      coef_data <- coef(summary_data)
      coef_row <- data.frame(Intensity_Column = intensity_col, coef_data)
      
      # Append to the summary data frame
      summary_df <- rbind(summary_df, coef_row)
    }
  }
  
  
  
  #Edit colnames of summary_df
  colnames(summary_df) <-  c("Feature_Intensity",	"Coefficient",	"Std_Error",	"t_value",	"p_value")
  
  #adjusted p-value, Perform FDR correction
  summary_df$p_adjusted  <- adjusted_p_values <- p.adjust(summary_df$p_value, method = "BH")
  
  #histogram
  hist(summary_df$p_adjusted, main= "Histogram of adjusted p-value linear model")
  
  
  # Save the summary data frame to a CSV file
  write.csv(summary_df, file = "linear_model_summary.csv", row.names = TRUE)
  
  
  
  #Anova test summary
  
  
  
  ---------------------------------------------------------------
    # .......................................................  ######
  ## Mixed Model (1|Monkey) ####
  # .......................................................  ######
  -------------------------------------------------------------------------
    
    
    #### Model1: Temporal Changes within the Same Monkey: Temporal Changes within the Same Monkey: To investigate how metabolites change over time within the same monkey, you should focus on modeling the temporal patterns. You can use a linear mixed-effects model with Monkey_ID as a random effect to account for the correlation between repeated measurements within the same monkey.
    
    
    Linear_model_data_noQC <- Linear_model_data[Linear_model_data$Not_a_sample == FALSE,]
  Linear_model_data_noQC <- data.frame(Linear_model_data_noQC)
  
  
  
  #1) Complex Model: 3 interaction terms####
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ PullDate + Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_Diet_complex.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + Sex  + Diet:Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_Age_complex.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~Diet + PullDate+ Diet:PullDate + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_Sex_complex.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_AgeAndSex_complex.csv", row.names = TRUE)
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_DietAndSex_complex.csv", row.names = TRUE)
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:Sex + PullDate:Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_ DietAndAge_complex.csv", row.names = TRUE)
  
  
  
  
  
  
  #2) Simple Model: no interaction terms####
  
  # a) Sex #####
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + PullDate + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_simple_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_simple_model, file = "lrt_results_Sex_simple.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  
  
  #b) Age####
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_simple_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_simple_model, file = "lrt_results_Age_simple.csv", row.names = TRUE)
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + (1 | Monkey_ID)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ PullDate + Sex + (1 | Monkey_ID)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_simple_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_simple_model, file = "lrt_results_Diet_simple.csv", row.names = TRUE)
  
  
  
  
  #Reformat Excel Results####
  
  Diet_Simple = read.csv("lrt_results_diet_simple.csv", header=TRUE)
  Sex_Simple=read.csv("lrt_results_Sex_simple.csv", header=TRUE)
  Age_Simple=read.csv("lrt_results_Age_simple.csv", header=TRUE)
  
  
  Diet_Complex=read.csv("lrt_results_Diet_complex.csv", header=TRUE)
  Sex_Complex=read.csv("lrt_results_Sex_complex.csv", header=TRUE)
  Age_Complex=read.csv("lrt_results_Age_complex.csv", header=TRUE)
  
  Diet_Age_Complex=read.csv("lrt_results_ DietAndAge_complex.csv", header=TRUE)
  Diet_Sex_Complex=read.csv("lrt_results_DietAndSex_complex.csv", header=TRUE)
  Age_Sex_Complex=read.csv("lrt_results_AgeAndSex_complex.csv", header=TRUE)
  
  
  
  Diet_Simple <- Diet_Simple[, -c(1,2,5,8)]
  Sex_Simple <- Sex_Simple[, -c(1,2,5,8)]
  Age_Simple <- Age_Simple[, -c(1,2,5,8)]
  
  Diet_Complex <- Diet_Complex[, -c(1,2,5,8)]
  Sex_Complex <- Sex_Complex[, -c(1,2,5,8)]
  Age_Complex <- Age_Complex[, -c(1,2,5,8)]
  
  Diet_Age_Complex <- Diet_Age_Complex[, -c(1,2,5,8)]
  Diet_Sex_Complex <- Diet_Sex_Complex[, -c(1,5,8)]
  Age_Sex_Complex <- Age_Sex_Complex[, -c(1,2,5,8)]
  
  
  # Reformat
  Diet_Simple <- Diet_Simple %>%
    filter(Model == "DietR")
  
  Sex_Simple <- Sex_Simple %>%
    filter(Model == "SexM")
  
  
  Age_Simple <- Age_Simple %>%
    filter(Model == "PullDate")
  
  
  
  Diet_Complex <- Diet_Complex %>%
    filter(Model == "DietR")
  
  Sex_Complex <- Sex_Complex %>%
    filter(Model == "SexM")
  
  
  Age_Complex <- Age_Complex %>%
    filter(Model == "PullDate")
  
  
  
  
  Diet_Age_Complex <- Diet_Age_Complex %>%
    filter(Model == "DietR:PullDate")
  
  Diet_Sex_Complex <- Diet_Sex_Complex %>%
    filter(Model == "DietR:SexM")
  
  
  Age_Sex_Complex <- Age_Sex_Complex %>%
    filter(Model == "PullDate:SexM")
  
  
  Diet_Simple$adjusted_p_value <- p.adjust(Diet_Simple$p_value, method = "BH")
  Sex_Simple$adjusted_p_value <- p.adjust(Sex_Simple$p_value, method = "BH")
  Age_Simple$adjusted_p_value <- p.adjust(Age_Simple$p_value, method = "BH")
  Diet_Complex$adjusted_p_value <- p.adjust(Diet_Complex$p_value, method = "BH")
  Sex_Complex$adjusted_p_value <- p.adjust(Sex_Complex$p_value, method = "BH")
  Age_Complex$adjusted_p_value <- p.adjust(Age_Complex$p_value, method = "BH")
  Diet_Age_Complex$adjusted_p_value <- p.adjust(Diet_Age_Complex$p_value, method = "BH")
  Diet_Sex_Complex$adjusted_p_value <- p.adjust(Diet_Sex_Complex$p_value, method = "BH")
  Age_Sex_Complex$adjusted_p_value1 <- p.adjust(Age_Sex_Complex$p_value, method = "BH")
  
  
  hist(Diet_Simple$p_value)
  hist(Diet_Complex$p_value)
  
  hist(Sex_Simple$p_value)
  hist(Sex_Complex$p_value)
  
  
  hist(Age_Simple$p_value)
  hist(Age_Complex$p_value)
  
  
  
  write.csv(Diet_Simple, file = "lrt_results_diet_simple_reformatted.csv", row.names = TRUE)
  write.csv(Sex_Simple, file = "lrt_results_Sex_simple_reformatted.csv", row.names = TRUE)
  write.csv(Age_Simple, file = "lrt_results_Age_simple_reformatteed.csv", row.names = TRUE)
  write.csv(Diet_Complex, file = "lrt_results_Diet_complex_reformatted.csv", row.names = TRUE)
  write.csv(Sex_Complex, file = "lrt_results_Sex_complex_reformatted.csv", row.names = TRUE)
  write.csv(Age_Complex, file = "lrt_results_Age_complex_reformatted.csv", row.names = TRUE)
  write.csv(Diet_Age_Complex, file = "lrt_results_ DietAndAge_complex_reformatted.csv", row.names = TRUE)
  write.csv(Diet_Sex_Complex, file = "lrt_results_DietAndSex_complex_reformatted.csv", row.names = TRUE)
  write.csv(Age_Sex_Complex, file = "lrt_results_AgeAndSex_complex_reformatted.csv", row.names = TRUE)
  
  
  
  
  
  ---------------------------------------------------------------
    # .......................................................  ######
  ## Mixed Model (1|Cohort) ####
  # .......................................................  ######
  -------------------------------------------------------------------------
    
    
    #1) Complex Model####
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Cohort)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + PullDate + Diet:PullDate + (1 | Cohort)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  lrt_results_df$adjusted_p_value <- p.adjust(lrt_results_df$p_value, method = "BH") #adjusted p_value
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_Sex_complex_RandomCohort.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Cohort)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~  PullDate + Sex  + PullDate:Sex + (1 | Cohort)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_Diet_complex_RandomCohort.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + Diet:PullDate + Diet:Sex + PullDate:Sex + (1 | Cohort)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~  Diet +  Sex + Diet:Sex + (1 | Cohort)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_complex_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_complex_model, file = "lrt_results_Age_complex_RandomCohort.csv", row.names = TRUE)
  
  
  
  
  
  
  
  
  
  #2) Simple Model####
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + (1 | Cohort)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + PullDate + (1 | Cohort)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_simple_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_simple_model, file = "lrt_results_Sex_simple_RandomCohort.csv", row.names = TRUE)
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + (1 | Cohort)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ Diet + Sex + (1 | Cohort)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_simple_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_simple_model, file = "lrt_results_Age_simple_RandomCohort.csv", row.names = TRUE)
  
  
  
  
  lrt_results <- lapply(28:ncol(Linear_model_data_noQC), function(i) {
    intensity_col <- colnames(Linear_model_data_noQC)[i]
    
    # Fit a model with random effect (1 | Monkey_ID)
    formula <- paste(intensity_col, "~ Diet + PullDate + Sex  + (1 | Cohort)")
    set.seed(107)
    temporal_model <- lmer(formula, data = Linear_model_data_noQC)
    
    # Fit a null model 
    null_formula <- paste(intensity_col, "~ PullDate + Sex + (1 | Cohort)")
    set.seed(107)
    null_model <- lmer(null_formula, data = Linear_model_data_noQC)
    
    # Perform the likelihood ratio test
    
    lrt_result <- anova(null_model, temporal_model)
    
    str(lrt_result)
    
    # Extract the p-value
    p_value <- lrt_result$"Pr(>Chisq)"[2]
    
    # Extract effect sizes (coefficient estimates) and their standard errors
    coefficients <- summary(temporal_model)$coefficients
    effect_sizes <- coefficients[, "Estimate"]
    std_errors <- coefficients[, "Std. Error"]
    
    # Create a data frame for LRT results, including effect sizes and standard errors
    lrt_df <- data.frame(
      Intensity_Column = intensity_col,
      p_value = p_value,
      Effect_Size = effect_sizes,
      Std_Error = std_errors
    )
    
    return(lrt_df)
  })
  
  
  
  # Combine the LRT results into a single data frame
  lrt_results_df <- do.call(rbind, lrt_results)
  
  
  
  #Add Model information in a column
  lrt_results_df$Model <- row.names(lrt_results_df)
  
  # Remove unique numbers from Model name
  lrt_results_df$Model <- gsub("[0-9]", "", lrt_results_df$Model)
  
  
  
  # Map Feature names
  result_simple_model <- lrt_results_df %>%
    left_join(Features, by = c("Intensity_Column" = "Feature_vector"))
  
  # The "feature_name" column from Features is now added to lrs_results_df
  
  
  # Save the ANOVA results to a CSV file
  write.csv(result_simple_model, file = "lrt_results_Diet_simple_RandomCohort.csv", row.names = TRUE)
  
  
  
  
  #Reformat Excel Results####
  
  Diet_complex_cohort = read.csv("lrt_results_Diet_complex_RandomCohort.csv", header=TRUE)
  Sex_complex_Cohort =read.csv("lrt_results_Sex_complex_RandomCohort.csv", header=TRUE)
  Age_complex_cohort=read.csv("lrt_results_Age_complex_RandomCohort.csv", header=TRUE)
  
  Diet_simple_cohort =read.csv("lrt_results_Diet_simple_RandomCohort.csv", header=TRUE)
  Sex_simple_cohort=read.csv("lrt_results_Sex_simple_RandomCohort.csv", header=TRUE)
  Age_simple_cohort=read.csv("lrt_results_Age_simple_RandomCohort.csv", header=TRUE)
  
  
  Diet_simple_cohort <- Diet_simple_cohort[, -c(1,2,5,8)]
  Sex_simple_cohort <- Sex_simple_cohort[, -c(1,2,5,8)]
  Age_simple_cohort <- Age_simple_cohort[, -c(1,2,5,8)]
  
  Diet_complex_cohort <- Diet_complex_cohort[, -c(1,2,5,8)]
  Sex_complex_Cohort <- Sex_complex_Cohort[, -c(1,2,5,8)]
  Age_complex_cohort <- Age_complex_cohort[, -c(1,2,5,8)]
  
  
  # Reformat
  Diet_simple_cohort <- Diet_simple_cohort %>%
    filter(Model == "DietR")
  
  Sex_simple_cohort <- Sex_simple_cohort %>%
    filter(Model == "SexM")
  
  
  Age_simple_cohort <- Age_simple_cohort %>%
    filter(Model == "PullDate")
  
  
  
  Diet_complex_cohort <- Diet_complex_cohort %>%
    filter(Model == "DietR")
  
  Sex_complex_Cohort <- Sex_complex_Cohort %>%
    filter(Model == "SexM")
  
  
  Age_complex_cohort <- Age_complex_cohort %>%
    filter(Model == "PullDate")
  
  
  Diet_simple_cohort$adjusted_p_value <- p.adjust(Diet_simple_cohort$p_value, method = "BH")
  Sex_simple_cohort$adjusted_p_value <- p.adjust(Sex_simple_cohort$p_value, method = "BH")
  Age_simple_cohort$adjusted_p_value <- p.adjust(Age_simple_cohort$p_value, method = "BH")
  Diet_complex_cohort$adjusted_p_value <- p.adjust(Diet_complex_cohort$p_value, method = "BH")
  Sex_complex_Cohort$adjusted_p_value <- p.adjust(Sex_complex_Cohort$p_value, method = "BH")
  Age_complex_cohort$adjusted_p_value <- p.adjust(Age_complex_cohort$p_value, method = "BH")
  
  
  hist(Diet_simple_cohort$p_value)
  hist(Diet_complex_cohort$p_value)
  
  hist(Sex_simple_cohort$p_value)
  hist(Sex_complex_Cohort$p_value)
  
  
  hist(Age_simple_cohort$p_value)
  hist(Age_complex_cohort$p_value)
  
  
  
  
  write.csv(Diet_simple_cohort, file = "lrt_results_diet_simple_COHORT_reformatted.csv", row.names = TRUE)
  write.csv(Sex_simple_cohort, file = "lrt_results_Sex_simple_COHORT_reformatted.csv", row.names = TRUE)
  write.csv(Age_simple_cohort, file = "lrt_results_Age_simple_COHORT_reformatteed.csv", row.names = TRUE)
  write.csv(Diet_complex_cohort, file = "lrt_results_Diet_complex_COHORT_reformatted.csv", row.names = TRUE)
  write.csv(Sex_complex_Cohort, file = "lrt_results_Sex_complex_COHORT_reformatted.csv", row.names = TRUE)
  write.csv(Age_complex_cohort, file = "lrt_results_Age_complex_COHORT_reformatted.csv", row.names = TRUE)
  
  
  
  
  
  
  read.csv("lrt_results_diet_simple_COHORT_reformatted.csv", header = TRUE)
  read.csv("lrt_results_Sex_simple_COHORT_reformatted.csv", header = TRUE)
  read.csv("lrt_results_Age_simple_COHORT_reformatteed.csv", header = TRUE)
  read.csv("lrt_results_Diet_complex_COHORT_reformatted.csv", header = TRUE)
  read.csv("lrt_results_Sex_complex_COHORT_reformatted.csv", header = TRUE)
  read.csv("lrt_results_Age_complex_COHORT_reformatted.csv", header = TRUE)
  
  
  
  ---------------------------------------------------------------------------------
    # .......................................................  ######
  #"Volcano plots after model"####
  # .......................................................  ######
  ------------------------------------------------------------------
    
    
    
    Volcano_data <- Diet_complex
  
  # Set1 significance threshold (e.g., q < 0.05) and create a significant column
  Volcano_data$significant <- ifelse(Volcano_data$adjusted_p_value < 0.05, "yes", "no")
  
  # Create new categorical column for Regulation (adjusted)
  Volcano_data <- Volcano_data %>%
    mutate(Regulation = case_when(Effect_Size >= 1 & adjusted_p_value <= 0.05 ~ "Up-regulated",
                                  Effect_Size <= -1 & adjusted_p_value <= 0.05 ~ "Down-regulated",
                                  Effect_Size > -1 & Effect_Size < 1 & adjusted_p_value <= 0.05 ~ "Significant",
                                  TRUE ~ "Not Significant"))
  
  
  
  
  
  
  cols <- c("Up-regulated" = "#B4B7E8", "Down-regulated" = "#67D2CA", "Significant" = "grey", "Not Significant" = "grey") 
  sizes <- c("Up-regulated" = 2, "Down-regulated" =2, "Not Significant" = 1, "Significant"= 1) 
  alphas <- c("Up-regulated" = 2, "Down-regulated" = 2, "Not Significant" = 0.5,  "Significant" = 0.5)
  
  
  library(forcats)
  
  Volcano_data <- Volcano_data %>%
    mutate(Regulation = fct_relevel(Regulation, "Up-regulated", "Down-regulated", "Significant", "Not Significant"))
  
  
  Volcano_data %>%
    ggplot(aes(x = Effect_Size,
               y = -log10(adjusted_p_value),
               fill = Regulation,
               size = Regulation,
               alpha = Regulation)) + 
    geom_point(shape = 21,
               colour = "grey") + 
    geom_vline(xintercept = c(-1, -1),
               linetype = "dashed", colour= "grey") +
    xlab("Log2(Fold Change)") +
    ylab("-log10(adj. p-value)") +
    ggtitle("Volcano Plot (Diet)") +
    scale_fill_manual(values = cols) +
    scale_size_manual(values = sizes) +
    scale_alpha_manual(values = alphas) +
    scale_x_continuous(breaks = c(seq(-3, 3, 1)),
                       limits = c(-3, 3)) + 
    scale_y_continuous(breaks = c(seq(0, 45, 5)),
                       limits = c(0, 40))+
    geom_label_repel(
      data = subset(Volcano_data, Regulation == "Up-regulated" | Regulation == "Down-regulated"),
      aes(label = feature),  
      force = 10,
      nudge_x = 0,  # Adjust the horizontal position of labels
      nudge_y = 5   # Adjust the vertical position of labels
    ) + 
    theme(
      panel.border = element_rect(colour = "Black", fill = NA, size = 0.5), 
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
  
  
  write.csv(Volcano_data, "Volcano_Sex_Complex.csv", row.names = TRUE)
  
  