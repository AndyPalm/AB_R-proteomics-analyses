# src/protein_ancova_analysis.R
library(dplyr)
library(stringr)
library(broom)
library(purrr)
library(tidyr)

# --- 1. Prepare Expression Reference (df2.protein) ---
# Goal: Calculate the Mean Expression Intensity for UF and DF states
# Logic: We filter for GROUP "UF" or "DF" (ignoring "static")
df_expr_ref <- df2.protein %>%
  filter(GROUP %in% c("UF", "DF")) %>%
  # Group by Protein and State to get the mean
  group_by(Protein, State = GROUP) %>%
  summarise(Mean_Expr = mean(LogIntensities, na.rm = TRUE), .groups = "drop")

# --- 2. Prepare Labeling Data (df1.protein) ---
# Goal: Parse "ALOD4DF" into Treatment="ALOD4" and State="DF"
df_label_clean <- df1.protein %>%
  # Filter only for the treatment groups (ignore the controls "UF"/"DF" in this file)
  filter(str_detect(GROUP, "ALOD4|OlyA")) %>%
  mutate(
    # Extract Treatment (Start of string)
    Treatment = str_extract(GROUP, "^(ALOD4|OlyA)"),
    # Extract State (End of string)
    State = str_extract(GROUP, "(UF|DF)$")
  ) %>%
  select(Protein, Treatment, State, Label_Intensity = LogIntensities)

# --- 3. Master Merge ---
# Join every Labeling Replicate with its corresponding Mean Expression
master_ancova_df <- df_label_clean %>%
  inner_join(df_expr_ref, by = c("Protein", "State"))

# --- 4. Define the Regression Function ---
# Model: Labeling ~ State + Mean_Expr
# We want to measure the effect of State (UF->DF) while controlling for Expression
fit_ancova <- function(df) {
  # Force UF as reference level so the coefficient represents shift TO DF
  df$State <- factor(df$State, levels = c("UF", "DF"))
  
  tryCatch({
    # The Model: Is there a Labeling shift driven by State, distinct from Expression?
    lm(Label_Intensity ~ State + Mean_Expr, data = df) %>%
      tidy()
  }, error = function(e) return(NULL))
}

# --- 5. Run Regression per Protein & Treatment ---
# Note: We filter for proteins that have data in BOTH states (UF and DF)
results_nested <- master_ancova_df %>%
  group_by(Protein, Treatment) %>%
  # Ensure we have both states represented before fitting
  filter(n_distinct(State) == 2) %>%
  nest() %>%
  mutate(model_fit = map(data, fit_ancova)) %>%
  unnest(model_fit)

# --- 6. Extract Proximity Scores ---
# We are looking for the term "StateDF" (The shift)
proximity_scores <- results_nested %>%
  filter(term == "StateDF") %>%
  group_by(Treatment) %>%
  mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
  arrange(adj.p.value) %>%
  select(Protein, Treatment, Proximity_Shift = estimate, p.value, adj.p.value)

# Print top hits
print(head(proximity_scores))

# Optional: Save results
# out_path <- here::here("output", "ANCOVA_Proximity_Scores.csv")
# write.csv(proximity_scores, out_path, row.names = FALSE)