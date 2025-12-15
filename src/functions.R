library(dplyr)
library(stringr)

filter_matched_comparisons <- function(df, state_negative, state_positive) {
  
  # OLD Regex: ^(.+)UF_vs_\1DF$  (Strict: Must be Same Treatment)
  # NEW Regex: ^(.+)UF_vs_(.+)DF$ (Relaxed: Can be Different Treatments)
  
  # Logic:
  # Group 1 (.+) = Treatment Numerator (e.g., ALOD4)
  # Group 2 (.+) = Treatment Denominator (e.g., OlyA)
  pattern <- paste0("^(.+)", state_negative, "_vs_(.+)", state_positive, "$")
  
  df_filtered <- df %>%
    filter(str_detect(Label, pattern)) %>%
    mutate(
      # Extract Treatment Names
      # str_match returns a matrix: [FullMatch, Group1, Group2]
      matches = str_match(Label, pattern),
      Treatment_Neg = matches[,2], # The Left side (Denominator)
      Treatment_Pos = matches[,3], # The Right side (Numerator)
      
      # Create a clean Condition Label
      # If they are same: "ALOD4"
      # If mixed: "ALOD4 / OlyA"
      Condition = ifelse(Treatment_Neg == Treatment_Pos, 
                         Treatment_Neg, 
                         paste0(Treatment_Pos, " vs ", Treatment_Neg)),
      
      # Explicitly record direction
      Direction = paste0(state_positive, " (Num) / ", state_negative, " (Denom)")
    ) %>%
    # Clean up the temp "matches" column if desired, though dplyr usually drops it
    select(-matches, -Treatment_Neg, -Treatment_Pos)
  
  return(df_filtered)
}