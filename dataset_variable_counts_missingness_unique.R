# Install these once if you haven't:
# install.packages(c("readxl", "dplyr", "tidyr", "purrr"))

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

#load the data--
# Path to your Excel file
file_path <- "your_file.xlsx"   # <-- change this

# Read the sheet (change sheet = if needed)
raw_df <- read_excel(
  path  = file_path,
  range = "A:DB"   # only columns A through DB
)

#treat blank strings as NA--
df <- raw_df %>%
  mutate(across(
    .cols = everything(),
    .fns  = ~ {
      if (is.character(.x)) {
        x <- trimws(.x)
        x[x == ""] <- NA
        x
      } else {
        .x
      }
    }
  ))
#per variable summary
summary_by_var <- df %>%
  summarise(
    across(
      .cols = everything(),
      .fns = list(
        non_missing = ~ sum(!is.na(.)),
        missing     = ~ sum(is.na(.)),
        unique_vals = ~ n_distinct(., na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to  = c("variable", ".value"),
    names_sep = "__"
  )

#write csv of summary
write.csv(summary_by_var, "variable_summary.csv", row.names = FALSE)

#builds a list of tibbles, one per variable, with value → count.
freq_list <- imap(
  df,
  ~ tibble(
    variable = .y,
    value    = as.character(.x)
  ) %>%
    mutate(is_na = is.na(value)) %>%
    group_by(variable, value, is_na) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
)

#makes frequency list one giant long table of all frequencies
freq_all <- bind_rows(freq_list)

# Example: look at frequencies for a specific variable
freq_list[["BF"]]      # if you have a column literally named "BF"
# or freq_list[[1]]    # first column


