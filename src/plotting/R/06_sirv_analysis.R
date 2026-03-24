# 06_sirv_analysis.R
# SIRV metrics calculation, time parsing utilities, and result export helpers

real_sirvs = c("SIRV101","SIRV102","SIRV103","SIRV105","SIRV106","SIRV107","SIRV108","SIRV109","SIRV201","SIRV202","SIRV203","SIRV204","SIRV205","SIRV206","SIRV301","SIRV302","SIRV303","SIRV304","SIRV305","SIRV306","SIRV307","SIRV308","SIRV309","SIRV310","SIRV311","SIRV403","SIRV404","SIRV405","SIRV406","SIRV408","SIRV409","SIRV410","SIRV501","SIRV502","SIRV503","SIRV504","SIRV505","SIRV506","SIRV507","SIRV508","SIRV509","SIRV510","SIRV511","SIRV512","SIRV601","SIRV602","SIRV603","SIRV604","SIRV605","SIRV606","SIRV607","SIRV608","SIRV609","SIRV610","SIRV611","SIRV612","SIRV613","SIRV614","SIRV615","SIRV616","SIRV617","SIRV618","SIRV701","SIRV702","SIRV703","SIRV704","SIRV705","SIRV706","SIRV708")

n_sirvs = length(real_sirvs)


calculate_sirv_metrics <- function(df){
  
  sirv_df <- df[grepl("^SIRV", df$chrom), ]
  
  sirv_tp <- sirv_df[(sirv_df$associated_transcript %in% real_sirvs & 
                       (sirv_df$subcategory == "reference_match" | 
                          (sirv_df$subcategory == "mono-exon" & sirv_df$diff_to_TSS < 50 & sirv_df$diff_to_TTS < 50))), ]
  
  sirv_tp <- sirv_tp %>%
    group_by(associated_transcript) %>%
    filter(diff_to_TSS + diff_to_TTS == min(diff_to_TSS + diff_to_TTS))
  
  sirv_ptp <- sirv_df %>%
    filter(sirv_df$structural_category %in% c("FSM", "ISM")) %>%
    dplyr::setdiff(sirv_tp)
  
  sirv_fn <- dplyr::setdiff(real_sirvs, sirv_df$associated_transcript)
  
  sirv_fp <- sirv_df %>%
    filter(sirv_df$structural_category %in% c("NIC", "NNC", "Genic\nGenomic", "Fusion"))
  
  SIRV_fsm_ism <- sirv_df[sirv_df$structural_category %in% c("FSM", "ISM"), ]
  
  non_redundant_precision <- NA
  positive_detection_rate <- NA
  false_discovery_rate <- NA
  false_detection_rate <- NA
  redundancy <- NA
  tp_recall <- NA
  ptp_recall <- NA
  
  non_redundant_precision <- nrow(sirv_tp) / nrow(sirv_df)
  
  redundant_precision <- (nrow(sirv_tp) + nrow(sirv_ptp)) / nrow(sirv_df)

  sensitivity <- nrow(sirv_tp) / length(real_sirvs)
  
  unique_detected_transcripts <- length(unique(c(sirv_tp$associated_transcript, sirv_ptp$associated_transcript)))
  positive_detection_rate <- unique_detected_transcripts / length(real_sirvs)

  false_discovery_rate <- (nrow(sirv_df) - nrow(sirv_tp)) / nrow(sirv_df)

  false_detection_rate <- nrow(sirv_fp) / nrow(sirv_df)

  unique_tp_ptp_transcripts <- length(unique(c(sirv_tp$associated_transcript, sirv_ptp$associated_transcript)))
  redundancy <- nrow(SIRV_fsm_ism) / unique_tp_ptp_transcripts
  
  metrics <- tibble(
    Metric = c(
      "Non-redundant Precision",
      "Precision",
      "Sensitivity",
      "Positive Detection Rate",
      "False Discovery Rate",
      "False Detection Rate",
      "Redundancy"
    ),
    Value = c(
      non_redundant_precision,
      redundant_precision,
      sensitivity,
      positive_detection_rate,
      false_discovery_rate,
      false_detection_rate,
      redundancy
    )
  ) %>%
    mutate(
      Value = round(Value * 100, 2),
      Value = paste0(Value, "%")
    )
  
  
  return(list(
    sirv_df = sirv_df,
    metrics = metrics,
    sirv_fn = sirv_fn
  ))
}


nested_to_dataframe <- function(results) {
  imap_dfr(results, function(platform_data, platform) {
    imap_dfr(platform_data, function(method_data, method) {
      imap_dfr(method_data, function(tissue_data, tissue) {
        imap_dfr(tissue_data, function(df, df_name) {
          metrics_wide <- df$metrics %>%
            pivot_wider(names_from = Metric, values_from = Value)
          
          tibble(
            platform = platform,
            method = method,
            tissue = tissue,
            sample = df_name,
            missing_sirvs = list(df$sirv_fn)
          ) %>%
            bind_cols(metrics_wide)
        })
      })
    })
  })
}


parse_time_to_seconds <- function(time_str) {
  if (grepl("-", time_str)) {
    parts <- unlist(strsplit(time_str, "-"))
    days <- as.numeric(parts[1])
    time_part <- parts[2]
  } else {
    days <- 0
    time_part <- time_str
  }
  
  time_components <- unlist(strsplit(time_part, ":"))
  
  hours <- as.numeric(time_components[1])
  minutes <- as.numeric(time_components[2])
  seconds <- as.numeric(time_components[3])
  
  total_seconds <- days * 86400 + hours * 3600 + minutes * 60 + seconds
  return(total_seconds)
}

seconds_to_hms <- function(total_seconds) {
  hours <- floor(total_seconds / 3600)
  remainder <- total_seconds %% 3600
  minutes <- floor(remainder / 60)
  seconds <- remainder %% 60
  
  hours <- as.integer(hours)
  minutes <- as.integer(minutes)
  seconds <- as.integer(seconds)
  
  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}
