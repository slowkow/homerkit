#' Read known and novel motif outputs from a HOMER output folder.
#'
#' @export
#' @importFrom stringr
#'   str_detect
#'   str_replace
#'   str_split_fixed
#' @importFrom magrittr
#'   %<>%
#' @importFrom readr
#'   read_tsv
#'   cols
#'   col_character
#'   col_character
#'   col_double
read_homer_output <- function(folder = "", max_files = 100) {

  # Put all of the data in a list.
  retval <- list()

  # Genes with novel motifs in their promoters --------------------------------
  novel_motif_files <- Sys.glob(
    paths = c(
      file.path(folder, "homerResults", "motif*.motif.tsv"),
      file.path(folder, "homerResults", "motif*.motif.tsv.gz")
    )
  )

  if (length(novel_motif_files) > 0) {

    # Exclude some of the results that are labeled "RV" or "similar".
    novel_motif_files <- novel_motif_files[
      !str_detect(novel_motif_files, "RV\\.") &
      !str_detect(novel_motif_files, "\\.similar")
    ]

    message(
      "Found ",
      length(novel_motif_files),
      " motif.tsv[.gz] files in ",
      file.path(folder, "homerResults")
    )
    n_files <- min(length(novel_motif_files), max_files)
    message(
      "Reading ", n_files, " of them..."
    )
    novel_motif_files <- head(novel_motif_files, n_files)

    dat_novel <- do.call(
      rbind.data.frame,
      lapply(novel_motif_files, read_annotatePeaks_tsv)
    )
    # dat_novel <- dat_novel[!is.na(dat_novel$distance_to_peak),]
    # Order by motif name like "motif1" and then by start position.
    order_novel <- order(
      as.numeric(substr(dat_novel$motif, 6, 10)),
      dat_novel$start
    )
    dat_novel <- dat_novel[order_novel,]

    retval[["novel_motif_peaks"]] <- dat_novel
  }

  # Novel motif transcription factors -----------------------------------------
  html_files <- Sys.glob(
    paths = file.path(folder, "homerResults", "motif*.info.html")
  )

  if (length(html_files) > 0) {
    message(
      "Found ",
      length(html_files),
      " motif*.info.html files in ",
      file.path(folder, "homerResults")
    )
    n_files <- min(length(html_files), max_files)
    message(
      "Reading ", n_files, " of them..."
    )
    html_files <- head(html_files, n_files)

    html_novel <- do.call(
      rbind.data.frame,
      lapply(html_files, read_findMotifs_html)
    )
    # Order by motif name like "motif1" and then by match rank.
    order_html <- order(
      as.numeric(substr(html_novel$motif, 6, 10)),
      as.numeric(html_novel$match_rank)
    )
    html_novel <- html_novel[order_html,]

    retval[["novel_motif_tfs"]] <- html_novel
  }

  # Genes with known motifs in their promoters ----------------------------------
  known_motif_files <- Sys.glob(
    paths = c(
      file.path(folder, "knownResults", "*.motif.tsv"),
      file.path(folder, "knownResults", "*.motif.tsv.gz")
    )
  )

  if (length(known_motif_files) > 0) {
    # Exclude some of the results that are labeled "RV" or "similar".
    known_motif_files <- known_motif_files[
      !str_detect(known_motif_files, "RV\\.") &
      !str_detect(known_motif_files, "\\.similar")
    ]

    message(
      "Found ",
      length(known_motif_files),
      " motif.tsv[.gz] files in ",
      file.path(folder, "knownResults")
    )
    n_files <- min(length(known_motif_files), max_files)
    message(
      "Reading ", n_files, " of them..."
    )
    known_motif_files <- head(known_motif_files, n_files)

    dat_known <- do.call(
      rbind.data.frame,
      lapply(known_motif_files, read_annotatePeaks_tsv)
    )
    # dat_known <- dat_known[!is.na(dat_known$distance_to_peak),]
    # Order by motif name like "motif1" and then by start position.
    order_known <- order(
      as.numeric(substr(dat_known$motif, 6, 10)),
      dat_known$start
    )
    dat_known <- dat_known[order_known,]

    retval[["known_motif_peaks"]] <- dat_known
  }

  retval[["known_motif_table"]] <- janitor::clean_names(
    read_tsv(
      file      = file.path(folder, "knownResults.txt"),
      col_types = cols(
        col_character(),
        col_character(),
        col_double(),
        col_double(),
        col_double(),
        col_double(),
        col_character(),
        col_double(),
        col_character()
      )
        #`Motif Name`                                   = col_character(),
        #`Consensus`                                    = col_character(),
        #`P-value`                                      = col_double(),
        #`Log P-value`                                  = col_double(),
        #`q-value (Benjamini)`                          = col_double(),
        #`# of Target Sequences with Motif(of 60)`      = col_double(),
        #`% of Target Sequences with Motif`             = col_character(),
        #`# of Background Sequences with Motif(of 490)` = col_double(),
        #`% of Background Sequences with Motif`         = col_character()
    )
  )

  retval[["novel_motif_table"]] <- read_homerResults_html(
    file.path(folder, "homerResults.html")
  )

  return(retval)
}

#' Read homerResults.html file produced by HOMER.
#'
#' @export
#' @importFrom stringr
#'   str_replace
#'   str_replace_all
#'   str_split_fixed
#'   str_match
#' @importFrom magrittr
#'   %<>%
#'   %>%
read_homerResults_html <- function(file) {
  dat <- janitor::clean_names(
    htmltab::htmltab(
      doc            = file,
      which          = 1,
      header         = 1,
      rm_nodata_cols = FALSE,
      rm_whitespace  = TRUE
    )
  )
  # Rename this column.
  colnames(dat)[colnames(dat) == "log_p_pvalue"] <- "log_pvalue"
  # Drop unused columns.
  dat$rank       <- NULL
  dat$motif_file <- NULL
  dat$motif      <- NULL
  # Drop the link text from this column.
  dat$best_match_details %<>% str_replace(
    ., "More Information \\| Similar Motifs Found", ""
  )
  # Extract the score for the best match.
  dat$best_match_score <- as.numeric(str_match(dat$best_match_details, "\\(([\\.\\d]+)\\)")[,2])
  dat$best_match_details <- str_split_fixed(dat$best_match_details, "\\(", 2)[,1]
  # Rename this column.
  colnames(dat)[colnames(dat) == "best_match_details"] <- "best_match"
  # Convert to numeric.
  dat$p_value               %<>% as.numeric
  dat$log_pvalue            %<>% as.numeric
  dat$percent_of_targets    %<>% str_replace(., "%", "") %>% as.numeric
  dat$percent_of_background %<>% str_replace(., "%", "") %>% as.numeric
  # Split the std bp column into two columns.
  std_bp         <- str_split_fixed(dat$std_bg_std, " ", 2)
  dat$std_bp     <- as.numeric(str_replace(std_bp[,1], "bp", ""))
  dat$bg_std_bp  <- as.numeric(str_replace_all(std_bp[,2], "[\\(\\)bp]", ""))
  dat$std_bg_std <- NULL
  # Set numeric rownames.
  rownames(dat) <- seq(nrow(dat))
  return(dat)
}

#' Read a HOMER motif file.
#'
#' @export
#' @importFrom stringr
#'   str_split_fixed
#'   str_match
#'   str_replace
#'   str_replace_all
#' @importFrom readr
#'   read_tsv
#'   cols
#'   col_character
#'   col_integer
#'   col_double
read_annotatePeaks_tsv <- function(motif_file) {
  dat <- read_tsv(
    file      = motif_file,
    col_types = cols(
      .default          = col_character(),
      `Start`           = col_integer(),
      `End`             = col_integer(),
      `Distance to TSS` = col_integer(),
      `Entrez ID`       = col_integer(),
      `CpG%`            = col_double(),
      `GC%`             = col_double()
    )
  )

  # Fix weird columns.
  command          <- str_split_fixed(colnames(dat)[1], " ", 2)[,2]
  colnames(dat)[1] <- str_split_fixed(colnames(dat)[1], " ", 2)[,1]

  # Discard unannotated rows.
  dat <- dat[!is.na(dat[[22]]),]

  # Column 22 is "<BESTGUESS> Distance From Peak(sequence,strand,conservation)"
  best_guess        <- str_split_fixed(colnames(dat)[22], " ", 2)[,1]
  colnames(dat)[22] <- str_split_fixed(colnames(dat)[22], " ", 2)[,2]

  x                     <- str_match(dat[[22]], "^(-?[0-9]+)\\(([A-Z]+),([+-]),([0-9.]+)\\)")
  dat$distance_to_peak  <- as.numeric(x[,2])
  dat$peak_sequence     <- x[,3]
  dat$peak_strand       <- x[,4]
  dat$peak_conservation <- as.numeric(x[,5])
  dat[[22]]             <- NULL

  dat$best_guess <- best_guess
  dat$motif      <- str_split_fixed(basename(motif_file), "\\.", 2)[,1]

  # Fix column names.
  colnames(dat) <- tolower(colnames(dat))
  colnames(dat) <- str_replace_all(colnames(dat), " ", "_")
  colnames(dat) <- str_replace(colnames(dat), "%", "_percent")
  colnames(dat) <- str_replace(colnames(dat), "/", "_over_")

  # Data types.
  if ("peak_score" %in% colnames(dat)) {
    dat$peak_score <- as.numeric(dat$peak_score)
  }

  #if (reduce) {
  #  cols <- c(
  #    "motif", "peakid", "distance_to_peak", "peak_sequence", "peak_strand"
  #  )
  #  return(list(
  #    best_guess = dat$best_guess[1],
  #    motif_peaks = dat[,cols]
  #  ))
  #}
  return(dat)
}

#' Read a HOMER HTML file.
#'
#' @export
#' @importFrom stringr
#'   str_replace
#'   str_split_fixed
read_findMotifs_html <- function(html_file) {
  doc <- XML::htmlParse(html_file)

  # Get the text inside h4 elements:
  o           <- XML::getNodeSet(doc, "//h4")
  match_names <- sapply(o, XML::xmlValue)

  o <- XML::getNodeSet(doc, "//table")
  x <- lapply(3:length(o), function(i) {
    tb <- o[[i]]
    tb <- tb[["tr"]][["td"]][["table"]]
    t(sapply(tb["tr"], function(x) sapply(x["td"], XML::xmlValue)))
  })

  x <- data.frame(do.call(rbind.data.frame, x[seq(1, length(x), by = 2)]))
  rownames(x) <- seq_len(nrow(x))
  colnames(x) <- c("variable", "value")

  x$variable   <- str_replace(x$variable, ":", "")
  x$variable   <- tolower(str_replace(x$variable, " ", "_"))
  x$match_name <- rep(match_names, each = 5)

  x       <- as.data.frame(reshape2::dcast(x, match_name ~ variable, value.var = "value"))
  x$motif <- str_split_fixed(basename(html_file), "\\.", 2)[,1]

  x$match_rank <- as.numeric(x$match_rank)
  x$offset     <- as.numeric(x$offset)
  x$score      <- as.numeric(x$score)

  # Split the alignment into two parts
  midpoint     <- nchar(x$alignment) / 2
  x$alignment1 <- substr(x$alignment, 1, midpoint)
  x$alignment2 <- substr(x$alignment, midpoint + 1, nchar(x$alignment))
  x$alignment  <- NULL

  return(tibble::as_data_frame(x))
}

