#' Read known and novel motif outputs from a HOMER output folder.
#'
#' @export
read_homer_output <- function(folder = "") {

  # Put all of the data in a list.
  retval <- list()

  # Genes with novel motifs in their promoters --------------------------------
  novel_motif_files <- Sys.glob(
    paths = file.path(folder, "homerResults", "motif*.motif.tsv")
  )

  if (length(novel_motif_files) > 0) {
    novel_motif_files <- novel_motif_files[
      !stringr::str_detect(novel_motif_files, "RV\\.") &
      !stringr::str_detect(novel_motif_files, "\\.similar")
    ]

    dat_novel <- do.call(
      rbind.data.frame,
      lapply(novel_motif_files, read_annotatePeaks_tsv)
    )
    dat_novel <- dat_novel[!is.na(dat_novel$distance_to_peak),]
    # Order by motif name like "motif1" and then by start position.
    order_novel <- order(
      as.numeric(substr(dat_novel$motif, 6, 10)),
      dat_novel$start
    )
    dat_novel <- dat_novel[order_novel,]

    retval[["novel_motif_targets"]] <- dat_novel
  }

  # Novel motif transcription factors -----------------------------------------
  html_files <- Sys.glob(
    paths = file.path(folder, "homerResults", "motif*.info.html")
  )

  if (length(html_files) > 0) {
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
    paths = file.path(folder, "knownResults", "*.motif.tsv")
  )

  if (length(known_motif_files) > 0) {
    dat_known <- do.call(
      rbind.data.frame,
      lapply(known_motif_files, read_annotatePeaks_tsv)
    )
    dat_known <- dat_known[!is.na(dat_known$distance_to_peak),]
    # Order by motif name like "motif1" and then by start position.
    order_known <- order(
      as.numeric(substr(dat_known$motif, 6, 10)),
      dat_known$start
    )
    dat_known <- dat_known[order_known,]

    retval[["known_motif_targets"]] <- dat_known
  }

  return(retval)
}

#' Read a HOMER motif file.
#'
#' @export
read_annotatePeaks_tsv <- function(motif_file) {
  dat <- readr::read_tsv(
    file = motif_file,
    col_types = readr::cols(
      .default = readr::col_character(),
      Start = readr::col_integer(),
      End = readr::col_integer(),
      `Not Used` = readr::col_integer(),
      `Distance to TSS` = readr::col_integer(),
      `Entrez ID` = readr::col_integer(),
      `CpG%` = readr::col_double(),
      `GC%` = readr::col_double()
    )
  )

  # Fix weird columns.
  command <- stringr::str_split_fixed(colnames(dat)[1], " ", 2)[,2]
  colnames(dat)[1] <- stringr::str_split_fixed(colnames(dat)[1], " ", 2)[,1]

  # Column 22 is "<BESTGUESS> Distance From Peak(sequence,strand,conservation)"
  best_guess <- stringr::str_split_fixed(colnames(dat)[22], " ", 2)[,1]
  colnames(dat)[22] <- stringr::str_split_fixed(colnames(dat)[22], " ", 2)[,2]

  x <- stringr::str_match(dat[[22]], "^(-?[0-9]+)\\(([A-Z]+),([+-]),([0-9.]+)\\)")
  dat$distance_to_peak <- as.numeric(x[,2])
  dat$peak_sequence <- x[,3]
  dat$peak_strand <- x[,4]
  dat$peak_conservation <- as.numeric(x[,5])
  dat[[22]] <- NULL

  dat$best_guess <- best_guess
  dat$motif <- stringr::str_split_fixed(basename(motif_file), "\\.", 2)[,1]

  # Fix column names.
  colnames(dat) <- tolower(colnames(dat))
  colnames(dat) <- stringr::str_replace_all(colnames(dat), " ", "_")
  colnames(dat) <- stringr::str_replace(colnames(dat), "%", "_percent")
  colnames(dat) <- stringr::str_replace(colnames(dat), "/", "_over_")

  return(dat)
}

#' Read a HOMER HTML file.
#'
#' @export
read_findMotifs_html <- function(html_file) {
  doc <- XML::htmlParse(html_file)

  # Get the text inside h4 elements:
  o <- XML::getNodeSet(doc, "//h4")
  match_names <- sapply(o, xmlValue)

  o <- XML::getNodeSet(doc, "//table")
  x <- lapply(3:length(o), function(i) {
    tb <- o[[i]]
    tb <- tb[["tr"]][["td"]][["table"]]
    t(sapply(tb["tr"], function(x) sapply(x["td"], XML::xmlValue)))
  })

  x <- data.frame(do.call(rbind.data.frame, x[seq(1, length(x), by = 2)]))
  rownames(x) <- seq_len(nrow(x))
  colnames(x) <- c("variable", "value")

  x$variable <- stringr::str_replace(x$variable, ":", "")
  x$variable <- tolower(stringr::str_replace(x$variable, " ", "_"))
  x$match_name <- rep(match_names, each = 5)

  x <- as.data.frame(reshape2::dcast(x, match_name ~ variable, value.var = "value"))
  x$motif <- stringr::str_split_fixed(basename(html_file), "\\.", 2)[,1]

  x$match_rank <- as.numeric(x$match_rank)
  x$offset <- as.numeric(x$offset)
  x$score <- as.numeric(x$score)

  # Split the alignment into two parts
  midpoint <- nchar(x$alignment) / 2
  x$alignment1 <- substr(x$alignment, 1, midpoint)
  x$alignment2 <- substr(x$alignment, midpoint + 1, nchar(x$alignment))
  x$alignment <- NULL

  return(tibble::as_data_frame(x))
}
