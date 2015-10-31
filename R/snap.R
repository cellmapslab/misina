
#slightly modified version of rsnps::LDSearch function
SNAP.LD <- function (SNPs, dataset = "onekgpilot", panel = "CEU", RSquaredLimit = 0.8,
                       distanceLimit = 500, GeneCruiser = TRUE, quiet = FALSE)
{
  library(RCurl)

  split_to_df <- function (x, sep, fixed = FALSE, perl = TRUE, useBytes = FALSE,
            names = NULL)
  {
    x <- as.character(x)
    if (fixed) {
      perl <- FALSE
    }
    tmp <- strsplit(x, sep, fixed = fixed, perl = perl, useBytes = useBytes)
    if (length(unique(unlist(lapply(tmp, length)))) > 1) {
      stop("non-equal lengths for each entry of x post-splitting")
    }
    tmp <- unlist(tmp)
    tmp <- as.data.frame(matrix(tmp, ncol = (length(tmp)/length(x)),
                                byrow = TRUE), stringsAsFactors = FALSE, optional = TRUE)
    if (!is.null(names)) {
      names(tmp) <- names
    }
    else {
      names(tmp) <- paste("V", 1:ncol(tmp), sep = "")
    }
    return(tmp)
  }

  tmp <- sapply(SNPs, function(x) {
    grep("^rs[0-9]+$", x)
  })
  if (any(sapply(tmp, length) == 0)) {
    stop("not all items supplied are prefixed with 'rs';\n",
         "you must supply rs numbers and they should be prefixed with ",
         "'rs', e.g. rs420358")
  }
  if (RSquaredLimit < 0 || RSquaredLimit > 1) {
    stop("RSquaredLimit must be between 0 and 1")
  }
  if (is.character(distanceLimit)) {
    n <- nchar(distanceLimit)
    stopifnot(substring(distanceLimit, n - 1, n) == "kb")
    distanceLimit <- as.integer(gsub("kb", "", distanceLimit))
  }
  valid_distances <- c(0, 10, 25, 50, 100, 250, 500)
  if (!(distanceLimit %in% valid_distances)) {
    stop("invalid distanceLimit. distanceLimit must be one of: ",
         paste(valid_distances, collapse = ", "))
  }
  distanceLimit_bp <- as.integer(distanceLimit * 1000)
  query_start <- "http://www.broadinstitute.org/mpg/snap/ldsearch.php?"
  SNP_query <- paste(sep = "", "snpList=", paste(SNPs, collapse = ","))
  dataset_query <- paste(sep = "", "hapMapRelease=", dataset)
  panel_query <- paste(sep = "", "hapMapPanel=", panel)
  RSquaredLimit_query <- paste(sep = "", "RSquaredLimit=",
                               RSquaredLimit)
  distanceLimit_query <- paste(sep = "", "distanceLimit=",
                               distanceLimit_bp)
  downloadType_query <- paste(sep = "", "downloadType=file")
  if (GeneCruiser) {
    columnList_query <- paste(sep = "", "columnList[]=DP,GA,MAF")
  }
  else {
    columnList_query <- paste(sep = "", "columnList[]=DP,MAF")
  }
  includeQuerySNP_query <- "includeQuerySnp=on"
  submit_query <- paste(sep = "", "submit=search")
  query_end <- paste(sep = "&", SNP_query, dataset_query,
                     panel_query, RSquaredLimit_query, distanceLimit_query,
                     downloadType_query, columnList_query, includeQuerySNP_query,
                     submit_query)
  query <- paste(sep = "", query_start, query_end)
  if (!quiet)
    cat("Querying SNAP...\n")
  dat <- getURL(query)
  if (length(grep("validation error", dat)) > 0) {
    stop(dat)
  }
  tmp <- unlist(strsplit(dat, "\r\n", fixed = TRUE))
  warning_SNPs <- grep("WARNING", tmp, value = TRUE)
  for (line in warning_SNPs) {
    warning(line)
  }
  bad_lines <- grep("WARNING", tmp)
  if (length(bad_lines) > 0) {
    tmp <- tmp[-bad_lines]
  }
  out <- split_to_df(tmp, sep = "\t", fixed = TRUE)
  names(out) <- unlist(unclass(out[1, ]))
  out <- out[2:nrow(out), ]
  return(out)
}
