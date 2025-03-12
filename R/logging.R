#' Initialize log to a file
#'
#' Initialize a log file with the current data and time.
#' All major operations run after this will be logged to the specified file.
#'
#' @param log_file Path to the log file
#'
#' @return None, the function is invoked for its side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' file_name <- "log.txt"
#' init_log(file_name)
#' # Print the contents of the file
#' scan(file_name, sep = "\n", what = "character")
#' \dontshow{setwd(.old_wd)}
#' @seealso \code{\link{log_text}}, \code{\link{finish_log}}
#'
#' @export
init_log <- function(log_file) {
  futile.logger::flog.appender(futile.logger::appender.tee(log_file), 
                               name = "notame")
  log_text("Starting logging")
  # Pass errors to log
  options(error = function() {
    futile.logger::flog.error(geterrmessage(), name = "notame")
  })
}

#' Log text to the current log file
#'
#' The specified text is printed and written to the current log file. Does not 
#' overwrite the file. Also used internally by many functions in the package. 
#'
#' @param text The text to be logged
#'
#' @return None, the function is invoked for its side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' file_name <- "log.txt"
#' init_log(file_name)
#' log_text("Hello World!")
#' # Print the contents of the file
#' scan(file_name, sep = "\n", what = "character")
#' \dontshow{setwd(.old_wd)}
#'
#' @seealso \code{\link{init_log}}, \code{\link{finish_log}}
#'
#' @export
log_text <- function(text) {
  futile.logger::flog.info(text, name = "notame")
}

#' Finish a log
#'
#' Logs the current date and time and session info, and switches logging off.
#'
#' @return None, the function is invoked for its side effect.
#'
#' @examples
#' finish_log()
#'
#' @seealso \code{\link{init_log}}, \code{\link{log_text}}
#'
#' @export
finish_log <- function() {
  # Return default option for error
  options(error = NULL)
  # Log end of session info
  futile.logger::flog.info(paste("Finished analysis. ", date(), 
                                 "\nSession info:\n", sep = ""))
  futile.logger::flog.info(utils::capture.output(utils::sessionInfo()))
  invisible(futile.logger::flog.appender(futile.logger::appender.console(),
                                         name = "notame"))
}
