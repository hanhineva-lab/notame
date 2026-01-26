# Initialize log to a file

Initialize a log file with the current data and time. All major
operations run after this will be logged to the specified file.

## Usage

``` r
init_log(log_file)
```

## Arguments

- log_file:

  Path to the log file

## Value

None, the function is invoked for its side effect.

## See also

[`log_text`](https://hanhineva-lab.github.io/notame/reference/log_text.md),
[`finish_log`](https://hanhineva-lab.github.io/notame/reference/finish_log.md)

## Examples

``` r
file_name <- "log.txt"
init_log(file_name)
#> INFO [2026-01-26 11:40:32] Starting logging
# Print the contents of the file
scan(file_name, sep = "\n", what = "character")
#> [1] "INFO [2026-01-26 11:40:32] Starting logging"
```
