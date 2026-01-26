# Log text to the current log file

The specified text is printed and written to the current log file. Does
not overwrite the file. Also used internally by many functions in the
package.

## Usage

``` r
log_text(text)
```

## Arguments

- text:

  The text to be logged

## Value

None, the function is invoked for its side effect.

## See also

[`init_log`](https://hanhineva-lab.github.io/notame/reference/init_log.md),
[`finish_log`](https://hanhineva-lab.github.io/notame/reference/finish_log.md)

## Examples

``` r
file_name <- "log.txt"
init_log(file_name)
#> INFO [2026-01-26 11:40:32] Starting logging
log_text("Hello World!")
#> INFO [2026-01-26 11:40:32] Hello World!
# Print the contents of the file
scan(file_name, sep = "\n", what = "character")
#> [1] "INFO [2026-01-26 11:40:32] Starting logging"
#> [2] "INFO [2026-01-26 11:40:32] Starting logging"
#> [3] "INFO [2026-01-26 11:40:32] Hello World!"    
```
