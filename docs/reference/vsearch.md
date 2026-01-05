# Test if `VSEARCH` can be executed

`vsearch` tests if the VSEARCH executable is a valid command.

## Usage

``` r
vsearch(verbose = TRUE)
```

## Arguments

- verbose:

  (Optional). Logical. Controls verbosity of the function. If `TRUE`,
  validation messages are printed to the console. Defaults to `TRUE`.

## Value

No return value, called for side effects (prints validation message to
console).

## Details

Use this function to test the command used to invoke the external
software VSEARCH on this computer.

## See also

[`set_vsearch_executable`](https://cassandrahjo.github.io/Rsearch/reference/set_vsearch_executable.md).
