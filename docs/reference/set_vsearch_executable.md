# Set the VSEARCH executable

`set_vsearch_executable` specifies the valid command to invoke
`VSEARCH`.

## Usage

``` r
set_vsearch_executable(vsearch_executable)
```

## Arguments

- vsearch_executable:

  (Required). Full path to the VSEARCH executable on your computer. See
  *Details* for more information on how to install `VSEARCH`.

## Value

Nothing is returned, but the option `Rsearch.vsearch_executable` is
updated. The string is also saved to a file for later R sessions, i.e.
you only need to update this once (or if you change how you run/install
VSEARCH).

## Details

Use this function to change the command used to invoke the external
software VSEARCH on this computer. When the `Rsearch` package is
installed this command is by default just `"vsearch"`.

If you have a windows computer and have copied the binary `vsearch.exe`
to the folder C:/Documents/ on your computer, you update R with this
information by `set_vsearch_executable("C:/Documents/vsearch")`.

You may use the function
[`vsearch`](https://cassandrahjo.github.io/Rsearch/reference/vsearch.md)
to test if the command is valid.

Visit <https://github.com/CassandraHjo/Rsearch> for more information on
how to install `VSEARCH`.

## See also

[`vsearch`](https://cassandrahjo.github.io/Rsearch/reference/vsearch.md).
