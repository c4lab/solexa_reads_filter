## solexa\_reads\_filter

solexa\_reads\_filter is a small command-line program to filter Solexa reads. Supported filters include:

* **s35**: From the 5' end of the read, the first 25 of 35 bases must have quality scores of at least 30, otherwise the read is discarded. (If the read length is lower than 35, the read will be discarded.)
* **Ns**: If there exists an ambiguous base (N call) then the read is discarded.
* **polyN**: If 85% of the read contains one type of base call (A, T, G, C), then the read deemed to have low complexity and is discarded.

## Usage

```
usage: solexa_reads_filter [-h] -r1 <file> -r2 <file> -Q (33|64) [-r <int>] -o
                           <dir> [-m] [-z] [-x] [-v] [--version]
                           
optional arguments:
  -h, --help  show this help message and exit
  -r1 <file>  read 1 in FASTQ format (also support gzip and bzip2 format)
  -r2 <file>  read 2 in FASTQ format (also support gzip and bzip2 format)
  -Q (33|64)  format of quality scores (33 or 64)
  -r <int>    minimum read length to be retained after trimming (default: 1)
  -o <dir>    output directory
  -m          merge read 1 and read 2 after filtering
  -z          turn off s35 filtering
  -x          turn off Ns filtering
  -v          turen off polyN filtering
  --version   show program's version number and exit
```

## License

solexa\_reads\_filter is released under the [MIT license](http://www.opensource.org/licenses/MIT).
