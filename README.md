# covbamic

### requirements

This tool has been tested with

python >= 3.9
pysam >= 0.17.0


### Installation

```git clone https://github.com/mjsull/covbamic.git```


### Usage

```python covbamic/covbamic.py -b sars_bam_file.primertrimmed.rg.sorted.bam -o output.svg -1 BA.4 -2 BA.5```

Where sars_bam_file.primertrimmed.rg.sorted.bam is a sorted and indexed bam file.

### options


#### required
```
  -o OUTPUT, --output OUTPUT
                        output svg file
  -b BAM_FILE, --bam_file BAM_FILE
                        sorted and indexed bam file
  -1 VARIANT_1, --variant_1 VARIANT_1
                        variant 1 (BA.2, BA.4 or BA.5)
  -2 VARIANT_2, --variant_2 VARIANT_2
                        variant 2 (BA.2, BA.4 or BA.5)


```

#### optional

```
  -h, --help            show this help message and exit

  -a, --all             List all sites different from reference (as opposed to only sites that differ between the two variants selected).
  -m, --all_minor       List all sites where the minor allele reaches defined threshold.
  -f MINOR_FRACTION, --minor_fraction MINOR_FRACTION
                        Fraction of reads with minor allele to report (when -m set)
  -d MINOR_DEPTH, --minor_depth MINOR_DEPTH
                        minimum depth to report minor allele site (when -m set)
  -p3, --panel3         Draw panel 3.
```


Example output

![covbamic](https://github.com/mjsull/covbamic/blob/main/example.svg?raw=true)