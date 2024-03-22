<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
The `--bam` input parameter for this workflow accepts a path to a single BAM file, or a folder containing multiple BAM files for the sample. A sample name can be supplied with `--sample`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```