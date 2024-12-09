# HiRes
Haemophilus influenzae Resistance Gene Extraction Tool

Takes reads and assemblies to determine if Haemophilus influenzae is resistant to either betalactams, fluoroquinolones or cephalosporins

```
pertpipe --R1 $PATH/$R1.fq.gz --R2 $PATH/$R2.fq.gz --outdir $OUTDIR
```

FLAGS

```
--outdir, -o [PATH]             optional folder to write output files to
--R1 [PATH]                     R1 fastq of sample (can be gzipped files)
--R2 [PATH]                     R2 fastq of sample (can be gzipped files)
--fasta, -f [PATH]              optional fasta file which will skip spades.             
--version, -v                   print version
```

## Dependencies
- SPades
- Snippy
- Abricate
