# pre_assembler
Visualize de bruijn assemble subgraph.

## keys
`-i, --input` Path to input .fasta  
`-k, --kmer` k-mer length (default=3)  
`-t, --type` type of output graph (full/not) (default=full)  
`-s, --strand` reverse (bw) or not (fw) (default=fw)    

## run example

### example1:
```
python3 assembler.py -i hw3_dataset11.fasta
```
### example2:
```
python3 assembler.py -i hw3_dataset.fasta -k 55 -t not -s bw
```
