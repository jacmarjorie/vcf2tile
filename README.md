## VCF Importer for GenomicsDB

Basic script to produce callsets map and variant id map [required for GenomicsDB vcf import](https://github.com/Intel-HLS/GenomicsDB/wiki/Importing-VCF-data-into-GenomicsDB), and then perform the (non-mpi) loading process. This script was designed to statisfy Q1 requirements for loading VCF's from the Sanger VC workflow, though this supports arbitrary vcf types. The documentation below specifically focuses on how to load the VCF tar balls that are produced as the final output from the Sanger workflow. The vcf loading process will be replacd by the vcf importing process of the variant store in future quarters.  

**NOTE**: 

1. The [incremental callset addition process](https://github.com/Intel-HLS/GenomicsDB/wiki/Incremental-import-into-GenomicsDB), currently, requires that the callset map file be udpated before the actual loading into GenomicsDB as well as some values accurately specified in the loader config to identify where to begin loading. The vcf2tile importer takes this into account, however, only one VCF importing process can run at once or else multiple processes will attempt to write to the same file.
2. A VCF file from somatic variant calling contains two samples: a normal and a tumor sample. This means that there will be two callsets added to the callset map per somatic variant calling generated VCF file. 
3. The outputs of the Sanger VC Workflow are tar balls. This means the -t option is required to import these vcfs.

### General Flow Through of Import Process

example command for python script:

```
python vcf2tile.py -l /cluster_share/GenomicsDB/CCC_Omics/vcf2tile/array_info/example_loader.config -i /cluster_share/GenomicsDB/vcfs/DO221503.somatic.indel.tar.gz -s -t -A testing1 -L
```

example command for shell script:

```
./vcf2tile.sh -i="/path/to/DO221503.somatic.indel.tar.gz /path/to/DO221503.somatic.snv_mnv.tar.gz" -t -s -L -A=test
```

`-A` if exists, appends vcfs to array, else creates new array with the specified vcfs

`-i` can be a list of vcf files or a list of tars, it cannot be a list of both vcf and tars.

`-t` if `-i` is a list of tars

`-L` load tiledb array as defined by `-A` using the config files created in array_info

`-s` use SAMPLE tag in header to define the callsets

Use example loader file to create a new callset and vid map file, which will be associated by the array name (identified by `-A` flag). 

### General Setup Information

In order to load from a compute node GenomicsDB (git clone git@github.com:Intel-HLS/GenomicsDB.git) must exist on a shared location for all compute nodes in the rack (ie. /cluster_share on both PoC and sparkdmz). 

sparkdmz g1 `/cluster_share/GenomicsDB/` contains:

1. compiled GenomicsDB (`/cluster_share/GenomicsDB/GenomicsDB`)
2. GenomicsDB workspace where TileDB stores data arrays (`/cluster_share/GenomicsDB/WS`)
3. samtools/bcftools required for [collapsing and indexing vcf files](https://github.com/Intel-HLS/GenomicsDB/wiki/Useful-external-tools)
4. test vcf files (`/cluster_share/GenomicsDB/GenomicsDB/vcfs`) from the Sanger workflow 
5. `CCC_Omics/vcf2tile` contains python and shell wrapper script
6. `CCC_Omics/vcf2tile/array_info` directory with example configs, and a place where all array configs should live. The script is designed to always use the example configs and the `-A (--array)` flag to be able to load the correct array. For this to work correctly, the paths for workspace, callset mapping file, and vid mapping file need to be set correctly. 

### Quick Validation

To quicly validate that the array was loaded go to `array_info` directory and run:

where arrayname is the name of the array you wish to query

`./query arrayname arrayname_loader.config`


