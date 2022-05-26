# HiFi assembly tools - *Leiopelma hochstetteri* genome assembly

Scripts associated with genome assembly using PacBio HiFi data.

[asm](https://github.com/natforsdick/hifi-assembly-tools/tree/main/asm) contains scripts used to produce draft assemblies from HiFi data.
These use the [hifiasm](https://github.com/chhylp123/hifiasm) and [HiCanu](https://canu.readthedocs.io/en/latest/quick-start.html#assembling-pacbio-hifi-with-hicanu) assemblers.

[purge_dups](https://github.com/natforsdick/hifi-assembly-tools/tree/main/purge_dups) implements the [purge_dups](https://github.com/dfguan/purge_dups) pipeline to remove haplotig and contig overlaps from the draft assembly prior to scaffolding. This pipeline is based on that implemented by Sarah Bailey (UoA) with modification.


