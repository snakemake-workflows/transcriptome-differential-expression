# Changelog

## [1.1.0](https://github.com/snakemake-workflows/transcriptome-differential-expression/compare/v1.0.0...v1.1.0) (2024-07-05)


### Features

* added wrappers ([#46](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/46)) ([9c3bfe4](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/9c3bfe49cde7118d841f95fcc65b36246e68d2c1))
* env update ([#45](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/45)) ([030e041](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/030e0411a261efad2178680997cf1b2e3de7e4d7))
* map qc ([#38](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/38)) ([4d5b8bc](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/4d5b8bcdd4fc98fc3c064f7f6d619cdf54b48e9c))
* qc rules ([#16](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/16)) ([3237e2c](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/3237e2c4e539d65bc49155ff18b27a4dc1df820a))
* read length filter ([#44](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/44)) ([10a8d93](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/10a8d93a21382ab4c8282e52d5a2539c06c394da))


### Bug Fixes

* added import sys for log ([#36](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/36)) ([0f8ac63](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/0f8ac63654ee6abb522d8590ae02f190be790c96))
* added new samples to config ([#39](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/39)) ([b993b01](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/b993b0168a63eed1ddbaf64f3aaf0e6d6b1061db))
* env update ([#42](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/42)) ([704b25c](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/704b25cb8a83a43d72b7ffe0775ab7081295720a))
* lookup [#18](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/18) ([#19](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/19)) ([ed226a3](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/ed226a3adc502694cc3ef9c30f3d2d5dae431935))
* read quantification  ([#22](https://github.com/snakemake-workflows/transcriptome-differential-expression/issues/22)) ([71422a1](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/71422a13319f0893e2494f432b0170bb116c08be))

## 1.0.0 (2024-04-30)


### Features

* added dummy 2nd contition to sample file ([812fc8a](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/812fc8a0a435f4b7cfd01d7cd3914c3cefc0de65))
* **customization:** mostly scriptizing for snakemake ([354a91c](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/354a91c70821816cfbff44e63e1aeae742390989))


### Bug Fixes

* **env:** fixed env for slurm and fs plugins ([65862e0](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/65862e0bace1d8e49a82e7c83bc597f178bde472))
* **param:** updated profile for human readable times ([a0395a2](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/a0395a2e67bdc7c4697f312e48170d9f8829d680))
* removed 'resdir' parameter from config file ([b28fdfe](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/b28fdfe0b72228302041e5ea1b990eb30d95b6ac))
* using configurable 'MAMBA_EXE' env var instead of a specific conda command ([a4b6e18](https://github.com/snakemake-workflows/transcriptome-differential-expression/commit/a4b6e18332a263cdd59edce28dd760656bba61bc))
