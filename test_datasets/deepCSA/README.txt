# deepCSA Test Datasets

This folder contains test data to be used for automated testing with the
bbglab/deepCSA pipeline (https://github.com/bbglab/deepCSA).

---

## Content of this repository

`deepCSA/testdata/maf/` : Somatic mutation calls (MAF format) for all samples.
`deepCSA/testdata/depth/` : Per-sample sequencing depth table.

Full tree:

    .
    └── deepCSA
        └── testdata
            ├── depth
            │   └── all_samples_indv.depths.tsv.gz
            └── maf
                └── all_samples.somatic.mutations.maf

---

## Dataset information

This is a small subset (3 samples) of the normal bladder urothelium dataset
published in:

    Calvet, F., Blanco Martinez-Illescas, R. et al.
    "Sex and smoking bias in the selection of somatic mutations in human bladder."
    Nature (2025). https://doi.org/10.1038/s41586-025-09521-x

The full dataset is publicly available at Zenodo:

    https://zenodo.org/records/15836679
    DOI: 10.5281/zenodo.15836679

### Samples included

| Sample ID       | Description              |
|-----------------|--------------------------|
| P19_0002_BDO_01 | Normal bladder urothelium |
| P19_0002_BTR_01 | Normal bladder urothelium |
| P19_0003_BDO_01 | Normal bladder urothelium |

### How to obtain the test files from Zenodo

Download and extract the Zenodo archive:

    wget https://zenodo.org/records/15836679/files/normal_bladder_urothelium.tar.gz
    tar -xzf normal_bladder_urothelium.tar.gz

The three samples in this test dataset are a subset of the full cohort.


---

## Reference files

Running deepCSA requires several external reference files. These are NOT
included in this repository due to size constraints. Below is a list of the
required files and where to obtain them, following the documentation at:

    https://github.com/bbglab/deepCSA/blob/main/docs/usage.md


## Usage in the pipeline (test configuration)

The test uses pre-computed mutation calls and depth tables instead of running
from raw BAM/VCF files. The relevant test parameters are:

    input_maf             = "path/to/deepCSA/testdata/maf/all_samples.somatic.mutations.maf"
    use_custom_depths     = true
    custom_depths_table   = "path/to/deepCSA/testdata/depth/all_samples_indv.depths.tsv.gz"
    use_custom_minimum_depth = 10

See the pipeline test suite documentation for more details:

    https://github.com/bbglab/deepCSA/blob/main/tests/README.md
    https://github.com/bbglab/deepCSA/pull/438
