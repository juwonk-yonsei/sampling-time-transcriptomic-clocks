# Clock-scoring packages

`scripts/confirm/score_clocks.py` and scripts `25`, `27`, `31`, `48` and `52` import two packages from the authors' unpublished repositories. The files they use are copied here unchanged, so the code that produced the results can be checked without those repositories. No other files of the two packages were used.

| Package | Version | Source repository and commit | Files | License |
|---|---|---|---|---|
| `clockspec` | 0.1.0 | clock-specificity `9d109499f0bf47a7869cd7b4bcffe311b8be9b75` | `__init__.py`, `clocks.py`, `preprocess.py`, `paths.py` | MIT (declared in the source `pyproject.toml`) |
| `scmodage` | 0.2.0 | sc-module-aging `a0e9d090cbabcdde441a599bb6443ed3f5ba0350` | `__init__.py`, `clocks.py`, `preprocess.py`, `paths.py` | MIT (`scmodage/LICENSE`) |

SHA-256 of the copied files (identical to the blobs at the commits above):

| File | SHA-256 |
|---|---|
| `clockspec/__init__.py` | `874a7f1fd586e813ad1818949b46b01577374224d46ceb902e006d0f350d4ad0` |
| `clockspec/clocks.py` | `9b0c05943d4f34bfbc35f700628ebd604cbb84477bed9b7a7c58703155c3d88c` |
| `clockspec/paths.py` | `4fb5236a63fac4daa5e64580852c1bd4ee252f8956e7bcaa38ce80faa1ec7c7c` |
| `clockspec/preprocess.py` | `0fa577daaa6c5fc9c7020d17d3b28e37885c7397f270a8a10a00fc466c175afd` |
| `scmodage/__init__.py` | `782372c5392dc964bd8f87a31a5b48d6598be4c6620d6558ff97e4a9407c24d3` |
| `scmodage/clocks.py` | `073f3bef2ac80455df32fd4118615178fa9bde57b3e72f7a647d96c328189877` |
| `scmodage/paths.py` | `cfb3bdfd0eba77271bc5a9bce49b7112eec7a2b1f937c2dad68dc1d1e65443ad` |
| `scmodage/preprocess.py` | `96384af955a5997107869dd120d32a1eec80aef47e1b195dedd5c47f745b01c7` |

The packages contain code only. Model files, gene lists, ortholog tables and module memberships are read at run time from the tAge resources, which are not redistributed (MGB Open Access License 1.0).

## Paths

- `scmodage/paths.py` reads tAge resources from `data/external/` two levels above itself, that is `<repository>/data/external/`.
- `clockspec/paths.py` reads them from the environment variable `CLOCKSPEC_TAGE`; its default is the analysis machine's path.

`scripts/setup_tage_layout.sh <tAge Zenodo folder> <tAge clone>` creates that folder with symbolic links:

| Link under `data/external/` | Target |
|---|---|
| `tAge` | clone of https://github.com/Gladyshev-Lab/tAge at `d048ba6016a7e5a114f45195e0b11071e2b8c424` (release 1.5.0; `inst/extdata` metadata) |
| `tage_small`, `tage_zenodo` | tAge Zenodo v2 folder (doi:10.5281/zenodo.22166800) |
| `tage_modules/rodent/Rodent module clocks` | `Rodent module clocks/` unpacked from the Zenodo archive |
| `tage_modules/multispecies/Multispecies module clocks` | `Multispecies module clocks/` unpacked from the Zenodo archive |

On the analysis machine these resources were read from `sc-module-aging/data/external/`. Its model files and unpacked module clocks are identical (SHA-256) to the Zenodo v2 files.

## Check

```bash
bash scripts/setup_tage_layout.sh /path/to/zenodo_22166800_v2 /path/to/tAge
PYTHON=/path/to/clock_env/bin/python bash scripts/verify_scoring.sh
```

`verify_scoring.sh` sets `CLOCKSPEC_TAGE` and `PYTHONPATH`, re-scores the committed GSE48113 (array) and GSE295758 (counts) gene matrices and compares the four outputs with `results/confirm/L1/` and `results/confirm/L3/`. On the analysis machine all four were byte-identical (about 20 s each).
