# ORBIT paper: code and data

Code and processed data to reproduce the figures in:

> **Rank-based integration identifies convergent disease mechanisms across omics**

> Zifeng Qiu, Duncan Palmer, Luke Jostins, Andrew Lewis, Katherine Bull, Jagdeep Nanchahal, Yang Luo.

**macOS** —  Some figures are plotted with the `quartz()` graphics device, which is
  only available on macOS. If error appears, please replace `quartz(type = "pdf", file = ...)` with other functions. 

The ORBIT method itself is implemented as an R package:
**<https://github.com/yang-luo-lab/ORBIT>**.
This repository contains only the analysis scripts and data used in the
manuscript.

## Repository structure

```
.
├── Code/
│   ├── Figure 2.R
│   ├── Figure 3.R
│   ├── Figure 4.R
│   ├── Figure 5.R
│   ├── Sup Fig 1.R
│   ├── Sup Fig 2.R
│   ├── Sup Fig 7.R
│   ├── Sup Fig 8.R
│   └── Sup Fig 9.R
└── Data/
    ├── Benchmark and simulation/   # FPR / power results from simulations
    ├── CKD/                        # Chronic kidney disease tubulointerstitial (TI) omics
    ├── DCM/                        # Dilated cardiomyopathy transcriptome and proteome
    └── P_merge/                    # Method-comparison results for ORBIT-P (FPR / power)
```

### Script-to-figure mapping

| Script | Figure(s) | Data used |
|---|---|---|
| `Figure 2.R` | Fig. 2, Sup Fig. 3 | `Data/Benchmark and simulation/` |
| `Figure 3.R` | Fig. 3, Sup Fig. 4 | `Data/CKD/` |
| `Figure 4.R` | Fig. 4, Sup Fig. 5, Sup Fig. 6 | `Data/DCM/DCM_Transcriptome.rds` |
| `Figure 5.R` | Fig. 5 | `Data/DCM/` |
| `Sup Fig 1.R` | Supplementary Fig. 1 | NA |
| `Sup Fig 2.R` | Supplementary Fig. 2 | NA |
| `Sup Fig 7.R` | Supplementary Fig. 7 | NA |
| `Sup Fig 8.R` | Supplementary Fig. 8 | NA |
| `Sup Fig 9.R` | Supplementary Fig. 9 | `Data/P_merge/` |


## Data description

All data are processed / summary-level. No individual-level data are included.

**`Benchmark and simulation/`** — false-positive rate (FPR) and power results
from simulations varying the number of omics (`K`), the number of missing
features (`NA`), and between-omic correlation (`rho`).

**`CKD/`** — tubulointerstitial (TI) omics inputs (`TI_omics_inputs.rds`),
ORBIT results (`TI_orbit.rds`) and gene-set enrichment results (`TI_gsea.rds`).

**`DCM/`** — dilated cardiomyopathy transcriptome (`DCM_Transcriptome.rds`),
proteome (`DCM_Prot_final_dat.rds`) and three omics integration (`DCM_3omics_ORBIT_result.rds`).

**`P_merge/`** — FPR and power of ORBIT-P compared with other p-value combination
methods.

## Requirements

- R >= 4.5.2 
- ORBIT package (`remotes::install_github("yang-luo-lab/ORBIT")`)

## How to reproduce

1. Clone the repository:
   ```bash
   git clone https://github.com/yang-luo-lab/ORBIT-manuscript.git
   cd ORBIT-manuscript
   ```
2. Open R with the repository root as the working directory (or open the
   `.Rproj` file if provided).
3. Run any script in `Code/`. Each script is self-contained and reads from
   `Data/` using relative paths.


## Citation

If you use this code or data, please cite:


## License

Code is released under the MIT License. Data are released under


## Contact

Zifeng Qiu — qiuzifeng4518@gmail.com; zifeng.qiu@stx.ox.ac.uk
Yang Luo lab — <https://github.com/yang-luo-lab>
