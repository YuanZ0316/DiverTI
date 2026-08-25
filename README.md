# DiverTI

**DiverTI: A Versatile Trajectory Inference Framework for Single-Cell and Spatial Transcriptomics**

DiverTI is a Python framework for versatile trajectory inference from single-cell and spatial transcriptomics data. It integrates graph feature modeling, incremental expression diversity analysis, and a graph attention autoencoder to reconstruct cell state transitions, developmental trajectories, and tissue dynamics, revealing biological directionality.

This repository provides installation instructions, example datasets, and Jupyter Notebooks demonstrating the complete analysis workflow.

---

## Installation

### 1. Create and activate a conda environment

We recommend using Python 3.8 (the environment comes with all dependencies pre-configured):

```bash
conda env create -f DiverTI_env.yml
conda activate DiverTI
```

### 2. Dependencies

Core dependencies (version-locked in `DiverTI_env.yml`):

- Python 3.8
- PyTorch 1.12+ (GPU version recommended)
- Scanpy 1.9
- CellRank 2.0
- Scvelo 0.3
- NumPy, Pandas, SciPy, Matplotlib, scikit-learn, etc.

See `DiverTI_env.yml` for the full list.

---

## Data Preparation

This repository currently contains the following local datasets (placeholder folders):

```
simulated_data/
realdata/
```

> **Note**: Due to file size limitations, the real datasets are not directly included in this repository. Please download the datasets yourself before running the corresponding analysis workflows.

### Dataset Download

All data is shared via Baidu Netdisk, as follows:

| Item | Content |
| --- | --- |
| Folder name | `DiverTI_data` |
| Link | https://pan.baidu.com/s/1vLHGqJsSIgpa59JZUw0EFQ?pwd=ej6v |
| Access code | `ej6v` |

After downloading, you will obtain the following two subdirectories:

- `simulated_data/` — contains simulated datasets for method validation and testing (including single-cell and spatial data);
- `realdata/` — contains real spatial transcriptomics data from actual sequencing or experiments (e.g., DLPFC, TBI model, ICC, etc.).

Place the downloaded `simulated_data` and `realdata` folders in the repository root directory (or the paths specified in the code) to run the corresponding Jupyter Notebooks.

---

## Repository Structure (Example)

```text
DiverTI/
├── simulated_data/
│   ├── scdata/
│   └── stdata/
├── realdata/
│   ├── DLPFC/
│   ├── real1/
│   ├── VLP4_C1_Visium/
│   └── zhang/
├── models/
├── util/
├── 01.DLPFC.ipynb
├── 02.TBI.ipynb
├── 03.real.ipynb
├── sc_linear1.ipynb
└── st_continuous.ipynb
```

---

## Usage Examples

The main example Notebook is `st_continuous.ipynb`, which demonstrates a complete spatial transcriptomics trajectory analysis workflow, including the following steps:

1. Load spatial transcriptomics data (supports both real and simulated data)
2. Data preprocessing (normalization, log transformation, highly variable gene selection)
3. Construct the cell spatial adjacency graph (based on PCA space or physical coordinates)
4. Construct the gene co-expression network (GeneNet)
5. Compute graph Fourier transform (GFT) features (gene domain and spatial domain)
6. Train the DiverTI model to obtain low-dimensional cell embeddings
7. Build an adaptive transition matrix from the embeddings to infer cell differentiation trajectories
8. Compute pseudotime and visualize the trajectory tree
9. Evaluate results against real developmental stage labels (if available)

In addition, you can use the following Notebooks for analysis of specific datasets:

| Notebook | Description |
| --- | --- |
| `sc_linear1.ipynb` | Trajectory inference example for single-cell data (linear differentiation) |
| `01.DLPFC.ipynb` | Hierarchical structure analysis of human dorsolateral prefrontal cortex (DLPFC) spatial transcriptomics data |
| `02.TBI.ipynb` | Spatial dynamic trajectories of the traumatic brain injury (TBI) regeneration process |
| `03.real.ipynb` | General analysis entry for other real single-cell or spatial transcriptomics data |

All Notebooks rely on the core functions in the `DiverTI/` module (utils, genenet, train, etc.) and support running on GPU or CPU. Simply organize your data into the `realdata/` and `simulated_data/` folders and modify the `DATA_PATH` and `ROOT_LABEL` parameters in the code to adapt to your own data.
