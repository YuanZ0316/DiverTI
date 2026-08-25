# DiverTI

**DiverTI** 是一个用于单细胞和空间转录组数据通用轨迹推断的 Python 框架。它整合了图特征建模、增量表达多样性分析和图注意力自编码器，能够重建细胞状态转换、发育轨迹和组织动态，呈现生物学方向性。

本仓库提供了安装说明、示例数据集以及展示完整分析流程的 Jupyter Notebook。

---

## 安装

### 1. 创建并激活 conda 环境

我们建议使用 Python 3.8（环境已预先配置好所有依赖）：

```bash
conda env create -f DiverTI_env.yml
conda activate DiverTI

依赖项
核心依赖（已在 DiverTI_env.yml 中锁定版本）：

Python 3.8

PyTorch 1.12+（推荐 GPU 版本）

Scanpy 1.9

CellRank 2.0

Scvelo 0.3

NumPy、Pandas、SciPy、Matplotlib、scikit‑learn 等

完整列表请参阅 DiverTI_env.yml 文件。

## 数据准备

本仓库目前包含以下本地数据集和 Jupyter Notebook 文件：
DiverTI/
├── simulated data/
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
