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

数据准备
本仓库当前包含以下本地数据集：

simulated_data/

realdata/

由于文件大小限制，真实数据集未直接包含在本仓库中。请在使用相应分析流程前，先自行下载数据集。

数据集下载
所有数据已通过百度网盘共享，具体信息如下：

文件夹名称：DiverTI_data
链接：https://pan.baidu.com/s/1vLHGqJsSIgpa59JZUw0EFQ?pwd=ej6v
下载后，您将获得以下两个子目录：

simulated_data/ —— 包含用于方法验证和测试的模拟数据集；

realdata/ —— 包含实际测序或实验获取的真实空间转录组数据

请将下载后的 simulated_data 和 realdata 文件夹放置于仓库根目录下（或按代码中指定的路径存放），即可运行对应的 Jupyter Notebook。
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
