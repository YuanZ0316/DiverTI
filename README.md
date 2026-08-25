# DiverTI

**DiverTI: 通用轨迹推断框架（单细胞与空间转录组）**

DiverTI 是一个用于单细胞和空间转录组数据通用轨迹推断的 Python 框架。它整合了图特征建模、增量表达多样性分析和图注意力自编码器，能够重建细胞状态转换、发育轨迹和组织动态，呈现生物学方向性。

本仓库提供了安装说明、示例数据集以及展示完整分析流程的 Jupyter Notebook。

---

## 安装

### 1. 创建并激活 conda 环境

我们建议使用 Python 3.8（环境已预先配置好所有依赖）：

```bash
conda env create -f DiverTI_env.yml
conda activate DiverTI
```

### 2. 依赖项

核心依赖（已在 `DiverTI_env.yml` 中锁定版本）：

- Python 3.8
- PyTorch 1.12+（推荐 GPU 版本）
- Scanpy 1.9
- CellRank 2.0
- Scvelo 0.3
- NumPy、Pandas、SciPy、Matplotlib、scikit-learn 等

完整列表请参阅 `DiverTI_env.yml` 文件。

---

## 数据准备

本仓库当前包含以下本地数据集（占位文件夹）：

```
simulated_data/
realdata/
```

> **注意**：由于文件大小限制，真实数据集未直接包含在本仓库中。请在使用相应分析流程前，先自行下载数据集。

### 数据集下载

所有数据已通过百度网盘共享，具体信息如下：

| 项目 | 内容 |
| --- | --- |
| 文件夹名称 | `DiverTI_data` |
| 链接 | https://pan.baidu.com/s/1vLHGqJsSIgpa59JZUw0EFQ?pwd=ej6v |
| 提取码 | `ej6v` |

下载后，您将获得以下两个子目录：

- `simulated_data/` —— 包含用于方法验证和测试的模拟数据集（含单细胞和空间数据）；
- `realdata/` —— 包含实际测序或实验获取的真实空间转录组数据（如 DLPFC、TBI 模型、ICC 等）。

请将下载后的 `simulated_data` 和 `realdata` 文件夹放置于仓库根目录下（或按代码中指定的路径存放），即可运行对应的 Jupyter Notebook。

---

## 仓库目录结构（示例）

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

## 使用示例

主要的示例 Notebook 为 `st_continuous.ipynb`，该文件演示了一套完整的空间转录组轨迹分析流程，具体步骤包括：

1. 加载空间转录组数据（支持真实数据或模拟数据）
2. 数据预处理（归一化、对数变换、筛选高变基因）
3. 构建细胞空间邻接图（基于 PCA 空间或物理坐标）
4. 构建基因共表达网络（GeneNet）
5. 计算图傅里叶变换（GFT）特征（基因域和空间域）
6. 训练 DiverTI 模型，获取细胞低维嵌入表示
7. 基于嵌入构建自适应转移矩阵，推断细胞分化轨迹
8. 计算拟时序（pseudotime）并可视化轨迹树
9. 对比真实发育阶段标签（若可用）评估结果

此外，您还可以使用以下其他 Notebook 进行特定数据集的分析：

| Notebook | 用途 |
| --- | --- |
| `sc_linear1.ipynb` | 单细胞数据的轨迹推断示例（线性分化） |
| `01.DLPFC.ipynb` | 人类背外侧前额叶皮层（DLPFC）空间转录组数据的层次结构分析 |
| `02.TBI.ipynb` | 创伤性脑损伤（TBI）再生过程的空间动态轨迹 |
| `03.real.ipynb` | 其他真实单细胞或空间转录组数据的通用分析入口 |

所有 Notebook 均依赖 `DiverTI/` 模块中的核心函数（utils、genenet、train 等），并支持在 GPU 或 CPU 环境下运行。您只需将数据按 `realdata/` 和 `simulated_data/` 文件夹组织好，并修改代码中的 `DATA_PATH` 和 `ROOT_LABEL` 参数即可适配自己的数据。
