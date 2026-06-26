# FormFactors 中文说明

本仓库用于核子电磁形状因子（nucleon electromagnetic form factors）的全局拟合、误差估计、结果汇总与后处理绘图。当前代码主要围绕质子和中子的 Sachs form factors：

- 质子电形状因子 `G_E^p`
- 质子磁形状因子 `G_M^p`
- 中子电形状因子 `G_E^n`
- 中子磁形状因子 `G_M^n`

两套拟合代码分别位于 `final_fit_proton/` 和 `final_fit_neutron/`，最终参数化结果、SBS projection、`F_1/F_2` 与 flavor-separated quark form factors 的整理和绘图主要位于 `ff_three_sets/`。`charge_density/` 是后续电荷密度相关分析目录，目前暂不作为本 README 的主要工作流。

## 目录结构

```text
FormFactors/
├── README.md
├── README_Ch.md
├── charge_density/
│   └── charge_distribution.ipynb
├── final_fit_proton/
│   ├── fitff
│   ├── FF_sumrules_leastsq.py
│   ├── FF_sumrules_coverr.py
│   ├── FF_funcs.py
│   ├── FF_loaddata.py
│   ├── run_coverr.py
│   ├── plot_com.py
│   ├── data/
│   └── */z*/   # 不同日期、数据集和拟合选项的历史输出
├── final_fit_neutron/
│   ├── fitff
│   ├── neut_FF_sumrules_leastsq.py
│   ├── neut_FF_sumrules_coverr.py
│   ├── FF_funcs.py
│   ├── FF_loaddata.py
│   ├── neut_run_coverr.py
│   ├── neut_plot_FF.py
│   ├── data/
│   └── */z*/   # 不同日期、数据集和拟合选项的历史输出
└── ff_three_sets/
    ├── GetFF.py
    ├── GetFF_SBSproj.py
    ├── GetFF_final_SBS.py
    ├── FF_newGMn_SBSproj.py
    ├── FF_newGp_SBSproj.py
    ├── print_gegm_table.py
    ├── plot_new.ipynb
    ├── SBSresult_Apr2026/
    │   ├── proton/
    │   │   └── SBS_fit_Gp_t0fix7.dat
    │   ├── neutron/
    │   │   └── out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat
    │   └── plots/
    ├── plotdata/
    │   ├── GMn_lookup_SBSproj.dat
    │   └── final_fit_SBS/
    │       ├── SBS_fit_Gp_t0fix7.dat
    │       ├── out_all_sumrules_leastsq_Q21000_z12_gb5_t0fix7_SBS.dat
    │       └── out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat
    ├── data/
    └── Quarks_* / Proton_* / Neutron_* 输出图
```

## 物理与数值方法

本项目采用的核心方法是 bounded z expansion。形状因子写为

```text
G(Q^2) = sum_k a_k z(Q^2)^k
```

其中

```text
z = (sqrt(t_cut + Q^2) - sqrt(t_cut - t0))
    / (sqrt(t_cut + Q^2) + sqrt(t_cut - t0))
```

代码中常见设置为：

- `tcut = 4 * mpi**2`
- `t0 = -0.7`，对应文件名中的 `_t0fix7`
- 质子默认 `kmax = 12`
- 中子当前主入口默认 `kmax = 11`
- 系数 bound 通常为 `bound = 5`

拟合中加入 sum-rule constraints，用于同时满足 `Q^2 = 0` 处归一化和大 `Q^2` 渐进行为。代码里通常使用 `nmax = kmax - 4`，最后若干 z-expansion 系数由 sum rules 反解得到，而不是全部自由拟合。

Sachs form factors 与 Dirac/Pauli form factors 的关系按项目 notes 中的约定为：

```text
G_E = F_1 - tau * kappa * F_2
G_M = F_1 + kappa * F_2
tau = Q^2 / (4 M_N^2)
```

在后续电荷密度分析中需要注意：在相对论体系下，不能简单把 `G_E(Q^2)` 或 `G_M(Q^2)` 当成三维静态电荷/磁化密度的 Fourier transform。更稳妥的解释通常使用 light-front 或 infinite-momentum-frame 下的二维 transverse density。

## 质子拟合：`final_fit_proton/`

`final_fit_proton/` 是 `G_E^p` 和 `G_M^p` 全局联合拟合的主目录。

主要文件：

- `fitff`：运行质子拟合的 shell wrapper。
- `FF_sumrules_leastsq.py`：主拟合脚本，使用 `scipy.optimize.leastsq`。
- `FF_sumrules_coverr.py`：读取 central fit 结果并生成误差/协方差传播后的输出。
- `run_coverr.py`：批量调用 `FF_sumrules_coverr.py`。
- `FF_loaddata.py`：读取 Mainz、world、polarization、fake high-Q2、SBS 等数据。
- `FF_funcs.py`：形状因子、截面、z expansion、TPE/radiative correction 等基础函数。
- `plot_com.py` / `plot_com.ipynb`：比较不同拟合版本并画图。

### 运行示例

请在 `final_fit_proton/` 目录下运行，因为数据文件路径多为相对路径：

```bash
cd final_fit_proton
bash fitff 1 1
```

`fitff` 的第一个参数选择数据集：

```text
1 -> all
2 -> world
3 -> Mainz
```

第二个参数选择拟合选项：

```text
1 -> default
2 -> no high-Q2 fake constraints，输出名含 _noHQ
3 -> no radius constraints，输出名含 _noRad
4 -> no extra high-Q2 TPE uncertainty，输出名含 _noExtraTPE
5 -> fixed radius option，输出名含 _RE8409_RM851
```

当前 `fitff` 中的重要默认值：

```text
com    = Apr2026
kmax   = 12
Q2max  = 1000
bound  = 5
mod    = _t0fix7
method = 1
```

输出目录形如：

```text
final_fit_proton/Apr2026_all_FF/z12/
```

拟合输出文件形如：

```text
all_sumrules_leastsq_Q21000_z12_gb5_t0fix7.dat
all_sumrules_leastsq_Q21000_z12_gb5_t0fix7.txt
```

注意：`FF_sumrules_leastsq.py` 会根据输出文件名中的关键词判断数据集和选项。例如 `world`、`all`、`Mainz` 决定 `csopt`，`_noHQ`、`_noRad`、`_noExtraTPE` 决定是否使用 high-Q2 fake points、半径约束和额外 TPE 项。因此修改输出文件命名时要非常小心。

### 误差输出

central fit 完成后，可运行：

```bash
cd final_fit_proton
python3 run_coverr.py 1 1
```

其中第一个参数选择数据集，第二个参数选择输出类型。`run_coverr.py` 会从 central fit `.dat` 生成 `out_*` 文件，供后处理和绘图使用。

## 中子拟合：`final_fit_neutron/`

`final_fit_neutron/` 是 `G_E^n` 和 `G_M^n` 拟合目录。与质子不同，中子代码中 `G_E^n` 和 `G_M^n` 通常分开拟合。

主要文件：

- `fitff`：运行中子拟合的 shell wrapper。
- `neut_FF_sumrules_leastsq.py`：主拟合脚本。
- `neut_FF_sumrules_coverr.py`：中子误差/协方差传播脚本。
- `neut_run_coverr.py`：批量调用 `neut_FF_sumrules_coverr.py`。
- `FF_loaddata.py`：读取 neutron world data、fake high-Q2 points 和 SBS `G_E^n` 数据。
- `FF_funcs.py`：中子拟合使用的基础函数。
- `neut_plot_FF.py` / `neut_plot_FF.ipynb`：中子形状因子绘图。

### 运行示例

请在 `final_fit_neutron/` 目录下运行：

```bash
cd final_fit_neutron
bash fitff 1   # GEn
bash fitff 2   # GMn
```

`fitff` 的第一个参数：

```text
1 -> GEn
2 -> GMn
```

当前 `fitff` 中的重要默认值：

```text
fitdata = Apr2026_world
kmax    = 11
Q2_Max  = 1000
bound   = 5
fmod    = _t0fix7
method  = 1
```

输出目录形如：

```text
final_fit_neutron/Apr2026_world_norm_bound5/z11/
```

输出文件形如：

```text
Apr2026_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn.dat
Apr2026_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GMn.dat
```

中子误差输出示例：

```bash
cd final_fit_neutron
python3 neut_run_coverr.py 1 1   # GEn default
python3 neut_run_coverr.py 2 1   # GMn default
```

## 汇总与后处理：`ff_three_sets/`

`ff_three_sets/` 汇总最终 form factor 参数化、SBS projection、`F_1/F_2`、flavor-separated quark form factors 和论文/报告级别图像。当前应优先把 `SBSresult_Apr2026/` 视为最新 Apr2026 SBS 拟合结果的归档目录。

主要文件：

- `GetFF.py`：基础 form factor 参数化接口。
- `GetFF_SBSproj.py`：SBS projection 版本的参数化接口，主要用于包含 Tyler `G_M^n` 更新的中间结果。
- `GetFF_final_SBS.py`：最终 SBS 更新版本，包含 Haocen 的 2026 proton/GEn SBS fit 和 Tyler 更新的 `G_M^n`。
- `FF_newGMn_SBSproj.py`：Tyler 更新 `G_M^n` 的 SBS projection 接口。
- `FF_newGp_SBSproj.py`：更新 proton form factors 的 SBS projection 接口。
- `plot_new.ipynb`：后处理 notebook；最新正式结果以 `SBSresult_Apr2026/` 中的 lookup table 和 plots 为准。
- `print_gegm_table.py`：调用 `GetFF(kID, Q2)` 输出 `mu_p G_E^p/G_M^p` 和 `mu_n G_E^n/G_M^n` 表格。

`GetFF(kID, Q2)` 的 `kID` 约定：

```text
1 -> GEp
2 -> GMp
3 -> GEn
4 -> GMn
```

示例：

```bash
cd ff_three_sets
python3 print_gegm_table.py
```

### Apr2026 SBS 结果

```text
ff_three_sets/SBSresult_Apr2026/
├── proton/SBS_fit_Gp_t0fix7.dat
├── neutron/out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat
└── plots/
```

其中两个 `.dat` 文件是最新 Apr2026 拟合结果的 lookup table：

- `proton/SBS_fit_Gp_t0fix7.dat`：proton FF lookup table，包含 `Q2`、`z`、`G_E^p`、`G_E^p/G_D`、`G_M^p`、`G_M^p/(mu_p G_D)`、`mu_p G_E^p/G_M^p` 及其误差。
- `neutron/out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat`：neutron `G_E^n` lookup table，包含 `Q2`、`z`、`G_E^n`、`dG_E^n` 和 `G_D`。

`plots/` 目录保存 Apr2026 拟合对应的检查图和汇总图。后续 charge density / transverse density 分析应优先读取 `SBSresult_Apr2026/` 中的 lookup table，而不是旧的临时输出目录。

该目录中的 `Quarks_*`、`Proton_*`、`Neutron_*` 文件是最终绘图输出，包含 Sachs form factors、Dirac/Pauli form factors、SBS projection 和 up/down quark flavor-separated 结果。

## 数据流概览

典型工作流可以理解为：

```text
data/proton, data/neutron
        |
        v
final_fit_proton/FF_sumrules_leastsq.py
final_fit_neutron/neut_FF_sumrules_leastsq.py
        |
        v
central fit .dat / .txt
        |
        v
FF_sumrules_coverr.py / neut_FF_sumrules_coverr.py
        |
        v
out_* form factor tables and uncertainties
        |
        v
ff_three_sets/GetFF*.py and plot scripts
        |
        v
SBSresult_Apr2026 lookup tables and plots
        |
        v
latest Apr2026 fit products for later density analysis
```

## 依赖

仓库没有固定的 environment 文件。根据当前脚本，常用依赖包括：

- Python 3
- `numpy`
- `scipy`
- `matplotlib`
- Jupyter Notebook，用于 `.ipynb`

部分脚本会读取相对路径下的数据文件，因此推荐在对应子目录内运行脚本，而不是从仓库根目录直接运行。

## 重要注意事项

- 不要随意删除历史输出目录。诸如 `mar30_*`、`may23_*`、`jun14_*`、`oct03_*`、`Apr2026_*` 的目录记录了不同日期、数据集、`kmax`、`t0`、约束条件和系统误差处理方式下的结果。

- 不要随意覆盖已有 `.dat` 输出。主拟合脚本使用写模式输出，重新运行同名文件会覆盖旧结果。

- `ff_three_sets/SBSresult_Apr2026/proton/SBS_fit_Gp_t0fix7.dat` 和 `ff_three_sets/SBSresult_Apr2026/neutron/out_world_sumrules_leastsq_Q21000_z11_gb5_t0fix7_GEn_SBS.dat` 是最新 Apr2026 拟合结果的 lookup table，后续密度分析会继续使用。

  
