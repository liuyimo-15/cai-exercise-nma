# 慢性踝关节不稳运动康复：六结局网状Meta分析复现项目

这是基于本研究实际分析数据和已运行脚本整理的 GitHub 项目。包含传统治疗节点 NMA、固定成分及交互规格的 CNMA、逐一剔除与高风险剔除敏感性分析、亚组分析、单变量网状Meta回归，以及可选的条件性 CINeMA/GRADE 复算。

**主分析代码的来源：目前未取得最初生成 R050 PDF 的原始脚本或模型对象。本项目保留的是此前根据 R050 森林图、异质性与 SIDE 数值核对形成的重建代码。它可以重新计算并核对已匹配的结果，但不代表已经恢复原始 R050 程序、筛选理由或全部60页图形。** 原始结果 PDF 和原始固定数值 Excel 均原样保留。

English instructions: [README_en.md](README_en.md).

## 快速运行

1. 解压项目，在 RStudio 中打开 `cai-exercise-nma.Rproj`，或在终端进入项目根目录。
2. 建议使用记录分析环境的 **R 4.4.0**。首次运行安装锁定依赖，需要网络连接；从源码安装部分依赖时还需要对应系统的编译工具。
3. 运行：

```r
# 在终端执行
Rscript setup.R
Rscript run_all.R
```

在 RStudio 控制台也可执行：

```r
source("setup.R")
source("run_all.R")
```

`setup.R` 使用 `renv.lock` 恢复到项目的 `.Rlibrary` 文件夹。锁文件依据本机实际安装包的版本及依赖元数据整理；未把本机编译的R包二进制文件打包。主要依赖为 netmeta 3.3-1、metafor 4.8-0、stringi 1.8.7。完整记录见 `dependencies/package_versions.csv`。`renv` 的用法见其[官方恢复说明](https://rstudio.github.io/renv/reference/restore.html)。

入口自动执行：输入检查 → 六结局主模型重建 → 网络诊断 → NMA/CNMA敏感性 → 亚组 → 回归 → CR1与代数核查 → 图形 → CSV导出 → 与冻结结果逐表比较。输入CSV会复制至工作目录；所有模型RDS均重新生成，入口不读取归档RDS或冻结结果来代替拟合。`expected/` 只用于计算完成后的核查。

```r
# 只省略图形，仍拟合全部模型
Rscript run_all.R --no-figures

# 主动允许不同依赖版本，并仍进行数值核查
Rscript run_all.R --allow-version-mismatch
```

如需在已有R包库运行，可设置 `NMA_R_LIBRARY` 为该库路径。该选项不替代版本检查。首次依赖恢复、Windows/Linux完整运行和GitHub Actions实际云端运行的验证状态，见 `validation/REPRODUCIBILITY_REPORT.md`，不能以本机通过推断所有平台均通过。

## 数据与结果在哪里

| 目录或文件 | 用途 |
|---|---|
| `data/source_snapshots/` | 原固定数值Excel、R050参考PDF、所用RoB2情景PDF；均保留原文件字节 |
| `data/analysis_inputs/` | 本次实际拟合使用的CSV、CNMA规格和阈值JSON；运行起点 |
| `data/extraction/` | 逐臂人口学及研究映射、协变量、来源页码、缺失与核查说明 |
| `data/study_bibliography.csv` | 报告ID、作者、标题、DOI及原文件名；不包含期刊全文 |
| `work/` | 实际R脚本及运行中间文件；保留旧目录名便于核对依赖 |
| `R/` | 新增的相对路径入口辅助、输入核查、主结果绘图、导出和数值核对 |
| `results/tables/` | 本次从CSV重新拟合产生的完整数值表 |
| `results/figures/` | 本次重新绘制的矢量PDF |
| `results/plot_data/` | 主图、敏感性、亚组、回归图对应的数值数据 |
| `results/model_objects/` | 本轮重算的RDS模型对象，便于检查；入口仍会从CSV重新拟合 |
| `results/previews/` | PNG预览图；矢量PDF用于放大和排版 |
| `results/validation/` | 依赖版本、输入一致性、27张结果表、R050显示精度的核对记录 |
| `reference_outputs/` | 此前交付的Excel和原版图形，仅归档，不由入口重建Excel版式 |
| `expected/` | 原分析的冻结数值，仅用于独立比较，不作为模型输入 |
| `archive/original_R_scripts/` | 整理前的R脚本原样副本；**不要从此目录直接运行** |
| `docs/` | 数据字典、方法、已知限制、脚本调整记录 |
| `SHA256SUMS.csv` | 发布文件完整性清单，可用于检查解压和后续修改 |

CSV均为UTF-8；空值/NA表示缺失，不代表0。`mean`和`sd`是**现有固定数值分析输入中的变化均值及变化SD/输入用SD**，不是本次新增提取的干预前后原始观测。原Excel个别行包含模型换算及SD含义注释，本项目保留其实际分析值，并未声称每个SD均为原文直接报告。不要把它们重新当作后测SD或再计算一次前后相关。

运行入口会检查冻结数据的完整性。若以后修订数据，应另存版本并更新清单；旧数值基准用于识别变化，不应用来迫使新结果贴近旧结果。`python3 scripts/check_integrity.py` 可在解压后、重跑前核验发布文件；运行日志和图形重新生成后，整包文件哈希自然可能改变。

## 可选：原Excel与CSV逐单元格核对、条件性GRADE

R分析本身只需R。以下两项另需 Python 3，安装记录版本的依赖后运行：

```bash
python3 -m pip install -r dependencies/requirements-grade.txt
python3 scripts/check_source_workbook.py
Rscript run_all.R --with-grade
```

若系统默认Python不可用，设置 `PYTHON` 为可用解释器路径。GRADE脚本使用实际贡献矩阵与诊断结果，不是把Excel静态评级复制到输出。其评级是**按用户指定RoB2情景和修改后的降级规则得到的条件性初评**：报告偏倚/间接性固定不降级，其他后三域取最高降级，不冒称完整标准CINeMA或已执行ROB-MEN；详见 `docs/METHODS_AND_LIMITATIONS.md`。

## 数值核查与图形对应

`validation/expected_files.csv` 列出待核对的原结果。数值比较同时使用绝对/相对容差；ID、分类、行列结构、缺失位置须完全一致。R050两位小数的显示匹配容差为0.00501。PDF字节、字体或绘图设备的差异不作为统计复现失败；以生成图表所用的CSV数值为准。

- `results/figures/main/primary_NMA_CNMA_networks_rankings.pdf`：新整理的4页展示，依次为网络、对NAT的NMA森林图、P-score、固定规格CNMA森林图；不是原R050全套版式的重制。
- `results/figures/sensitivity/Figures/`：高风险剔除效应、排名、逐一剔除范围、6张热图及48/43页补充森林图。
- `results/figures/subgroups/subgroup_forests_clean.pdf`：疗程、测量工具、人群亚组及成年限定敏感性。
- `results/figures/meta_regression/meta_regression_core_forests_and_curves.pdf`：疗程与成人基线CAIT斜率、BAL对NAT的拟合曲线。

整理仅调整库路径、输出路径、通用PDF设备、按操作系统选择PNG设备和字体；敏感性PNG改为180 dpi预览，矢量PDF保留。原版PDF也已归档。数值分析和实际图形页数的核查结果见复现报告。

## 解释边界

主网络沿用R050所对应的报告/分析ID。k不能直接写成已确认的独立RCT数；潜在共享队列、测试侧、聚类随机化以及方向映射仍须按原文裁决。RoB2输入包含用户此前调整的评级情景，原始Overall另列保留。数字可复算不等于这些科学判断已通过独立验证；本项目不通过改数据来消除既有问题。

## 上传GitHub

解压后，将**项目文件夹内的全部内容**上传到新仓库根目录，包括 `.github`、`.gitignore`、`.gitattributes`。保留目录层级。只上传ZIP会将项目作为压缩附件保存，不能直接浏览代码或使用工作流。此处已准备上传文件，尚未替你创建仓库、推送或公开。

已附带手动触发的GitHub Actions复现工作流；上传后在 Actions 中运行 `Reproduce analyses`，即可在云端再次核查并下载运行结果。首次云端运行可能需要较长时间安装依赖。

未替研究作者指定开源许可或虚构作者信息；权利说明见 `LICENSE-NOTICE.md`。所有文献数据均为研究层汇总数据，不包含个体参与者资料。

GitHub版本遵循`.gitignore`，不跟踪本机运行日志；本仓库的SHA256清单只覆盖已纳入版本控制的文件。重跑后日志会重新生成。
CSV保留源文件字节与换行符，避免Git自动转换换行导致数据完整性校验失败。
