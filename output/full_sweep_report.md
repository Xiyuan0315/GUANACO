# Brooks-Lint — Full Sweep Report

Mode: Full Sweep | Scope: GUANACO_v3 全项目静态扫描及数据读取、缓存、绘图消费路径重点整理
Date: 2026-09-06
Health Score: 79/100（基于本轮已识别问题的估算，不是覆盖率或完整正确性证明）
Trend: 上次 Full Sweep 96 → 本次 79；上次仅针对 linking 等局部，本次范围不同，不应解读为代码质量下降。

## Summary

保留 Python、MuData、Dash 和 notebook 的公开接口，在现有模块中统一数据访问。
本轮修复 5 组问题；全量测试从 394 passed / 2 failed 变为 412 passed。
103 个生产 Python 文件总行数从 33,524 变为 33,525（净增加 1 行，含注释和空行）。
未修改用户已有的 linking 功能、示例或 abstract 文件；未提交或推送。

## Scope and limitations

- 枚举 269 个 tracked 文件；另外纳入已有未跟踪的 feature_stats.py 及相关测试。
- AST/导入/重复逻辑及风险模式扫描覆盖 103 个生产 Python 文件；语法解析同时覆盖 37 个测试文件（含 conftest）、1 个 scripts Python 文件、9 个 examples Python 文件及 docs/conf.py。
- 对读取层、缓存、绘图消费者、回调入口、配置/分享边界及测试进行重点阅读。不是对所有文档、notebook、图片、生成文件逐行语义审计。
- vendored/generated dash-draggable 运行时仅做静态检查，不重建；没有运行下载或覆盖输出的构建脚本、示例预处理脚本。
- 未启动公网分享、访问真实云端数据、执行桌面打包或浏览器端到端测试；本报告不构成这些路径的验证。
- 无 .brooks-lint.yaml，使用默认风险项。原有工作区修改全部保留。

## Dimension Summary

| Dimension | Scanned | Safe Applied | Extended Applied | Reverted | Residual |
|---|---|---:|---:|---:|---:|
| Review (R1–R6) | 全局静态扫描，读取/绘图热点人工检查 | 0 | 3 | 0 | 1 |
| Test (T1–T6) | 36 个 test 文件，全套 412 项测试 | 0 | 1 | 0 | 1 |
| Debt | 重复读取、缓存策略、遗留接口 | 0 | 1 | 0 | 2 |
| Audit | 103 模块静态导入图及关键入口 | 0 | 0 | 0 | 1 |

Extended 修改按不超过 5 文件的小批次验证；全量基线的两项既有失败先单独记录，各修复使用原本通过的相关测试集，再运行全套回归。

## Findings / Fix Log

### 1. R3 — 延迟稀疏 embedding 的重复物化逻辑（Warning，applied）

Symptom: utils/embeddings.py 先判断 sparse 再 compute，dask 返回 scipy sparse 时会落入 np.asarray(sparse)，不能得到坐标矩阵；计算调度逻辑也与共享读取层分叉。
Source: A Philosophy of Software Design — Information Hiding；Refactoring — Duplicate Code。
Consequence: notebook 和多模态嵌入对同一种延迟稀疏输入行为不一致。
Remedy: embedding_to_numpy 复用 densify_matrix，保留维度验证和 to_memory 兼容处理。
验证：dense/CSR/CSC × eager/lazy 的 6 个数值用例、3 个非法维度用例及现有多模态测试通过。

### 2. R6 — 同一文件的等长筛选视图混用缓存（Critical，applied）

Symptom: gene cache、heatmap bin cache 和 violin cache 用文件名或文件名加形状作为数据标识。
Source: Domain-Driven Design — Identity；Working Effectively with Legacy Code — Characterization Tests。
Consequence: 同一个 h5ad 的前 3 个与后 3 个细胞可能返回相同表达值，影响科学结果正确性。
Remedy: 共用 dataset_cache_token：每个活跃数据对象有进程内唯一标识，视图相互区分；弱引用不延长数据生命周期，单调标识避免 Python id 回收造成碰撞。原导入和函数签名保留。
验证：先复现错误（期望 [6,8,10]，实际 [0,2,4]），修复后 gene/heatmap/violin 三条真实输出路径通过；另测数据对象能被回收。
边界：这不提供原地修改 X/obs 后的自动缓存失效机制；当前仍按准备好且稳定的数据集使用。

### 3. R4 — violin 为单列分组读取全部注释（Warning，applied）

Symptom: plots/violin1.py 的筛选分支对整个 obs 行切片调用 to_memory，再取 groupby。
Source: A Philosophy of Software Design — Pull Complexity Downwards。
Consequence: 分组只需一列，却从磁盘读取不相关注释列；注释越宽浪费越大。
Remedy: 先通过 obs_col 读取分组列，再对 Series 应用行筛选，删除原来的整表物化分支。
验证：真实 lazy Zarr 用例禁止 Dataset2D.to_memory，输出仍为预期 [2,4]。没有把未测量的网络收益写成加速倍数。

### 4. T3 — loader 测试与按需注释契约脱节（Warning，applied）

Symptom: 两个旧测试断言 backed Zarr 的 obs 必须是 pandas DataFrame，文档和状态消息也写“metadata in memory”；实际实现有意保留 lazy obs。
Source: The Art of Unit Testing — Tests of Observable Behavior；xUnit Test Patterns — Fragile Test。
Consequence: 正常的按需读取行为被报告为失败，维护者可能为迎合测试恢复整表加载。
Remedy: 保留 X/obs 延迟与 var 已物化断言，检查读取层不主动 eager-read obs，并验证注释值、细胞索引和可选分组。同步 loader 注释、状态消息及 DATA_LOADING.md。
验证：两项旧失败已消除；保留原有读取行为，没有为了过测试改回 eager obs。

### 5. R3/R5 — 多个入口重复定义注释读取适配器（Warning，applied）

Symptom: loader、widget、marimo、capabilities 各自实现 obs[col] → to_series；多个绘图模块仅为读取一列注释依赖完整 loader。
Source: The Pragmatic Programmer — DRY；A Philosophy of Software Design — Information Hiding。
Consequence: 存储兼容性需要多处同步，基础注释工具向加载层反向依赖。
Remedy: obs_col 移入已有 utils/obs_utils.py，网页、能力识别、绘图、回调及 notebook 均复用；loader.obs_col 和 notebook 内部别名保留。
验证：全局搜索仅剩一个 obs_col 实现；项目内部不再从 loader 导入 obs_col；pandas 分类列及 lazy Zarr 的行为测试通过。
Debt score: Pain 2 × Spread 3 = 6。

## Verification

- 基线：394 passed, 2 failed（两项为上述旧注释契约断言）。
- 最终：412 passed, 38 warnings，75.21 秒。
- 命令：.pixi/envs/default/bin/python -m pytest -o addopts='' -q
- .pixi/envs/default/bin/ruff check src tests scripts：通过。
- git diff --check：通过。
- 三个 shell 脚本的 bash -n 语法检查通过；没有执行构建或生物数据转换。
- 保留第三方弃用警告和非有限值测试产生的数值警告，没有全局屏蔽来制造“干净”输出。

数值与性能复测：200,000 cells × 80 features，20 groups，稀疏矩阵密度 1%，3 次交错运行的中位数。对照为前一轮公共读取器整理前保存的实现，已经包含此前稀疏统计优化，不是最初未优化版本。

| Storage | 对照耗时 ms | 本轮耗时 ms | 对照/本轮临时峰值 MiB |
|---|---:|---:|---:|
| CSC | 59.80 | 59.58 | 9.116 / 9.116 |
| CSR | 84.18 | 85.18 | 9.115 / 9.116 |
| Dense | 290.71 | 288.46 | 19.340 / 19.342 |

dot plot 颜色与点大小数组完全一致；耗时变化在约 1.2% 内，不据此声称新加速。tracemalloc 峰值不含已加载源矩阵，不是进程 RSS。此前优化的较大收益见 optimization_results.md，不能与本轮维护性收益重复计算。

## Architecture notes

注释消费者现在直接依赖 obs_utils；embedding 物化使用共享读取器；gene、heatmap、violin 复用数据身份规则。保留现有包布局，没有新建缓存框架或存储抽象层。

静态图纳入绝对、相对及函数内导入。唯一发现的 guanaco ↔ guanaco.marimo 环是函数内延迟访问顶层 pl 的已知调用方式，不是已复现的初始化循环，因此不为消除图上的环而改写 API。utils/memory_utils.py 的遗留 load_adata 延迟依赖见 residual。

## Residual Items (5 not applied)

### Warning — 缓存总内存预算尚未统一（R4，Debt 2×3=6）

Symptom: callbacks/register.py 的 FigureMemoCache 按 24 项限制，plots/violin1.py 的数据缓存按 50 项限制；每项可包含随细胞数增长的数组。
Source: Release It! — Capacity；A Philosophy of Software Design — Information Hiding。
Consequence: gene/ATAC 的字节预算不能保证整个多用户进程的总内存上限。
Remedy: 先确定进程级预算、对象共享的计费方式及多会话压力基线，再统一可复用的限额策略。
Not applied because: remedy ambiguous；直接引入通用缓存层或估算字节数会改变缓存命中与延迟，需要独立容量测试。本轮没有声称全项目内存已严格封顶。

### Warning — import-time 运行时注册表（R5，Architecture）

Symptom: data/registry.py 在模块导入时 load_config / initialize_data；callback register 通过颜色配置依赖此模块。
Source: Clean Architecture — Dependency Rule；Working Effectively with Legacy Code — Seams。
Consequence: 单独复用部分回调仍隐含依赖全局配置和初始化顺序。
Remedy: 后续把运行时 registry 显式传入应用组装入口，先补启动顺序与多配置实例测试。
Not applied because: cross-module structural change；涉及应用组装契约，不在这轮安全自动整理范围。

### Warning — 桌面、浏览器与公网路径缺乏本轮端到端验证（T2）

Symptom: 分享测试使用 Flask test_client 并替换 tunnel 函数；桌面打包和真实浏览器交互未在本轮执行。
Source: How Google Tests Software — Test Pyramid；The Art of Unit Testing — Integration Tests。
Consequence: Python 回归通过不能证明打包资源齐全、真实隧道可达或浏览器拖拽/联动正确。
Remedy: 用合成公开数据增加桌面启动、浏览器操作和可选受控隧道的 smoke tests。
Not applied because: no test infrastructure；本轮没有为验证而公开用户数据或启动外部服务。

### Warning — motif shell 流程可能在失败后继续（R4）

Symptom: scripts/motif_extraction/motif_atac_intersect.sh 串行调用外部转换工具，无 fail-fast；中间文件使用固定工作目录文件名。
Source: Release It! — Fail Fast；The Pragmatic Programmer — Design by Contract。
Consequence: 工具失败后后续步骤仍可能运行，多次并行任务也可能争用中间文件。
Remedy: 用独立临时目录、显式输出策略与失败退出，并以小型基因组 fixture 验证外部命令链。
Not applied because: no test coverage / output contract；未执行缺乏输入的生物预处理，更未修改已有输出路径。

### Suggestion — 遗留 memory_utils 接口（R4，Debt 1×1=1）

后续状态（2026-09-06）：已按用户要求完成本项收尾。移除未使用的私有元数据读取及状态；旧基因入口改为共享读取器的兼容包装；保留公开名称以免外部调用失效。以下为原始发现，保留作历史记录。

Symptom: memory_utils.py 的 LazyAnnData、memory_efficient_gene_expression 等在仓库内没有找到消费者，且仍维护另一套读取路径。
Source: Refactoring — Dead Code；The Pragmatic Programmer — DRY。
Consequence: 新代码可能误用未接入共享 backed/raw/layer 逻辑的旧接口。
Remedy: 确认外部调用情况后弃用/删除，或只保留薄兼容包装；不要继续往这个模块增加独立读取实现。
Not applied because: public API removal / compatibility ambiguity；仓库内无引用并不能证明外部 notebook 未使用，故未直接删除。

## Iteration History

- Round 1: mixed；完成 5 组安全修复，按小批次跑相关测试。
- Round 2: 复查变更、同模块及直接消费者；补齐全项目注释导入迁移，未发现新的可直接安全修复项；最终全量测试通过。
- Stopped at: no outstanding criticals；保留上述 5 项需要策略或验证条件的问题。
- 无 3-retry 退休项，无破坏性 API 修改。测试编写中曾误假设 lazy Series 保留 categorical dtype，已改为验证真实的值、索引和分组契约，未以此改变生产语义。

## Health Score Delta

Before: 44/100 → After: 79/100（估算）
计分为已识别问题：修复 1 Critical + 4 Warning；剩余 4 Warning + 1 Suggestion。
这些分数依赖本轮检查范围，不代表未检查路径无问题。

## Scope Inventory

### Legacy cleanup follow-up — 2026-09-06

Mode: PR Review（按用户要求实施修复）
Scope: 本轮仅检查 gene_extraction_utils.py、violin1.py、memory_utils.py 的旧入口及 test_memory_utils.py，不重审其他工作区修改。
Health Score: 100/100（仅此局部变更未发现待修复项，不代表全项目评分）。
Trend: 上次 PR Review 82 → 本次 100；范围不同，不用于判断全项目趋势。

#### Findings / Fix Summary

- R4 — Middle Man，已修复。
  Symptom: _adata_id / _get_adata_id 仅转发共享标识，或附带调用处已有的 shape。
  Source: Refactoring — Middle Man。
  Consequence: 多一层跳转，没有隐藏独立规则。
  Remedy: 删除两个内部函数，调用处直接使用 dataset_cache_token，缓存键内容保持不变。
- R3 — 重复表达读取，已修复。
  Symptom: memory_efficient_gene_expression 自行读取 X 并判断稀疏类型。
  Source: The Pragmatic Programmer — DRY。
  Consequence: 旧入口未受益于统一的磁盘/延迟数组读取兼容逻辑。
  Remedy: 调用 extract_gene_expression(use_cache=False, dtype=None)，保留旧的 log、无裁剪总体标准差 z-score 及常量行为。
- R4 — Dead Code，已修复。
  Symptom: LazyAnnData._load_metadata 和 _metadata 没有仓库内调用；sparse_safe_slice 重复按类型分支。
  Source: Refactoring — Dead Code / Duplicate Code。
  Consequence: 保留另一套未用的 HDF5 读取假设与不必要分支。
  Remedy: 删除私有方法和字段、简化切片；LazyAnnData 仅保留首次访问时调用 load_adata 的功能。

Summary: 3 个生产文件净减少 28 行，删除 3 个内部函数。公开兼容入口保留；embedding_to_numpy 的维度检查、分块转换及监控工具有独立职责，不作为纯包装删除。没有新增读取框架。
验证：新增 41 个兼容性/读取测试，先运行修改前行为测试，再检查 dense/CSR/CSC、磁盘/延迟筛选、dtype、变换、缺失基因和延迟加载。最终全量 453 passed，38 个既有警告；Ruff 和 git diff --check 通过。没有修改用户数据或删除整个文件；删除的旧代码可从版本控制恢复。

### Original inventory

以下是检查开始时的 tracked 文件清单；生成/二进制文件仅纳入枚举，不声称逐行审查。另纳入未跟踪的 src/guanaco/utils/feature_stats.py、tests/test_feature_stats.py、tests/test_callback_render_state.py。

```text
.brooks-lint-history.json
.gitattributes
.gitignore
Dockerfile
INSTALL.md
LICENSE
README.md
docs/.buildinfo
docs/.nojekyll
docs/Makefile
docs/assets/ConfigJson.png
docs/assets/Figure21.png
docs/assets/Figure22.png
docs/assets/Figure23.png
docs/assets/Figure3.png
docs/assets/GUANACOs.png
docs/assets/bar.png
docs/assets/configguanaco.png
docs/assets/dotplot.png
docs/assets/download.png
docs/assets/footprint.png
docs/assets/gb_motif.png
docs/assets/guanaco-interface.png
docs/assets/heatmap.png
docs/assets/interface1.png
docs/assets/pseudotime.png
docs/assets/selection.png
docs/assets/violin1.png
docs/assets/violin2.png
docs/assets/violin2_mode3.png
docs/assets/z-score_heatmap.png
docs/conf.py
docs/deployment/configuring-json.html
docs/deployment/configuring-json.rst
docs/deployment/index.html
docs/deployment/index.rst
docs/deployment/installation.html
docs/deployment/installation.rst
docs/deployment/private-s3.rst
docs/deployment/running-python.html
docs/deployment/running-python.rst
docs/development/colormap_previews/continuous_colormap_overview.png
docs/development/colormap_previews/discrete_colormap_overview.png
docs/development/colors.rst
docs/development/continuous_colormaps.md
docs/development/discrete_colormaps.md
docs/development/expression_trend_plot.md
docs/development/index.rst
docs/development/paga_plot.md
docs/development/volcano_plot.md
docs/getting-started.rst
docs/index.html
docs/index.rst
docs/linked_views.md
docs/make.bat
docs/marimo_guide.md
docs/notebook_api.md
docs/objects.inv
docs/searchindex.js
examples/CLL_CLINICAL_MULTIOMICS.md
examples/configs/cll_clinical_multiomics.json
examples/configs/cloud_zarr_config.json
examples/configs/cluster_1_cloud_zarr_peakbrowser.json
examples/configs/cluster_2_three_datasets.json
examples/configs/pbmc_atac_browser_test.json
examples/configs/pbmc_liana.json
examples/configs/pbmc_lr_visualization.json
examples/configs/pbmc_unpaired_rna_atac.json
examples/configs/visium_hne_spatial.json
examples/data/pbmc_ligand_receptor_demo.csv
examples/liana/requirements.txt
examples/linked_views/README.md
examples/linked_views/demo_data.py
examples/marimo/__marimo__/session/guanaco_explorer.py.json
examples/marimo/__marimo__/session/guanaco_panels.py.json
examples/marimo/guanaco_explorer.py
examples/marimo/guanaco_panels.html
examples/marimo/guanaco_panels.py
examples/marimo/guanaco_template.py
examples/notebooks/Linked_views_demo.ipynb
examples/notebooks/TACoWig_Template.ipynb
examples/scripts/build_unpaired_pbmc_example.py
examples/scripts/prepare_cll_clinical_multiomics.py
examples/scripts/prepare_cll_rds_export.R
examples/scripts/prepare_pbmc_lr_visualization_demo.py
examples/scripts/prepare_visium_hne_spatial.py
examples/scripts/run_pbmc_liana.py
examples/widgets/README.md
examples/widgets/Widget_example.ipynb
examples/widgets/Widget_example_offline.html
examples/widgets/umap_example.html
guanaco.json
pixi.lock
pixi.toml
pyproject.toml
scripts/build_dash_draggable.sh
scripts/build_desktop.sh
scripts/desktop_entry.py
scripts/motif_extraction/README_motif_atac.md
scripts/motif_extraction/motif_atac_intersect.sh
src/dash_draggable/LICENSE
src/dash_draggable/ResponsiveGridLayout.py
src/dash_draggable/__init__.py
src/dash_draggable/_imports_.py
src/dash_draggable/dash_draggable.min.js
src/dash_draggable/dash_draggable.min.js.map
src/dash_draggable/metadata.json
src/dash_draggable/package-info.json
src/guanaco/__init__.py
src/guanaco/app.py
src/guanaco/assets/configguanaco.png
src/guanaco/assets/favicon.ico
src/guanaco/assets/footprint.png
src/guanaco/assets/lamp_guanaco.png
src/guanaco/assets/logo.png
src/guanaco/assets/scientific_style.css
src/guanaco/cli.py
src/guanaco/config_wizard.py
src/guanaco/dash_app.py
src/guanaco/data/DATA_LOADING.md
src/guanaco/data/__init__.py
src/guanaco/data/capabilities.py
src/guanaco/data/capability_schema.py
src/guanaco/data/hdf5_capabilities.py
src/guanaco/data/ligand_receptor.py
src/guanaco/data/loader.py
src/guanaco/data/multiomics.py
src/guanaco/data/registry.py
src/guanaco/desktop.py
src/guanaco/layouts/__init__.py
src/guanaco/layouts/app_layout.py
src/guanaco/linking/__init__.py
src/guanaco/linking/base.py
src/guanaco/linking/data.py
src/guanaco/linking/engine.py
src/guanaco/linking/model.py
src/guanaco/linking/native_adapters.py
src/guanaco/linking/registry.py
src/guanaco/linking/runtime.py
src/guanaco/linking/table_adapters.py
src/guanaco/main.py
src/guanaco/marimo.py
src/guanaco/pages/__init__.py
src/guanaco/pages/matrix/__init__.py
src/guanaco/pages/matrix/analysis/__init__.py
src/guanaco/pages/matrix/analysis/composition_da.py
src/guanaco/pages/matrix/callbacks/__init__.py
src/guanaco/pages/matrix/callbacks/atac_browser_callbacks.py
src/guanaco/pages/matrix/callbacks/cross_modal_concordance_callbacks.py
src/guanaco/pages/matrix/callbacks/dotplot_callbacks.py
src/guanaco/pages/matrix/callbacks/guanaco.json
src/guanaco/pages/matrix/callbacks/heatmap_callbacks.py
src/guanaco/pages/matrix/callbacks/ligand_receptor_callbacks.py
src/guanaco/pages/matrix/callbacks/multiomics_composition_callbacks.py
src/guanaco/pages/matrix/callbacks/paga_callbacks.py
src/guanaco/pages/matrix/callbacks/pseudotime_callbacks.py
src/guanaco/pages/matrix/callbacks/register.py
src/guanaco/pages/matrix/callbacks/scatter_callbacks.py
src/guanaco/pages/matrix/callbacks/spatial_relationships_callbacks.py
src/guanaco/pages/matrix/callbacks/stacked_bar_callbacks.py
src/guanaco/pages/matrix/callbacks/unpaired_multiomics_callbacks.py
src/guanaco/pages/matrix/callbacks/violin_callbacks.py
src/guanaco/pages/matrix/callbacks/volcano_callbacks.py
src/guanaco/pages/matrix/layouts/__init__.py
src/guanaco/pages/matrix/layouts/atac_browser_layout.py
src/guanaco/pages/matrix/layouts/cross_modal_concordance_layout.py
src/guanaco/pages/matrix/layouts/dotplot_layout.py
src/guanaco/pages/matrix/layouts/embedding_layout.py
src/guanaco/pages/matrix/layouts/heatmap_layout.py
src/guanaco/pages/matrix/layouts/ligand_receptor_layout.py
src/guanaco/pages/matrix/layouts/multiomics_composition_layout.py
src/guanaco/pages/matrix/layouts/paga_layout.py
src/guanaco/pages/matrix/layouts/pseudotime_layout.py
src/guanaco/pages/matrix/layouts/spatial_relationships_layout.py
src/guanaco/pages/matrix/layouts/stacked_bar_layout.py
src/guanaco/pages/matrix/layouts/violin_layout.py
src/guanaco/pages/matrix/layouts/volcano_layout.py
src/guanaco/pages/matrix/plots/__init__.py
src/guanaco/pages/matrix/plots/atac_browser.py
src/guanaco/pages/matrix/plots/cross_modal_concordance.py
src/guanaco/pages/matrix/plots/dotmatrix.py
src/guanaco/pages/matrix/plots/embedding.py
src/guanaco/pages/matrix/plots/gene_annotation.py
src/guanaco/pages/matrix/plots/guanaco.json
src/guanaco/pages/matrix/plots/heatmap.py
src/guanaco/pages/matrix/plots/ligand_receptor.py
src/guanaco/pages/matrix/plots/multiomics_composition.py
src/guanaco/pages/matrix/plots/paga.py
src/guanaco/pages/matrix/plots/pseudotime.py
src/guanaco/pages/matrix/plots/spatial_relationships.py
src/guanaco/pages/matrix/plots/stacked_bar.py
src/guanaco/pages/matrix/plots/violin1.py
src/guanaco/pages/matrix/plots/violin2.py
src/guanaco/pages/matrix/plots/volcano.py
src/guanaco/pages/visualizations/__init__.py
src/guanaco/pages/visualizations/callbacks.py
src/guanaco/pages/visualizations/layout.py
src/guanaco/pages/visualizations/plots/__init__.py
src/guanaco/pages/visualizations/plots/igv/__init__.py
src/guanaco/pages/visualizations/plots/igv/callbacks.py
src/guanaco/pages/visualizations/plots/igv/layout.py
src/guanaco/pages/visualizations/plots/igv/motif.py
src/guanaco/pages/visualizations/registry.py
src/guanaco/share.py
src/guanaco/utils/__init__.py
src/guanaco/utils/colors.py
src/guanaco/utils/cvd_color.json
src/guanaco/utils/cytoscape.py
src/guanaco/utils/embeddings.py
src/guanaco/utils/gene_extraction_utils.py
src/guanaco/utils/memory_utils.py
src/guanaco/utils/obs_utils.py
src/guanaco/utils/plot_config.py
src/guanaco/utils/plot_style.py
src/guanaco/utils/progress_utils.py
src/guanaco/utils/render_guard.py
src/guanaco/utils/search.py
src/guanaco/utils/ui_helpers.py
src/guanaco/utils/volcano_utils.py
src/guanaco/widget.py
tests/conftest.py
tests/test_atac_browser.py
tests/test_atac_browser_callbacks.py
tests/test_callback_boundaries.py
tests/test_composition.py
tests/test_composition_da.py
tests/test_config_wizard.py
tests/test_cross_modal_concordance.py
tests/test_dotmatrix.py
tests/test_embedding.py
tests/test_heatmap.py
tests/test_ligand_receptor.py
tests/test_linked_view_demo_data.py
tests/test_linking_data.py
tests/test_linking_engine.py
tests/test_linking_integration.py
tests/test_linking_model.py
tests/test_linking_native_adapters.py
tests/test_linking_public_api.py
tests/test_linking_runtime.py
tests/test_loader.py
tests/test_multiomics.py
tests/test_multiomics_composition.py
tests/test_peak_browser.py
tests/test_pseudotime.py
tests/test_scatter_callbacks.py
tests/test_search.py
tests/test_share.py
tests/test_spatial_relationships.py
tests/test_table_adapters.py
tests/test_ui_helpers.py
tests/test_violin1.py
tests/test_violin2.py
tests/test_visualization_sections.py
tests/test_widget_helpers.py
uv.lock
vendor/dash_draggable/.babelrc
vendor/dash_draggable/.npmrc
vendor/dash_draggable/LICENSE
vendor/dash_draggable/README.md
vendor/dash_draggable/README.upstream.md
vendor/dash_draggable/package-lock.json
vendor/dash_draggable/package.json
vendor/dash_draggable/src/lib/components/ResponsiveGridLayout.react.js
vendor/dash_draggable/src/lib/components/style.css
vendor/dash_draggable/src/lib/constants.js
vendor/dash_draggable/src/lib/index.js
vendor/dash_draggable/src/lib/localStorage.js
vendor/dash_draggable/webpack.config.js
```
