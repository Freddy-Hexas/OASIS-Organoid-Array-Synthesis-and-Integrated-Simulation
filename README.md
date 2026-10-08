# OASIS: Organoid Array Synthesis and Integrated Simulation

[English](README.en.md) · 中文

项目仓库：[OASIS](https://github.com/Freddy-Hexas/OASIS-Organoid-Array-Synthesis-and-Integrated-Simulation)

OASIS 读取 GDS 中的支撑结构，自动放置圆形电极，生成导线，并安排外部四边 Pad。你可以在网页中查看输入、几何分析和输出。每个 GDS 对应一个任务。不同 GDS 可以同时运行。

v4 是从 v2 整理出的独立运行目录。它保留现有求解方法，统一输入和输出路径。它不包含历史任务、研究文档、日志或缓存。首次启动时，所有案例均显示为待运行。

## 先理解这个工具

把原始支撑结构想成一组弯曲的道路。电极是采集信号的小圆盘。导线沿这些道路走向外部。Pad 是较大的金属接点，用于连接外部设备。

一个有效输出必须包含完整连接：**一个电极 → 一条独立的金属网络 → 一个专属 Pad**。导线不能在没有支撑的位置直接跨过去。工具可以在电极位置增加小圆形支撑，也可以从允许的外端出口增加承载桥。

### 术语表

| 术语 | 本项目中的含义 | 通俗解释 |
| --- | --- | --- |
| GDS | 存储几何轮廓、单元和图层的版图文件 | 描述形状的文件；像素图片不能代替它 |
| 支撑结构 / substrate | 选定 GDS 图层中的支撑材料轮廓 | 导线下面的道路 |
| 电极 / electrode | 用于采集信号的圆形金属区域 | 放在目标区域内的小圆盘 |
| 圆岛 / island | 在电极附近添加的局部圆形支撑 | 贴在原结构上的圆形底座 |
| 导线 / wire | 连接电极与 Pad 的金属区域 | 沿支撑结构铺设的连接线 |
| 网络 / net | 属于同一个电极连接的全部金属 | 必须连通，并与其他网络分开 |
| Pad | 外部矩形金属接点 | 连接外部设备的较大接点 |
| 承载桥 / bridge | 从允许出口到外部 Pad 区的新增支撑 | 结构外部的连接道路 |
| 图层 / layer | GDS 的数字层号与 datatype | 例如 `10/0`；数字本身不表示材料 |
| 走廊 / corridor | 可用于导线通过的支撑通道 | 一段可通行的道路 |
| 净距 / clearance | 两个边界之间的最短距离 | 用边缘计算，不用中心计算 |
| 下界 L | 经规定检查确认的完整连接数 | 已经做出来并检查过的数量 |
| 上界 U | 声明模型内任何合法布局都不能超过的数量界 | 理论上限，可能比较宽松 |

### 文档的 ASD-STE100 写作方式

ASD-STE100 是一种受控技术英语写作标准。它包含写作规则和受控词典。本说明采用其清晰表达原则：短句、主动语态、固定术语和直接操作步骤。中文说明采用对应的表达方式。英文版也解释专业术语。参见 [ASD-STE100 介绍](https://www.asd-ste100.org/about_STE.html)和[官方问答](https://www.asd-ste100.org/STE_faq.html)。

例如：“如果你要添加输入文件，把 GDS 放到 `data/`，然后点击刷新。”这句话先给条件，再给动作。每个操作步骤说明动作和预期结果。

本文没有进行完整的受控词典审查，不声明已获得 ASD-STE100 符合性认证。

## 目录

```text
.
├── README.md                 # 中文说明
├── README.en.md              # 英文说明
├── start_workbench.py        # 推荐启动入口
├── start_workbench.ps1       # PowerShell 启动包装
├── requirements.txt          # 运行依赖
├── requirements-dev.txt      # 可选浏览器验证依赖
├── .gitignore                # 排除生成结果和本机环境
├── data/                     # 原始输入 GDS
├── gds_frontend/             # 几何、布局、布线和审计源码
│   ├── process_rules.json    # 默认工艺规则
│   ├── verify_*.py           # 必要的验证和证书复核工具
│   └── web_app/              # 主页面、批量页面和 HTTP 服务
│       └── vendor/katex/     # 本地公式资源、字体与第三方许可证
├── tests/
│   └── verify_package.py     # 独立目录、网页和完整任务验证
└── outputs/                  # 运行时自动创建，不随源码发布
    └── runs/                 # 任务结果与运行缓存
```

`data/` 包含以下 9 份结构。GDS 内容与整理时的 v2 输入相同。

```text
01_curved_hex_reference_like.gds
A_spiral_basket.gds
B_serpentine_honeycomb.gds
C_open_petal_mesh.gds
XingMing1.gds
organoid_hexagon_lattice.gds
organoid_rhombus_lattice.gds
organoid_square_lattice.gds
organoid_triangle_lattice.gds
```

## 安装与启动

### 环境

使用 Python 3.11 或更新版本。Windows / Python 3.13.9 已用于本次验证。Linux 和 macOS 启动命令列在下方，但本次没有验证这两个系统。

前端直接使用 HTML、CSS 和 JavaScript，无须 Node.js 或 npm 构建。KaTeX 在目录中提供。安装 Python 依赖后，页面和公式渲染可以在本机离线使用。

### Windows PowerShell

1. 打开 PowerShell。
2. 进入项目根目录。移动项目后，用新的路径替换下方路径。

```powershell
Set-Location -LiteralPath 'A:\Electrode_design\结构\布线布电极v4'
```

3. 创建虚拟环境并安装依赖。

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

如果本机没有 `py` 命令，使用已安装的 Python 执行 `python -m venv .venv`。以下命令直接调用虚拟环境解释器，不需要修改 PowerShell 的激活脚本策略。

4. 启动服务。

```powershell
.\.venv\Scripts\python.exe -B start_workbench.py
```

5. 打开主页面：<http://127.0.0.1:8769/>。
6. 打开批量页面：<http://127.0.0.1:8769/admin>。

启动窗口会显示输入目录、输出目录和并发数量。保持窗口打开。按 `Ctrl+C` 停止 HTTP 服务；已经提交的任务会等待执行结束。关闭浏览器不会取消后台任务。

已有运行依赖的 Python 也可以直接启动：

```powershell
python -B start_workbench.py --port 8769 --workers 2
```

`start_workbench.ps1` 优先使用项目内的 `.venv`，否则使用 PATH 中的 `python`。

### Linux / macOS

在项目根目录执行：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B start_workbench.py
```

### 启动参数

| 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `--host` | `127.0.0.1` | 监听地址；默认仅本机访问 |
| `--port` | `8769` | 网页端口 |
| `--workers` | `2` | 同时计算的任务数，可选 1–4 |
| `--input-dir` | 项目内 `data/` | 输入 GDS 目录 |
| `--output-dir` | 项目内 `outputs/runs/` | 结果、任务状态和缓存目录 |

默认路径根据脚本位置确定，不依赖启动时的工作目录。显式传入的相对路径按当前工作目录解析。

```powershell
python -B start_workbench.py --input-dir 'D:\GDS inputs' --output-dir 'D:\GDS results' --workers 2
```

如果需要局域网访问，执行：

```powershell
python -B start_workbench.py --host 0.0.0.0 --port 8769
```

其他设备使用服务器的实际 IP 地址访问。`0.0.0.0` 是监听地址，不能作为远程访问地址。系统防火墙和路由配置决定设备之间是否可达。当前 `/admin` 是批量控制台，没有用户登录或访问权限管理，适用于本机或受控网络。

直接运行 `gds_frontend/web_app/server.py` 时，还可用 `GDS_WORKBENCH_HOST`、`GDS_WORKBENCH_PORT`、`GDS_WORKBENCH_WORKERS`、`GDS_WORKBENCH_DATA_DIR` 和 `GDS_WORKBENCH_RUNS_DIR` 环境变量配置。推荐入口使用上述命令行参数覆盖这些变量。

## 操作一个任务

1. 在左侧选择一个 GDS。页面显示原始矢量轮廓、图层、文件校验值和结构读数。
2. 确认“用于分析的支撑层”。工具优先建议 `10/0`，否则建议多边形最多的层。你必须根据文件来源确认材料含义。
3. 选择“连接四边 Pad”。此模式执行电极到实体 Pad 的完整流程。
4. 填写工艺参数。每个输入框标有单位。
5. 点击“开始分析与试布线”。页面显示计算阶段与队列状态。
6. 等待完成。查看上方的“已布置电极”，以及下方的下界 L、上界 U 和审计状态。
7. 点击阶段卡片查看原始支撑、锚点可放区、走廊与出口、电极与金属、解析曲线。
8. 下载 `routing.gds` 和 `summary.json`。需要复核线宽时，同时保留 `routing.width_witness.json`。

主页面的大数字是经 GDS 回读的电极数。严格审计未通过时，这个数字不能作为认证下界 L。读取结果时，同时检查 `capacity_interval.integer_polygon_lower_verified`。

“仅几何终端试布线”用于分析电极到几何出口的路线。该模式没有认证实体 Pad 连接，其数量不能当作完整 Pad 模式的电极数。

### 查看细节

用滚轮放大或缩小。按住鼠标拖动平移。点击“适合窗口”重置视图。原始结构与最终布局根据 GDS 多边形生成矢量视图，缩放不受 PNG 像素限制。走廊提取图使用 PNG。

解析曲线卡片仅在结果含曲线段时可用。直线路线可能没有该视图。GDS 输出仍采用多边形；曲线视图展示解析中心线，不能据此认定加工曲率已获独立认证。

## 参数与单位

`1 mm = 1000 μm`。例如 `0.07 mm = 70 μm`。API 和 JSON 中以 `_um` 结尾的参数使用 μm。

| 参数 | 默认值 | 对结果的影响 |
| --- | --- | --- |
| 电极直径 | 30 μm | 电极金属圆盘的直径 |
| 最小电极中心距 | 0.07 mm | 任意两个电极中心之间的距离约束 |
| 电极中心允许半径 | 3 mm | 电极中心必须位于该半径圆内；默认直径为 6 mm |
| 导线宽度 | 5 μm | 本次生成和检查的导线宽度 |
| 网络间距 | 4 μm | 不同金属网络之间的边缘净距；也用于圆岛间距 |
| 支撑余量 | 4 μm | 金属边缘与支撑边缘之间的余量 |
| 起始外框边长 | 32 mm | 外部方形 Pad 区的起始尺寸 |
| Pad 沿边宽度 | 0.5 mm | Pad 沿方框一条边的尺寸 |
| Pad 径向长度 | 3 mm | Pad 沿方框内外方向的尺寸 |
| Pad 中心节距 | 1 mm | 相邻 Pad 中心沿边的距离 |
| Pad 调整模式 | 保留尺寸，自动扩框 | Pad 需求增加时扩大外框 |
| 最小 Pad 宽 / 长 | 0.5 / 3 mm | 选择缩小模式时允许的下限 |

电极位置圆以原支撑的最小外接圆中心为参考，并按 GDS 原生网格取整。此规则适用于平移后的结构。导线和 Pad 可以位于电极位置圆外。

如果类器官直径是 2–3 mm，可将电极中心允许半径改为 1–1.5 mm。这里填写的是半径，不是直径。当前约束检查电极中心，不强制整个电极圆盘都位于该圆内。

尺寸必须大于 0。间距、支撑余量和额外中心距可以设为 0。零间距仍须保持网络独立。圆岛与金属的实际大小可能使有效中心距比输入的中心距更大。

Pad 节距必须不小于 Pad 宽度加网络间距。Pad 宽度和长度必须容纳导线与两侧余量。每边使用的 Pad 连续排列。每边 Pad 数量由需求和几何决定。

“先缩小到下限，再扩框”采用统一尺寸缩放策略。如果下限与参考尺寸相同，Pad 不缩小。最终尺寸、节距和外框以运行报告为准。

完整 Pad 模式使用原结构的外端窗口。“宽区识别阈值”仅用于几何终端模式，不控制完整 Pad 模式的接出位置。

## 批量与并发

1. 打开 `/admin`。
2. 选择多个 GDS，或点击“选择所有可运行”。
3. 按需确认每个 GDS 的支撑层。
4. 设置本批统一工艺规则与 Pad 尺寸。
5. 点击“批量提交所选任务”。
6. 查看队列，点击单个结果返回主页面。

默认两个任务同时计算。其余任务排队。`--workers` 调整并发数量，不是网页线程数量。最多 32 个任务可处于运行或排队状态。同一个 GDS 运行或排队时不能重复提交。批次在提交前统一校验；一个参数或图层错误会阻止整批入队。

当前默认没有求解墙钟时间预算。复杂结构可能运行较长时间。候选数量、搜索轮次等离散限制仍然存在。不限时不等于穷尽全部连续位置与路线。并发数越大，通常需要越多内存。

## 添加输入与保存结果

把原始结构 GDS 放到 `data/` 顶层，然后点击刷新。支持 `.gds` 和 `.GDS` 扩展名。不扫描子目录。输入只需包含支撑几何，不需要生成器脚本、预设中心线或连接关系。

输入 GDS 的单元、引用和原生单位由读取器处理。你仍需选对支撑层。零宽点接触不等于有限线宽可通行的通道。默认 5 μm 线宽和两侧各 4 μm 支撑余量，在直通截面至少需要约 13 μm 支撑宽度；转弯和多轨可能需要更多空间。

不要把生成的 `routing.gds` 放回输入目录作为原始支撑。输出包含新增圆岛、桥和外框，重新输入会改变问题定义。输入变更后，已有任务可标为过期，需要重新运行。

每次运行使用独立目录 `outputs/runs/<job_id>/`：

| 文件 | 内容 |
| --- | --- |
| `status.json` | 任务状态、规则、进度和结果 |
| `summary.json` | 完整几何、布局、布线和容量报告 |
| `routing.gds` | 支撑、电极、金属、Pad 等输出层；成功构造后生成 |
| `routing.width_witness.json` | 金属线宽复核所需的中心线证据 |
| `integer_polygon_audit.json` | 导出 GDS 的整数网格几何审计 |
| `navigation_diagnostics.json` | 导航余量、拓扑检查和分量排除记录 |
| `graph.json.gz` / `regions.json.gz` | 压缩的结构图与矢量可放区 |
| `graph.png` / `routing.png` | 预览图片 |
| `error.log` | 失败时的错误详情 |

没有构造出完整路线时，部分文件不会生成。报告中的 `output_layers` 记录实际层号；不要假设所有输入都使用相同输出层号。

显示与导航缓存也存于 `outputs/runs/`。`.gitignore` 排除整个 `outputs/`。服务重启后可以查看已完成任务；运行中任务不会自动恢复。清理结果前，先等待任务结束并停止服务。然后删除需要清理的任务目录。删除整个 `outputs/` 会同时删除结果和缓存，工具会在下次运行时重新创建。

## 方法与结果边界

整个流程保持统一：

```text
GDS + 工艺规则
  → 读取支撑层并检查几何拓扑
  → 构造导航域、走廊、节点窗口和外端出口
  → 生成圆岛附着位置与多轨路线候选
  → 联合选择电极、路线、出口和 Pad
  → 在保持数量的前提下改进中心位置
  → 导出金属与支撑 GDS
  → 回读检查 + 整数网格审计
  → 报告电极数量及容量区间 [L, U]
```

几何前端从实际 GDS 多边形提取信息。数值导航余量按几何特征调整，检查连通块和孔洞，记录在声明模型下可证明排除的微小分量。无法可靠提取的几何会报告错误。

求解器先增加可接通电极数，再在数量相同时优先靠近中心。它将候选完整路线作为整数选择变量，用冲突约束表示中心距、金属净距、圆岛、出口和 Pad 占用。宽支撑允许候选多轨；几何检查决定哪些轨道能同时使用。圆岛附着搜索、残余路线增广和中心重布线进一步改进构造结果。

外部桥遵守从原支撑外端接出的政策。桥和四边 Pad 的布局也参与几何检查。曲线中心线使用切向连续的圆角提案，输出前按误差要求采样为 GDS 多边形。

对声明的模型，真实最大数量记为 N*。可信界满足：

$$
L \le N^* \le U.
$$

L 来自已导出并通过规定审计的完整电极—Pad 网络。U 来自独立于候选库的几何装填界和适用的割界。若 L = U 且全部证明前提通过，可以声明**当前几何模型内数量最多**。若 L < U，只能报告已构造的数量和未闭合的区间。有限候选库的最优解本身不能证明连续几何最多。

检查这些结果字段：

```text
result.capacity_interval.lower_bound
result.capacity_interval.upper_bound
result.capacity_interval.integer_polygon_lower_verified
result.capacity_interval.declared_geometric_model_optimality_proven
result.routing.gds_roundtrip_audit.passed
result.routing.integer_polygon_audit.passed
```

数量证明与完整制造验证范围不同。额外桥形状要求、未使用金属毛刺、曲率工艺下限和完整制造 DRC 仍需专门验证。报告保留这些范围说明。当前工程不宣称可无条件求解任意 GDS 的连续全局最优布局。

## 电极分布评价

页面提供电极数、最近邻间距变异系数、最大未覆盖距离、覆盖率曲线和距离热图。点击指标可查看公式，公式用本地 KaTeX 渲染。

评价只读取最终布局，不改变电极位置或路线。目标圆与电极允许位置圆是两个独立参数。目标圆心默认 `(0, 0)`；需要时点击“取当前结构圆心”。比较多个布局时，固定目标圆、覆盖距离和计算精度，同时报告电极数。

| 指标 | 通俗含义 |
| --- | --- |
| 电极数 N | 布置了多少个电极，同时显示目标圆内数量 |
| 最近邻 CV | 每个电极到最近其他电极的距离变化有多大；越小表示这些距离越接近 |
| 最大未覆盖距离 h | 目标区域中离最近电极最远的位置有多远 |
| 覆盖率 C(ℓ) | 目标区域中，到最近电极中心距离不超过 ℓ 的面积比例 |

目标圆包含没有支撑的空白区域。目标圆外的已布置电极也参与最近距离计算。CV 至少需要两个电极；没有电极时 h 无界。热图与覆盖率显示数值估计和误差区间，不能直接解释为生物信号的实际采集范围。

## 命令行与源码入口

常规使用推荐网页完整 Pad 模式。以下命令为高级分析入口，从项目根目录执行。

只进行几何提取：

```powershell
python -B gds_frontend/run_frontend.py --no-route
```

默认读 `data/`，写 `outputs/geometry/`。这个入口保留早期几何演示逻辑；不加 `--no-route` 时也不等同于网页完整 Pad 流程。

构造一例并生成独立容量证书：

```powershell
python -B gds_frontend/certified_maximum_pipeline.py --input data/C_open_petal_mesh.gds --output-dir outputs/proof/C_open_petal_mesh --center-spacing-um 70 --support-layer 10 --support-datatype 0
```

复核该证书及绑定的 GDS：

```powershell
python -B gds_frontend/verify_proof_result.py --proof-result outputs/proof/C_open_petal_mesh/proof_result.json
```

证明流水线目前使用默认物理规则和 Pad 尺寸，允许选择中心距及支撑层。它有独立的搜索设置，与网页的一次运行不保证相同结果。使用网页自定义尺寸得到的报告不能自动当作此 CLI 的有效输入。命令名不保证 L = U，须读取实际证明结果。默认不传求解时间限值。

| 源码 | 职责 |
| --- | --- |
| `web_app/server.py` | HTTP、文件发现、参数校验、队列和结果接口 |
| `workspace_paths.py` | 可移动的输入、输出路径 |
| `frontend.py` / `island_router.py` | GDS 读取、几何拓扑、导航图和候选路线 |
| `attachment_placement.py` / `center_compaction.py` | 圆岛附着位置与中心优化 |
| `joint_path_flow.py` / `joint_port_augment.py` | 联合路线搜索和出口增广 |
| `pad_router.py` / `ordered_pad_fanout.py` / `pad_sizing.py` | 四边 Pad、外部连接和尺寸策略 |
| `outer_exit_policy.py` | 原支撑外端出口规则 |
| `curved_centerline.py` | 曲线中心线提案与采样 |
| `exact_gds_audit.py` / `capacity_bounds.py` | 导出整数几何审计与容量界 |
| `distribution_metrics.py` | 固定目标区域中的分布指标 |

## 验证

核心检查只使用合成结构，不重跑提供的 9 例布线：

```powershell
python -B tests/verify_package.py
python -B gds_frontend/web_app/verify_parallel_batch.py
python -B gds_frontend/verify_no_time_budget.py
python -B gds_frontend/verify_navigation_geometry.py
```

`verify_package.py` 将项目复制到临时目录，从不同工作目录启动，用随机空闲端口检查输入、前端资源、真实合成任务、输出下载和重启读取结果。临时输入、结果和服务在检查结束后清理。它不启动 8767、8768 或默认 8769。

如果需要真实浏览器检查，先安装可选依赖：

```powershell
python -m pip install -r requirements-dev.txt
python -m playwright install chromium
python -B tests/verify_package.py --ui
```

浏览器检查覆盖主页面、`/admin`、电极数、缩放平移和公式渲染。其他 `verify_*.py` 保留在对应模块附近。其中 `verify_proof_result.py` 是证书复核工具；`verify_frontend.py` 会对输入目录执行较多几何检查。使用 `--help` 查看带参数的入口。

本次整理的依赖验证环境：

```text
Python 3.13.9; NumPy 2.4.1; SciPy 1.17.1; Shapely 2.1.2;
gdstk 1.0.0; pyclipper 1.4.0; NetworkX 3.6.1;
scikit-image 0.26.0; OpenCV 4.13.0.92; Matplotlib 3.10.8.
```

`requirements.txt` 指定兼容范围，不是完整的依赖锁文件。

## 故障处理

| 现象 | 检查与处理 |
| --- | --- |
| `ModuleNotFoundError` | 用启动服务的同一个 Python 安装 `requirements.txt` |
| 端口已占用 | 用 `--port 8770` 等空闲端口启动 |
| 文件列表为空 | 检查当前输入目录；GDS 必须在顶层；点击刷新 |
| 页面无预览 | 查看提示和浏览器控制台；确认静态资源与矢量接口可读取 |
| 支撑层错误 | 在图层列表中确认原始支撑材料，不用电极或金属层代替 |
| 导航拓扑错误 | 阅读 `navigation_diagnostics.json`；检查极窄颈部、点接触和本次工艺尺寸 |
| 电极数少或为 0 | 阅读候选、出口、间距和审计记录；本次构造数量不构成连续不可行证明 |
| 任务运行较久 | 查看阶段和终端输出；不限时搜索可能耗时；减少同时运行数量以降低资源竞争 |
| 重启后显示已中断 | 本工具不自动恢复运行中的求解；重新提交任务 |
| 数量显示但 L 未定义 | 检查整数审计失败原因；回读数量与认证下界范围不同 |

## 发布内容与第三方资源

本目录包含源码、必要验证、输入 GDS 和使用说明。历史研究材料与运行数据不在发布目录内。输入 GDS 是项目数据，未被 `.gitignore` 排除。

KaTeX 0.19.0 的资源随项目提供。其 MIT 许可证在 [`vendor/katex/LICENSE`](gds_frontend/web_app/vendor/katex/LICENSE)，来源与包校验值在 [`provenance.json`](gds_frontend/web_app/vendor/katex/provenance.json)。Python 依赖保留各自的许可证。

主项目及输入 GDS 的发布许可证尚未指定。本次只整理文件并提供忽略规则，没有初始化 Git、配置远程地址或上传仓库。
