# GDS Electrode Layout and Routing Workbench

English · [中文](README.md)

This tool reads a support structure from a GDS file. It places circular electrodes, makes metal routes, and connects each electrode to an external Pad. The browser shows the input, geometry analysis, and output. Each GDS file is one task. Different files can run at the same time.

Version v4 is a separate package made from v2. It keeps the existing solver method and changes the input and output paths. It contains no previous task results, research reports, logs, or caches. All sample cases start without results.

The browser interface currently uses Chinese labels. This guide gives the main label names and their functions.

## Understand the tool

Think of the support structure as a set of curved roads. An electrode is a small metal disk that collects a signal. A wire follows a road to an external connection. A Pad is a larger metal contact for an external device.

One valid connection contains **one electrode, one separate metal net, and one dedicated Pad**. A wire needs support below it. The tool can add a circular support island at an electrode. It can add a support bridge from an allowed outer exit to the Pad area.

### Terms

| Term | Meaning in this project |
| --- | --- |
| GDS | A file that stores layout shapes, cells, references, and layers. It stores geometry, rather than image pixels. |
| Support / substrate | The material geometry read from the selected input layer. It carries the wires. |
| Electrode | A circular metal area that collects a signal. |
| Island | A small circular support area added near an electrode. |
| Wire | A metal area that connects an electrode to a Pad. |
| Net | All metal for one connection. It must stay connected and separate from other nets. |
| Pad | A rectangular metal contact for an external device. |
| Bridge | Added support from an allowed exit to the external Pad area. |
| Layer | A GDS layer number and datatype, such as `10/0`. These numbers do not identify the material. |
| Corridor | A support passage that can carry a wire. |
| Clearance | The shortest distance between two boundaries. It is an edge distance. |
| Lower bound L | The number of complete connections that pass the stated checks. |
| Upper bound U | A number that no valid layout can exceed in the stated model. The bound can be loose. |

### ASD-STE100 writing approach

ASD-STE100 is a standard for controlled technical English. It has writing rules and a controlled dictionary. This guide uses clear writing principles from the standard. It uses short sentences, active instructions, and consistent terms. See the [official introduction](https://www.asd-ste100.org/about_STE.html) and [official FAQ](https://www.asd-ste100.org/STE_faq.html).

For example: “If you want to add an input, put the GDS file in `data/`. Then select Refresh.” The condition comes before the action. The glossary defines project terms before you use them.

This guide has not passed a full controlled-dictionary review. It does not claim certified ASD-STE100 compliance.

## Package structure

```text
.
├── README.md                 # Chinese guide
├── README.en.md              # English guide
├── start_workbench.py        # Recommended start command
├── start_workbench.ps1       # PowerShell wrapper
├── requirements.txt          # Runtime dependencies
├── requirements-dev.txt      # Optional browser-check dependencies
├── .gitignore                # Excludes generated and local files
├── data/                     # Input GDS structures
├── gds_frontend/             # Geometry, placement, routing, and audit code
│   ├── process_rules.json    # Default process rules
│   ├── verify_*.py           # Verification and certificate tools
│   └── web_app/              # Browser pages and HTTP server
│       └── vendor/katex/     # Local formulas, fonts, and license
├── tests/
│   └── verify_package.py     # Package and complete-task check
└── outputs/                  # Created at runtime; not part of the release
    └── runs/                 # Task results and runtime caches
```

The package includes these nine input structures. Their GDS contents match the v2 inputs used for this package.

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

## Install and start

### Requirements

Use Python 3.11 or later. This package was checked on Windows with Python 3.13.9. The Linux and macOS commands below are provided for use, but those systems were not checked in this release preparation.

The frontend uses HTML, CSS, and JavaScript directly. You do not need a Node.js or npm build. The package includes KaTeX. After you install Python dependencies, the pages and formulas can operate locally without an internet connection.

### Windows PowerShell

1. Open PowerShell.
2. Go to the package folder. Replace the example path with your package location.

```powershell
Set-Location -LiteralPath 'A:\Electrode_design\结构\布线布电极v4'
```

3. Make a virtual environment and install the dependencies.

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

If the `py` command is not available, use `python -m venv .venv`. The commands use the virtual environment interpreter directly. You do not need to activate a PowerShell script.

4. Start the server.

```powershell
.\.venv\Scripts\python.exe -B start_workbench.py
```

5. Open the main page: <http://127.0.0.1:8769/>.
6. Open the batch page: <http://127.0.0.1:8769/admin>.

The terminal shows the input folder, output folder, and worker count. Keep the terminal open. Press `Ctrl+C` to stop the HTTP service. Already submitted jobs wait to finish. If you close the browser, the jobs continue.

If your Python environment already has the dependencies, you can start directly:

```powershell
python -B start_workbench.py --port 8769 --workers 2
```

The `start_workbench.ps1` wrapper uses the package `.venv` if it exists. Otherwise, it uses `python` from PATH.

### Linux / macOS

From the package root, run:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -B start_workbench.py
```

### Start options

| Option | Default | Function |
| --- | --- | --- |
| `--host` | `127.0.0.1` | Bind address. The default permits access from this computer. |
| `--port` | `8769` | HTTP port. |
| `--workers` | `2` | Number of simultaneous jobs. Select 1 to 4. |
| `--input-dir` | Package `data/` | Folder of input GDS files. |
| `--output-dir` | Package `outputs/runs/` | Folder for results, status, and caches. |

Default paths use the script location. You can start from a different working folder. Explicit relative paths use the current working folder.

```powershell
python -B start_workbench.py --input-dir 'D:\GDS inputs' --output-dir 'D:\GDS results' --workers 2
```

If you need access from your local network, run:

```powershell
python -B start_workbench.py --host 0.0.0.0 --port 8769
```

Use the server computer's actual IP address on another device. `0.0.0.0` is a bind address. Network routes and firewall settings determine whether that device can connect. The current `/admin` page is a batch console. It has no login or access-control system. Use it locally or on a controlled network.

For direct use of `gds_frontend/web_app/server.py`, the equivalent environment variables are `GDS_WORKBENCH_HOST`, `GDS_WORKBENCH_PORT`, `GDS_WORKBENCH_WORKERS`, `GDS_WORKBENCH_DATA_DIR`, and `GDS_WORKBENCH_RUNS_DIR`. The recommended launcher sets these variables from its command-line options.

## Run one task

1. Select a GDS file in the left list. The page shows the vector support outline, layers, file hash, and geometry facts.
2. Check **用于分析的支撑层** (support layer for analysis). The tool suggests `10/0` first. Otherwise, it suggests the layer with the most polygons. Confirm the material from the file source.
3. Select **连接四边 Pad** (connect four-side Pads). This mode makes complete electrode-to-Pad connections.
4. Set the process parameters. Each field shows its unit.
5. Select **开始分析与试布线** (start analysis and routing). The page shows the stage and queue state.
6. Wait for completion. Read **已布置电极** (placed electrodes), lower bound L, upper bound U, and the audit state.
7. Select the stage cards to inspect the input, anchor domain, corridors and exits, metal layout, and analytic curves.
8. Download `routing.gds` and `summary.json`. Keep `routing.width_witness.json` with the GDS for a wire-width check.

The large electrode count shows connections confirmed by GDS readback. If the integer audit fails, that count is not a certified lower bound L. Check `capacity_interval.integer_polygon_lower_verified` before you use L.

The **仅几何终端试布线** mode routes to geometry terminals. It does not certify physical Pad connections. Its count has a different scope.

### Inspect a view

Use the mouse wheel to zoom. Hold and drag the mouse to move the view. Select **适合窗口** (fit to window) to reset the view.

The input and final layout use vectors from GDS polygons. Their zoom detail does not depend on PNG resolution. The extracted corridor image is a PNG.

The analytic curve card is available when the result contains curve segments. A straight route can have no curve view. GDS output uses polygons. The analytic view shows the centerline; it does not certify a manufacturing bend-radius limit.

## Parameters and units

`1 mm = 1000 μm`. Thus, `0.07 mm = 70 μm`. API and JSON fields with the `_um` suffix use μm.

| Parameter | Default | Effect |
| --- | --- | --- |
| Electrode diameter | 30 μm | Diameter of the metal electrode disk. |
| Minimum center distance | 0.07 mm | Minimum distance between any two electrode centers. |
| Allowed electrode radius | 3 mm | Electrode centers must be inside this disk. Its default diameter is 6 mm. |
| Wire width | 5 μm | Width used for route construction and checks. |
| Net spacing | 4 μm | Edge clearance between different nets. It also sets island spacing. |
| Support margin | 4 μm | Margin from the metal edge to the support edge. |
| Initial frame side | 32 mm | Initial size of the square external Pad area. |
| Pad width along side | 0.5 mm | Pad size parallel to a square side. |
| Pad radial length | 3 mm | Pad size along the inward/outward direction. |
| Pad center pitch | 1 mm | Distance between adjacent Pad centers along a side. |
| Pad size policy | Keep size and expand frame | Makes the frame larger when more Pads are needed. |
| Minimum Pad width / length | 0.5 / 3 mm | Lower dimensions for the shrink policy. |

The placement disk uses the center of the original support's minimum enclosing circle. That center is rounded to the native GDS grid. This definition also applies to translated structures. Wires and Pads can extend outside the disk.

If an organoid diameter is 2–3 mm, you can set the allowed radius to 1–1.5 mm. Enter a radius. The rule checks electrode centers; it does not require the full metal disk to stay inside the placement disk.

Dimensions must be positive. Spacing, support margin, and the extra center-distance rule can be zero. Nets must remain separate with a zero spacing setting. Island and metal dimensions can require a larger effective center distance than the entered value.

Pad pitch must be at least Pad width plus net spacing. Pad dimensions must contain the wire and its margins. The used Pads form a continuous bank on each side. The tool derives the Pad count from demand and geometry.

The shrink policy first reduces a common size factor, within the specified lower dimensions. It then expands the frame if necessary. If the lower dimensions equal the reference dimensions, Pads retain their size. Read the final dimensions, pitch, and frame size in the report.

Complete Pad mode uses allowed outer exits from the original support. The wide-region threshold belongs to geometry-terminal mode. It does not set the exits in complete Pad mode.

## Batch jobs and parallel work

1. Open `/admin`.
2. Select several GDS files, or select **选择所有可运行** (select all available).
3. Check the support layer for each file.
4. Set the shared process rules and Pad dimensions.
5. Select **批量提交所选任务** (submit selected tasks).
6. Read the queue. Open a task result on the main page.

Two jobs run simultaneously by default. Other jobs wait in the queue. `--workers` controls solver jobs, rather than browser requests. At most 32 jobs can be running or queued. An active GDS cannot enter the queue twice. The server validates the entire batch before submission. An invalid layer or parameter prevents that batch from entering the queue.

The default solver has no wall-clock time budget. Complex structures can take a long time. Candidate counts, iteration counts, and other discrete search limits still apply. Unlimited runtime does not search every possible continuous layout. More workers usually need more memory.

## Inputs and outputs

Put original support GDS files at the top level of `data/`. Select Refresh. Both `.gds` and `.GDS` are accepted. The server does not scan subfolders. An input needs support geometry, but no generator source, predefined centerline, or connection metadata.

The reader handles cells, references, and native GDS units. Select the correct support layer. A point contact has no width and cannot carry a finite-width wire. With a 5 μm wire and 4 μm margins on both sides, a straight passage needs approximately 13 μm of support width. Turns and multiple wires can need more space.

Keep generated `routing.gds` files outside the input folder. They contain added islands, bridges, and the outer frame. Using them as a support input changes the problem. If an input changes, a previous task can become stale and need another run.

Each run has a folder: `outputs/runs/<job_id>/`.

| File | Contents |
| --- | --- |
| `status.json` | State, rules, progress, and results. |
| `summary.json` | Full geometry, placement, routing, and bound report. |
| `routing.gds` | Output support, electrodes, metal, and Pads when a construction succeeds. |
| `routing.width_witness.json` | Centerline evidence for the exported wire-width check. |
| `integer_polygon_audit.json` | Integer-grid audit of the exported GDS. |
| `navigation_diagnostics.json` | Navigation inset, topology checks, and component exclusions. |
| `graph.json.gz` / `regions.json.gz` | Compressed graph and vector domains. |
| `graph.png` / `routing.png` | Preview images. |
| `error.log` | Error details if the task fails. |

Some files are absent when no complete route is constructed. Read `output_layers` in the report for actual output layer numbers.

Display and navigation caches also use `outputs/runs/`. `.gitignore` excludes all of `outputs/`. Completed tasks remain available after a restart. Running solves do not resume automatically.

Before you remove results, wait for active jobs to finish and stop the service. Then remove the task folders you no longer need. Removing all of `outputs/` removes both results and caches. The tool creates the folder again when needed.

## Method and result scope

The tool uses one workflow:

```text
GDS + process rules
  → Read support geometry and check topology
  → Build navigation domains, corridors, node windows, and outer exits
  → Propose island attachments and multiple-track routes
  → Select electrodes, routes, exits, and Pads jointly
  → Improve center positions while retaining the connection count
  → Export metal and support GDS
  → Check GDS readback and integer geometry
  → Report the electrode count and capacity interval [L, U]
```

The geometry frontend reads actual GDS polygons. It adapts navigation insets to geometric features. It checks connected components and holes. It records exclusions that have a proof within the stated model. Geometry that cannot be extracted reliably produces a diagnostic error.

The solver first increases the number of complete connections. For equal counts, it prefers positions near the center. Complete candidate routes become integer selection variables. Conflict constraints describe center distance, metal spacing, islands, exits, and Pad occupancy.

Wide support can carry multiple route candidates. Geometry checks determine which tracks can coexist. Attachment search, residual route augmentation, and center rerouting improve the constructed layout.

External bridges follow the outer-exit policy. The bridge and Pad geometry also enter the checks. Curved-centerline proposals use tangent-continuous rounded turns. The exporter samples them into GDS polygons under an error tolerance.

Let N* be the true maximum count in the stated model. Valid bounds satisfy:

$$
L \le N^* \le U.
$$

L is the complete connection count that passes the stated exported-GDS audits. U comes from candidate-independent geometric packing bounds and applicable cut bounds.

If L = U and all proof premises pass, the result proves the maximum **within the stated geometric model**. If L < U, report the constructed count and the open interval. An optimum in a finite candidate library does not prove a continuous maximum.

Read these fields:

```text
result.capacity_interval.lower_bound
result.capacity_interval.upper_bound
result.capacity_interval.integer_polygon_lower_verified
result.capacity_interval.declared_geometric_model_optimality_proven
result.routing.gds_roundtrip_audit.passed
result.routing.integer_polygon_audit.passed
```

A count proof and a complete manufacturing check have different scopes. Extra bridge-shape requirements, unused metal protrusions, manufacturing bend-radius limits, and full design-rule checks need separate validation. The reports record these limits. The package does not claim an unconditional continuous global optimum for every GDS.

## Electrode distribution assessment

The result page shows the electrode count, nearest-neighbor spacing variation, maximum uncovered distance, a coverage curve, and a distance heatmap. Select a metric to open its formulas. Local KaTeX resources render the formulas.

Assessment reads the final layout. It does not move electrodes or change routes. The assessment disk and placement disk are separate settings. The assessment center defaults to `(0, 0)`. Select **取当前结构圆心** to use the current structure center.

For comparisons, keep the target disk, coverage distance, and resolution fixed. Always report the electrode count with the metrics.

| Metric | Plain meaning |
| --- | --- |
| Electrode count N | Number of electrodes, with the count inside the target disk. |
| Nearest-neighbor CV | Variation in each electrode's distance to its closest other electrode. A smaller value means these distances are more similar. |
| Maximum uncovered distance h | The largest distance from a point in the target to its nearest electrode center. |
| Coverage C(ℓ) | Fraction of target area within distance ℓ of an electrode center. |

The full target disk includes gaps without support. Accepted electrodes outside the target can contribute to nearest distances. CV requires at least two electrodes. With no electrodes, h is unbounded.

Coverage and heatmap values include numerical estimates and uncertainty intervals. These geometric distances are not a calibrated biological recording range.

## Command-line tools and code map

Use complete Pad mode in the browser for normal operation. Run the following advanced commands from the package root.

Extract geometry only:

```powershell
python -B gds_frontend/run_frontend.py --no-route
```

This command reads `data/` and writes `outputs/geometry/`. It keeps an earlier geometry demonstration flow. Its route mode does not perform the full browser Pad workflow.

Construct one case and make an independent bound certificate:

```powershell
python -B gds_frontend/certified_maximum_pipeline.py --input data/C_open_petal_mesh.gds --output-dir outputs/proof/C_open_petal_mesh --center-spacing-um 70 --support-layer 10 --support-datatype 0
```

Verify the certificate and its bound GDS files:

```powershell
python -B gds_frontend/verify_proof_result.py --proof-result outputs/proof/C_open_petal_mesh/proof_result.json
```

The certificate pipeline currently uses default physical rules and Pad dimensions. It permits a selected center distance and support layer. Its search settings are separate from a browser run, so the counts can differ. A browser report with custom dimensions is not automatically valid input for this CLI.

The command name does not guarantee L = U. Read the actual proof result. The default commands do not supply a solver time limit.

| Module | Function |
| --- | --- |
| `web_app/server.py` | HTTP, file discovery, parameter checks, queue, and result endpoints. |
| `workspace_paths.py` | Portable input and output locations. |
| `frontend.py` / `island_router.py` | GDS reading, topology, navigation, and route candidates. |
| `attachment_placement.py` / `center_compaction.py` | Island attachments and center improvement. |
| `joint_path_flow.py` / `joint_port_augment.py` | Joint route search and exit augmentation. |
| `pad_router.py` / `ordered_pad_fanout.py` / `pad_sizing.py` | Four-side Pads, external connections, and size policy. |
| `outer_exit_policy.py` | Rules for exits from the original support. |
| `curved_centerline.py` | Curve proposals and sampling. |
| `exact_gds_audit.py` / `capacity_bounds.py` | Exported integer-geometry audits and capacity bounds. |
| `distribution_metrics.py` | Metrics in a fixed target region. |

## Verification

These checks use synthetic geometry. They do not reroute the nine supplied layouts.

```powershell
python -B tests/verify_package.py
python -B gds_frontend/web_app/verify_parallel_batch.py
python -B gds_frontend/verify_no_time_budget.py
python -B gds_frontend/verify_navigation_geometry.py
```

The package check copies the project to a temporary folder. It starts from another working folder and uses a random available port. It checks input reading, assets, a real synthetic task, artifact downloads, and restored results after restart.

Temporary inputs, results, and servers are removed after the check. The check does not start ports 8767, 8768, or the default 8769.

For an actual browser check, install the optional tools:

```powershell
python -m pip install -r requirements-dev.txt
python -m playwright install chromium
python -B tests/verify_package.py --ui
```

This check covers the main page, `/admin`, electrode counts, zoom and pan, and formula rendering. Other `verify_*.py` files remain next to their modules. `verify_proof_result.py` is a certificate verifier. `verify_frontend.py` performs broader geometry checks on its input folder. Use `--help` for tools that accept arguments.

Dependency environment used in release preparation:

```text
Python 3.13.9; NumPy 2.4.1; SciPy 1.17.1; Shapely 2.1.2;
gdstk 1.0.0; pyclipper 1.4.0; NetworkX 3.6.1;
scikit-image 0.26.0; OpenCV 4.13.0.92; Matplotlib 3.10.8.
```

`requirements.txt` gives compatibility ranges. It is not a full dependency lock file.

## Troubleshooting

| Symptom | Action |
| --- | --- |
| `ModuleNotFoundError` | Install `requirements.txt` with the same Python that starts the service. |
| Port already in use | Select an available port, such as `--port 8770`. |
| No files in the list | Check the input folder. Put GDS files at its top level. Select Refresh. |
| No preview | Read the page message and browser console. Check static and vector endpoints. |
| Incorrect support layer | Select the original support material layer. |
| Navigation topology error | Read `navigation_diagnostics.json`. Check narrow necks, point contacts, and the process dimensions. |
| Few or zero electrodes | Read candidate, exit, clearance, and audit records. A construction count does not prove continuous infeasibility. |
| Long-running job | Read the current stage and terminal output. Unlimited searches can take time. Reduce concurrent jobs to reduce resource competition. |
| Interrupted after restart | Running solves do not resume. Submit the task again. |
| Visible count but undefined L | Read the integer audit failure. The readback count and certified lower bound have different scopes. |

## Release contents and third-party resources

This folder contains runtime code, necessary verification tools, input GDS, and user guides. Historical research and run data are excluded. `.gitignore` keeps the supplied GDS inputs eligible for a future repository.

The package includes KaTeX 0.19.0. Its MIT license is in [`vendor/katex/LICENSE`](gds_frontend/web_app/vendor/katex/LICENSE). Source and integrity data are in [`provenance.json`](gds_frontend/web_app/vendor/katex/provenance.json). Python dependencies keep their respective licenses.

A release license for the project and its input GDS has not been selected. This preparation did not initialize Git, configure remotes, or upload a repository.
