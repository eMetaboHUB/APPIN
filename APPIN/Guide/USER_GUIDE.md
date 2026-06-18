<div align="center">

# 🧪 APPIN

### Automated Peak Picking and INtegration for 2D NMR

**User Guide v3.0**

[![R](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://cran.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-App-green.svg)](https://shiny.posit.co/)
[![License](https://img.shields.io/badge/License-CeCILL%202.1-blue.svg)](https://cecill.info/licences/Licence_CeCILL_V2.1-en.html)

*A comprehensive solution for metabolomics 2D NMR data analysis*

[Getting Started](#-getting-started) •
[Workflow](#-workflow) •
[ID Conventions](#-id-naming-conventions) •
[Troubleshooting](#-troubleshooting)

---

</div>

## 📖 Overview

**APPIN** is a Shiny application designed for loading, visualizing, processing, and exporting 2D NMR spectra from Bruker instruments. It provides a complete interactive interface optimized for metabolomics research, from raw Bruker files to publication-ready intensity tables.

### ✨ Key Features

| Feature | Description |
|---------|-------------|
| 🔬 **Multi-format support** | TOCSY, COSY, HSQC, and UFCOSY experiments |
| 🎯 **Two detection methods** | Local maxima with DBSCAN clustering OR CNN neural network |
| 🧠 **Smart CNN detection** | Deep learning model auto-selected per spectrum type |
| ✏️ **Interactive editing** | Click-based or form-based add / edit / delete / fuse |
| 🔄 **Pending changes workflow** | Preview modifications before applying them |
| 📐 **Advanced integration** | Sum (AUC), Gaussian fit, Voigt fit with R² fallback |
| 📈 **Quality diagnostics** | R² distribution, residuals, fit visualization |
| 💾 **Flexible I/O** | Import / export peaks, boxes, full sessions (RDS) |
| ⚡ **Batch processing** | Multi-spectra export with dynamic shift recentering |

---

## 🚀 Getting Started

### Prerequisites

Before installing APPIN, ensure you have:

- **R** (version ≥ 4.0) — [Download from CRAN](https://cran.r-project.org/)
- **RStudio** (recommended) — [Download RStudio Desktop](https://posit.co/download/rstudio-desktop/)

> 💡 **CNN support:** The CNN peak detection requires `tensorflow` and `keras` R packages. They are installed automatically on first launch, but the initial setup may take a few minutes.

### Installation

#### Option A: Direct Download

1. Visit the [GitHub repository](https://github.com/JulienGuibertTlse3/APPIN)
2. Click the green **"Code"** button → **"Download ZIP"**
3. Extract the archive to your preferred location

#### Option B: Git Clone

```bash
git clone https://github.com/JulienGuibertTlse3/APPIN.git
cd APPIN
```

### Launching the Application

1. Open **RStudio**
2. Open the `run_app.R` file
3. Click **"Source"** or press `Ctrl+Shift+Enter`

> 💡 **Note:** The script automatically installs all required packages on first run. Expect a longer startup the first time.

---

## 🗂️ Project Structure

```
APPIN/
│
├── run_app.R                 # 🚀 Entry point — run this file
├── Shine.R                   # Main application
│
├── R/                        # Shiny modules (UI + reactive logic)
│   ├── utils.R
│   ├── mod_load_data.R
│   ├── mod_peak_picking.R    # Local Max + CNN detection
│   ├── mod_manual_editing.R  # Click modes, fuse, delete, edit, add
│   ├── mod_integration.R
│   └── mod_save_export.R
│
├── Function/                 # Core scientific functions
│   ├── Read_2DNMR_spectrum.R # Bruker data reading
│   ├── Vizualisation.R       # Contour plots & DBSCAN
│   ├── Peak_picking.R        # Local maxima detection
│   ├── Peak_fitting.R        # Gaussian / Voigt fitting
│   │
│   ├── CNN_shiny.R           # 🧠 CNN entry point
│   ├── CNN_model.R           # CNN architecture & weights loading
│   ├── CNN_detection.R       # Peak detection (batch + sequential)
│   ├── CNN_filtering.R       # Peak filtering
│   ├── CNN_clustering.R      # DBSCAN clustering & bounding boxes
│   └── CNN_main.R            # Main CNN pipeline function
│
├── saved_model/              # 🧠 CNN pre-trained weights
│   ├── weights/              # Generic homonuclear (TOCSY/COSY/UFCOSY)
│   └── weights_hsqc/         # HSQC-specific model
│
└── www/                      # Web assets (CSS, JS)
```

---

## 🧭 Application Layout

The application has two main tabs accessible from the left sidebar:

| Tab | Purpose |
|-----|---------|
| 📖 **Guide** | Built-in introduction and contextual help |
| 📊 **Analysis** | Main workspace for spectrum visualization and analysis |

The **Analysis** tab is split in two:

- **Left Panel (collapsible accordion)** — Six numbered sections matching the workflow steps
- **Right Panel (tabbed visualization)** — Three tabs:

| Right Panel Tab | Contents |
|-----------------|----------|
| 🌈 **Spectrum** | Interactive Plotly plot with contours, boxes, and centroids |
| 📋 **Data** | Tables of peaks, boxes, and pending changes |
| 📈 **Fit Quality** | R² distribution, fitted boxes details, 2D fit view, residuals |

> 💡 **Tip:** Only one accordion section can be open at a time. This keeps the interface clean and guides you through the workflow in order.

---

## 📋 Workflow

The workflow follows a linear path through the six numbered accordion sections, with pending changes acting as a safety buffer between steps 4 and 5.

### Step 1: 📂 Load Data

<details open>
<summary><b>Loading Bruker spectra into APPIN</b></summary>

#### Bruker File Requirements

Your Bruker data must follow this structure:

```
/<experiment_folder>/<expno>/
    ├── acqus          # Acquisition parameters
    ├── ser  or  fid   # Raw data
    └── pdata/1/
        ├── 2rr        # Processed 2D data
        ├── procs      # Processing parameters (F2)
        ├── proc2s     # Processing parameters (F1)
        └── ...
```

APPIN automatically detects valid Bruker folders by looking for `acqus` + (`ser` or `fid`).

#### Loading Steps

1. Navigate to the **📊 Analysis** tab
2. Expand **"📂 1. Load Data"**
3. Click **"Select Folder"** and browse to the parent folder containing your experiments
4. APPIN scans recursively and lists all valid Bruker folders found
5. Check the spectra you want to analyze (or use **✅ All** / **❌ None**)
6. Click **"📥 Load Selected"**

> 💡 **Tip:** Spectra are cached after first load. Reloading the same folder is much faster.

</details>

---

### Step 2: 📈 Plot Settings

<details>
<summary><b>Configuring the contour plot</b></summary>

#### Spectrum Types

Choose the experiment type — this affects default contour thresholds and DBSCAN epsilon values.

| Type | Description | Default threshold |
|------|-------------|-------------------|
| **TOCSY** | Homonuclear correlations through multiple bonds | 80 000 |
| **COSY** | Direct homonuclear correlations (2-3 bonds) | 80 000 |
| **HSQC** | Heteronuclear ¹H-¹³C single-bond correlations | 20 000 |
| **UFCOSY** | Ultrafast COSY acquisition | 30 000 |

> ⚠️ **Important:** Always verify the spectrum type matches your actual experiment before generating the plot. Wrong defaults will give unusable results.

#### Threshold Parameters

| Parameter | Description |
|-----------|-------------|
| **Threshold** | Minimum intensity for contour display (manual value) |
| **Auto** | Automatically calculates threshold from spectrum noise or maximum |
| **Advanced** | Two automatic methods: `% of max` or `Noise × multiplier` |

> ⚠️ **HSQC note:** The `% of max` method can give inconsistent thresholds for HSQC due to varying ¹³C intensities. The app displays a warning suggesting `Noise ×` instead.

Click **"📊 Generate Plot"** to display the spectrum.

</details>

---

### Step 3: 🎯 Peak Picking

<details>
<summary><b>Automatic peak detection — Local Max or CNN</b></summary>

APPIN offers two detection methods, side by side in the interface. Both work independently and can be re-run as many times as needed — each run replaces the previous results.

#### Method Comparison

| Method | Description | Best For |
|--------|-------------|----------|
| 🟢 **Local Max** | Local maxima detection with optional DBSCAN clustering | Fast detection, well-resolved spectra |
| 🟡 **CNN** | Convolutional neural network with sliding window scan | Complex spectra, overlapping peaks, t1 noise |

> ⚠️ **Important:** Running Local Max after CNN (or vice versa) overwrites the previous results. There is no automatic merging between methods.

---

#### 🟢 Local Max Method

Click **"Local Max"** to run traditional detection. Uses DBSCAN clustering by default to group nearby maxima into single peaks.

**Options (in ⚙️ Options collapsible):**

| Option | Effect |
|--------|--------|
| **No clustering** | Each local maximum becomes a separate peak (no DBSCAN) |
| **Epsilon (eps)** | Maximum distance between points in a cluster (lower = more peaks) |
| **Delete ranges** | Exclude peaks from specific ppm regions (e.g., water, solvent) |

**Recommended epsilon values (set automatically when changing spectrum type):**

- TOCSY / HSQC / COSY: `0.0068`
- UFCOSY: `0.014`

**Delete ranges format:** Comma-separated coordinate pairs, semicolon-separated:

```
0.5,-0.5; 1,0.8; 1.55,1.45; 5.1,4.6;
```

This excludes the F2 ranges `[-0.5, 0.5]`, `[0.8, 1]`, etc. Useful for water suppression artifacts (typically `5.1,4.6`) or known impurities.

---

#### 🟡 CNN Method (Neural Network)

Click **"CNN"** to run deep learning-based detection. The model is automatically selected based on the spectrum type:

| Spectrum type | Model used |
|---------------|-----------|
| TOCSY / COSY / UFCOSY | `saved_model/weights` (generic homonuclear) |
| HSQC | `saved_model/weights_hsqc/weights` (HSQC-specific) |

The CNN is particularly effective for:

- Spectra with overlapping peaks
- Low signal-to-noise conditions
- Complex multiplet patterns
- TOCSY spectra with t1 noise artifacts (vertical streaks)

**CNN Parameters (in 🧠 CNN Parameters sub-section):**

| Parameter | Description | Default |
|-----------|-------------|---------|
| **Prediction threshold** | Minimum CNN confidence to detect a peak (0–1) | 0.3 |
| **Trace filter** | Removes t1 noise artifacts along F2 lines (% of line max) | 50% |

##### How CNN Detection Works

1. **Spectrum normalization** — Normalized using the 99.9th percentile to avoid outlier scaling issues
2. **Sliding window scan** — The CNN scans rows and columns using 2048-point windows with 256-point overlap
3. **Peak classification** — Each point is classified as background / peak edge / peak center
4. **DBSCAN clustering** — Detected points are grouped into peak clusters (same algorithm as Local Max)
5. **CNN-cluster matching** — Only clusters containing a CNN-detected peak are kept (50% adaptive margin)
6. **Bounding box generation** — Each kept cluster gets a bounding box, renamed `cnn_1`, `cnn_2`, ...

> 💡 **HSQC tip:** For HSQC, an additional intensity filter automatically removes peaks below 3% of the maximum intensity to reduce false positives.

##### Progress Indicator

When running CNN detection, a progress bar displays the current step:

| Step | Progress |
|------|----------|
| Preparing spectrum (normalization) | 5% |
| Running neural network | 20% |
| Clustering contour data | 45% |
| Filtering by CNN detections | 70% |
| Updating plot | 90% |
| Complete | 100% |

##### CNN Tuning Tips

| Goal | Action |
|------|--------|
| Detect more peaks | ↓ Lower prediction threshold (e.g., 0.2) |
| Reduce false positives | ↑ Higher prediction threshold (e.g., 0.5) |
| Remove t1 noise traces | ↑ Higher trace filter (e.g., 70%) |
| Keep weak correlations | ↓ Lower trace filter (e.g., 30%) |

> 💡 **Note:** CNN detection typically takes longer than Local Max (10–60 s for typical spectra), but usually produces cleaner results on complex matrices.

</details>

---

### Step 4: ✏️ Manual Editing

<details>
<summary><b>Click modes, fuse, delete, edit, and manual add</b></summary>

This step contains six collapsible sub-sections plus the Apply / Discard buttons at the bottom. **All operations except Fusion and Delete-by-lasso are staged as pending changes** — see Step 5 for the validation workflow.

#### 🖱️ Click Mode

Choose one of four mutually exclusive modes:

| Mode | Action |
|------|--------|
| **Off** | Clicks have no effect (default) |
| **Add peak (1 click)** | Click once on the spectrum to add a centroid → staged as pending |
| **Add box (selection)** | Use the rectangular **Box Select** tool from the Plotly toolbar to draw a box → staged as pending |
| **Delete box on click** | Click inside an existing box to mark it for deletion → staged as pending |

> 💡 **Add box mode:** This mode does NOT use a 2-click sequence. Instead, activate **Box Select** in the Plotly toolbar (square icon at the top right of the plot), then drag a rectangle on the spectrum. The drawn rectangle becomes a new box.

> 💡 **Delete-on-click with overlapping boxes:** If you click in a region covered by multiple overlapping boxes, the smallest one is selected.

#### 🔗 Fusing Peaks and Boxes

Merge multiple peaks into one weighted centroid, and combine their boxes into a single enclosing box.

1. In the **Spectrum** tab, use Plotly's **Box Select** or **Lasso Select** (toolbar icons) to enclose the peaks you want to fuse
2. Return to **Manual Editing**
3. Click **"🔗 Fuse Selected"**

**Result:**

- The selected peaks are removed and replaced by a single fused centroid named `fused_point<N>`
- Position = volume-weighted average (or simple average if no volumes)
- Volume = sum of source volumes
- Their enclosing boxes are merged into `bbox_fused_point<N>` (min/max of all corners)

> ⚠️ **Note:** Fusion modifies your data **immediately** — not via the pending workflow. The fused peak is logged in the "Pending Fusions" table for traceability only. Clicking **Discard** at the bottom does not undo a fusion.

#### 🗑️ Deleting Peaks and Boxes

Bulk delete via lasso/box selection — mirror function of Fusion.

1. In the **Spectrum** tab, use **Box Select** or **Lasso Select** to enclose the peaks to delete
2. Return to **Manual Editing**
3. Click **"🗑️ Delete Selected"**

> ⚠️ **Same as Fusion:** Deletion is applied immediately. The deleted items are logged in "Pending Deletions" for traceability only. Cannot be undone via Discard.

> 💡 **Alternative — table-based delete:** In the **Data** tab, select rows via `Ctrl+Click` (multi-selection) and use the "Delete Selected" button below the table. This staging method goes through the pending workflow and IS reversible via Discard.

#### 📦 Edit Selected Box

Select a row in the boxes table (Data tab) to load it for editing. A dashed green preview appears on the plot.

| Control | Function |
|---------|----------|
| **xmin / xmax / ymin / ymax** | Enter exact coordinates in ppm |
| **↑ ↓ ← →** | Move box by the Step value |
| **+ / −** | Expand or shrink box by Step in all directions |
| **Step** | Movement increment (default 0.01 ppm) |

The preview updates in real time. Click **"Apply Edit"** to stage the change in pending boxes.

> ⚠️ **Validation:** APPIN rejects boxes where xmin ≥ xmax, ymin ≥ ymax, or any side < 0.001 ppm. An error notification will appear.

#### ➕ Add Manually

Form-based addition (no clicking required).

| Type | Inputs | Generated ID |
|------|--------|--------------|
| **Peak** | F2 (ppm), F1 (ppm) | `man1`, `man2`, ... |
| **Box** | xmin, xmax, ymin, ymax | `manual_box1`, `manual_box2`, ... |

Both are staged as pending changes.

</details>

---

### Step 5: 🔄 Pending Changes — Apply or Discard

<details>
<summary><b>The safety buffer for your manual edits</b></summary>

Most manual operations (add peak, add box, edit box, delete box) **do not immediately modify** your data. Instead they go into pending tables, visible in the **Data** tab, so you can review everything before committing.

#### Why Pending Changes?

This staging mechanism protects you from accidental modifications and lets you:

- Review your accumulated changes before applying them
- Discard all pending changes at once if you made a mistake
- Cancel selected pending items individually

#### What goes through pending vs. not?

| Operation | Goes through pending? |
|-----------|----------------------|
| Add peak by click (mode `add_peak`) | ✅ Yes |
| Add box by selection (mode `box_select`) | ✅ Yes |
| Delete box on click (mode `delete_click`) | ✅ Yes |
| Manual add (peak or box via form) | ✅ Yes |
| Edit selected box (coordinates, move, resize) | ✅ Yes |
| Delete peaks/boxes from tables (Ctrl+Click + Delete Selected) | ✅ Yes |
| **Fuse Selected (lasso)** | ✅ Yes  |
| **Delete Selected (lasso)** | ✅ Yes |

#### How to Apply or Discard

At the bottom of the Manual Editing section:

| Button | Effect |
|--------|--------|
| **✅ Apply** | Commits all pending changes to the main data, refreshes the plot |
| **❌ Discard** | Cancels all pending changes (peaks, boxes, fusions, deletions logs) |

> 💡 **Selective discard:** In the Data tab, you can also select specific pending rows and click "Discard Selected" to cancel only those items.

#### Pending Tables in the Data Tab

The **Data** tab shows the current state plus four pending tables:

| Table | Contents |
|-------|----------|
| **Pending Centroids** | Peaks to add or delete (column `status`) |
| **Pending Boxes** | Boxes to add, edit, or delete (column `status`) |
| **Pending Fusions** | Log of recent fusions (already applied) |
| **Pending Deletions** | Log of recent lasso deletions (already applied) |

A yellow warning bar above the plot also displays a summary like "⏳ Pending: 3 peaks, 2 boxes" whenever you have unapplied changes.

</details>

---

### Step 6: 📐 Integration

<details>
<summary><b>Computing peak intensities</b></summary>

Compute the integrated intensity (volume) for each box, using one of three methods.

#### Integration Methods

| Method | Description | Best For |
|--------|-------------|----------|
| **Sum (AUC)** | Direct summation of intensities inside the box | Fast, robust, no parameter tuning |
| **Gaussian** | 2D Gaussian peak fitting | Symmetric peaks, well-resolved |
| **Voigt** | Pseudo-Voigt (Gaussian + Lorentzian mix) | Asymmetric peaks, broad shoulders |

Sum and Peak Fitting are mutually exclusive — selecting one automatically deselects the other.

#### Fitting Options (Gaussian / Voigt only)

| Option | Description | Default |
|--------|-------------|---------|
| **Include R² in export** | Add quality metrics to output CSV | ✅ Enabled |
| **Min R² threshold** | Below this R², the peak falls back to Sum integration | 0.85 |

#### R² Fallback Behavior

When you select Gaussian or Voigt, APPIN fits each box independently. If a fit fails (e.g. too few points, multiplet too complex) or has R² below the threshold:

- The box automatically falls back to Sum integration
- Its method column is marked `sum_fit_failed` or `sum_r2_below_<threshold>`
- A warning notification appears summarizing how many peaks fell back

For multiplets (multiple maxima in one box), each peak is fitted separately and summed:

- If all sub-fits succeed → method = `multiplet_fit`
- If at least one sub-fit fails → method = `multiplet_sum`

Click **"▶️ Run Integration"** to process all boxes. The summary appears below the button, and the **Fit Quality** tab fills with diagnostics.

> ⚠️ **Negative intensities:** Some peaks may show negative volumes. This can be expected (CH₂ groups in multiplicity-edited HSQC have inverted phase) or indicate phase / baseline issues. APPIN warns you when this happens. Exports replace negatives with 0 by default.

</details>

---

### Step 7: 📈 Fit Quality (read-only diagnostics)

<details>
<summary><b>Visualizing the quality of your fits</b></summary>

This panel is in the right-hand visualization area (third tab). It only contains data when you have run an integration with Gaussian or Voigt.

| Section | Content |
|---------|---------|
| **Fit Summary** | Statistics grouped by fitting method (n_boxes, mean / median / min / max R²) |
| **Fitted Boxes Details** | Sortable table with R² per box, color-coded |
| **R² Distribution** | Histogram showing the spread of fit quality |
| **2D Fit View** | Experimental contours (black) vs. fitted model (red) for the selected box |
| **Residuals** | Histogram of `observed − fitted` for the selected box |

To inspect a specific peak, click its row in **Fitted Boxes Details**. The 2D Fit View and Residuals update accordingly.

#### R² Interpretation Guide

| R² Value | Quality | Recommendation |
|----------|---------|----------------|
| **> 0.95** | ✅ Excellent | Use the fit as-is |
| **0.90 – 0.95** | ✅ Very good | Use the fit |
| **0.80 – 0.90** | ⚠️ Good | Verify visually in 2D Fit View |
| **0.70 – 0.80** | ⚠️ Moderate | Consider adjusting the box or switching to Sum |
| **< 0.70** | ❌ Poor | Fallback already triggered (default 0.85 threshold) |

</details>

---

### Step 8: 💾 Save & Export

<details>
<summary><b>Persisting your work and exporting results</b></summary>

The last accordion section bundles four sub-modules:

#### 💼 Session Management (RDS format)

| Action | Description |
|--------|-------------|
| **💾 Save** | Save complete session (peaks, boxes, pending changes, fit results, UI parameters) as `.rds` |
| **📂 Load** | Restore a previously saved session |

> 💡 **What's saved:** Data and parameters only. The original spectra path is recorded but spectra are **not** included in the RDS — you must reload them after restoring a session.

#### 📥 Import (CSV)

| Format | Required columns |
|--------|------------------|
| **Peaks CSV** | `stain_id`, `F2_ppm`, `F1_ppm`, `Volume` |
| **Boxes CSV** | `stain_id`, `xmin`, `xmax`, `ymin`, `ymax` (`Volume` optional) |

APPIN tries `;` separator first (French Excel format), then `,` if only one column is detected. Decimal commas are auto-converted to dots.

#### 📤 Export (CSV)

| Export | Description |
|--------|-------------|
| **Peaks** | Centroids as CSV |
| **Boxes** | Boxes with integrated intensities (current spectrum) |
| **📤 Batch Export (all spectra)** | Apply boxes to **all loaded spectra**, output one intensity column per spectrum |

All exports use `;` separator (compatible with French Excel). Negative intensities are replaced with 0.

#### 🔄 Dynamic Shift Recentering (Batch Export)

The Batch Export now supports **per-box dynamic recentering** to compensate for chemical shift variations between samples (pH, temperature, ionic strength differences).

| Parameter | Description |
|-----------|-------------|
| **Shift Tolerance (ppm)** | Half-width of the search window around each box (slider, 0 – 0.1 ppm) |

##### How It Works

For each box and each spectrum independently:

1. APPIN expands the box by `±tolerance` in both F2 and F1 dimensions
2. The local maximum is found within this expanded window
3. The box is recentered on this maximum, **but only if** this new position captures more intensity than the original
4. Integration is then performed at the (possibly recentered) position

This conservative criterion ("only recenter if it improves intensity") avoids drifting onto noise or neighboring peaks.

##### When to Use

| Scenario | Recommendation |
|----------|----------------|
| Stable peak positions across samples | Tolerance = 0 (disabled) |
| Variable pH / temperature | 0.01 – 0.03 ppm |
| Complex matrices, biological fluids | 0.03 – 0.05 ppm |
| Highly variable conditions | Up to 0.1 ppm (use with caution) |

> 💡 **Best practice:** Start with 0.01 ppm and increase only if you see systematic intensity drops in some samples for known stable peaks.

#### 🔄 Reset All

A single button at the bottom clears everything: peaks, boxes, fits, pending changes, caches, and UI parameters. **Use with care** — there is no confirmation dialog.

#### CSV File Formats

**Peaks:**
```csv
stain_id;F2_ppm;F1_ppm;Volume
peak1;1.234;3.456;123456
```

**Boxes (single spectrum):**
```csv
stain_id;F2_ppm;F1_ppm;xmin;xmax;ymin;ymax;Intensity_<spectrum>
peak1;1.234;3.456;1.200;1.268;3.400;3.512;123456
```

**Batch export (multi-spectra):**
```csv
stain_id;F2_ppm;F1_ppm;xmin;xmax;ymin;ymax;Intensity_<spec1>;Intensity_<spec2>;...
peak1;1.234;3.456;1.200;1.268;3.400;3.512;123456;120345;...
```

> 💡 **Batch processing capacity:** For optimal performance, limit to ~25 TOCSY or 50 COSY/HSQC spectra per batch. Larger batches may run but become slow.

</details>

---

## 🏷️ ID Naming Conventions

Each peak and box has a unique `stain_id`. Understanding the prefix tells you where it came from:

| Prefix | Source | Example |
|--------|--------|---------|
| `peak<N>` | Detected by Local Max | `peak1`, `peak42` |
| `cnn_<N>` | Detected by CNN | `cnn_1`, `cnn_17` |
| `man<N>` | Manually added via form (peak) | `man1`, `man2` |
| `manual_box<N>` | Manually added via form (box) | `manual_box1` |
| `click<N>` | Added by click in `add_peak` mode | `click1`, `click2` |
| `sel_box_<HHMMSS>` | Drawn via Box Select on the spectrum | `sel_box_143052` |
| `fused_point<N>` | Result of fusing multiple peaks | `fused_point1` |
| `bbox_fused_point<N>` | Enclosing box of a fused peak | `bbox_fused_point1` |

> 💡 **Tip:** When importing externally generated CSVs, you can use any string as `stain_id`. Just make sure they are unique within the file — APPIN will auto-rename duplicates with `_dup1`, `_dup2`, ... suffixes.

---

## 💡 Tips & Best Practices

### Workflow Optimization

1. **Start with a QC sample** — Optimize parameters on your most intense, well-resolved spectrum before batch analysis. The same parameters often transfer well to the rest.

2. **Use Auto threshold first** — Then refine manually if needed. The automatic estimates are usually within 20% of the optimal value.

3. **Review the Pending tables before Apply** — A quick scroll through the pending boxes/centroids in the Data tab catches misclicks before they're committed.

4. **Validate with Fit Quality** — When using Gaussian or Voigt, always check the R² distribution histogram. A bimodal distribution often signals two populations of peaks (clean singlets vs. multiplets).

5. **Use lasso operations for bulk cleanup** — Fuse Selected and Delete Selected via Plotly's selection tools are much faster than clicking individually.

6. **Enable shift recentering for biological samples** — Even at 0.01 ppm, dynamic recentering noticeably improves batch consistency for serum / urine / cell extract studies.

7. **Save your session regularly** — Especially after manual editing. An RDS save takes 1-2 seconds and protects against accidental Reset or crashes.

### Parameter Tuning at a Glance

| Goal | Action |
|------|--------|
| Group more peaks together | ↑ Increase epsilon |
| Separate overlapping peaks | ↓ Decrease epsilon |
| Reduce noise detection | ↑ Increase threshold |
| Detect weak signals | ↓ Decrease threshold |
| Stricter CNN detection | ↑ Increase prediction threshold |
| Remove TOCSY t1 noise streaks | ↑ Increase trace filter to 60-80% |
| Allow more peak position drift in batch | ↑ Increase shift tolerance |

---

## 🔧 Troubleshooting

### Data loading

| Problem | Solution |
|---------|----------|
| **No spectrum detected in folder** | Verify Bruker folder structure (`acqus` + `ser`/`fid` + `pdata/1/2rr`) |
| **Wrong ppm range displayed** | Check that `procs` / `proc2s` files are present and not corrupted |
| **Spectra load slowly the first time** | Normal — they are cached in memory for subsequent reads |

### Plot generation

| Problem | Solution |
|---------|----------|
| **Empty or error plot** | Adjust threshold (try Auto), check spectrum type matches the experiment |
| **Plot generation very slow** | Reduce contour count in Advanced settings (Plot Settings) |
| **Wrong axis range** | Click "Autoscale" in the Plotly toolbar to reset |
| **Axes look inverted (numbers)** | Normal — APPIN follows NMR convention (chemical shift decreases left-to-right) |

### Peak detection

| Problem | Solution |
|---------|----------|
| **Too many false positives (Local Max)** | ↑ Increase threshold, ↓ decrease epsilon, or use Delete Ranges for known artifacts |
| **Too few peaks detected (Local Max)** | ↓ Decrease threshold, ↑ increase epsilon |
| **CNN model not loaded** | Verify `saved_model/weights` (or `weights_hsqc/weights` for HSQC) exists |
| **CNN very slow on first run** | TensorFlow initialization (1–2 min). Subsequent runs are fast |
| **CNN finds no peaks** | ↓ Lower prediction threshold (e.g., 0.2) |
| **CNN: too many false positives** | ↑ Increase prediction threshold (0.4–0.6) and trace filter (60–80%) |
| **CNN: TOCSY shows vertical streaks** | ↑ Increase trace filter to 70–80% to remove t1 noise |

### Manual editing

| Problem | Solution |
|---------|----------|
| **Box Select tool doesn't create a box** | Confirm click mode is set to "Add box (selection)" AND the Plotly Box Select tool is active in the toolbar |
| **Click-add-peak doesn't register** | Confirm click mode is "Add peak (1 click)" — modes are mutually exclusive |
| **Pending changes don't appear in plot** | Pending changes are previewed in the Data tab but not on the plot until you Apply |
| **Fuse / Delete Selected does nothing** | Make a selection on the spectrum first (Box Select or Lasso Select), then click the button |
| **Box preview disappears after editing** | Normal — the preview clears after Apply or when you deselect the row |

### Integration

| Problem | Solution |
|---------|----------|
| **All peaks marked as `sum_r2_below_*`** | Threshold too high for your data, or boxes too tight. Lower R² threshold or expand boxes |
| **Many `sum_fit_failed`** | Boxes contain too few points. Expand them or use Sum directly |
| **Negative intensities** | Expected for CH₂ in HSQC multiplicity-edited; otherwise check phase / baseline of the spectrum |
| **Run Integration is slow on Voigt** | Normal — Voigt has more parameters than Gaussian. Use Gaussian if speed matters |

### Save / Load / Batch

| Problem | Solution |
|---------|----------|
| **Import CSV fails** | Check column names match required ones, separator is `;` or `,`, encoding is UTF-8 |
| **Session won't load** | Ensure spectra are reloaded first — session restores annotations only |
| **Batch export gives 0 for some peaks** | Box may fall outside ppm range in some spectra. Check `shift_tolerance_ppm` is appropriate |
| **Shift recentering jumps to wrong peak** | ↓ Reduce tolerance, or shift drift is too large — consider per-spectrum manual review |
| **All intensities are 0 in batch export** | Box coordinates likely outside the spectra's ppm range. Verify with Load Data spectra metadata |

---

## 📚 Additional Resources

- **GitHub Repository:** [JulienGuibertTlse3/APPIN](https://github.com/JulienGuibertTlse3/APPIN)
- **Issue Tracker:** Report bugs and request features via the GitHub issues page
- **Developer Guide:** A separate `Guide_Developpeur_APPIN_v4.0.docx` is available for those who need to extend or maintain APPIN

---

<div align="center">

## 📬 Contact

**Author:** Julien Guibert
**Email:** julien.guibert@inrae.fr

**Project Maintainer:** Marie TREMBLAY-FRANCO
**Email:** marie.tremblay-franco@inrae.fr

**Institution:** INRAe Toxalim / MetaboHUB

---

*APPIN v3.0 — Developed for metabolomics research*

*License: [CeCILL 2.1](https://cecill.info/licences/Licence_CeCILL_V2.1-en.html) — French free software license (GPL-compatible)*

*Last updated: May 2026*

</div>
