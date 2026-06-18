# DESCRIPTION

```
Package: APPIN
Title: Automated Peak Picking and INtegration for 2D NMR Metabolomics
Version: 2.0.0
Author: Julien Guibert
Maintainer: Julien Guibert <your.email@aviwell.com>
Description: An R/Shiny application for the analysis of 2D NMR spectra in
    metabolomics. APPIN supports HSQC, TOCSY, COSY and UF-COSY spectra.
    It provides peak detection via either local maxima or a convolutional
    neural network (CNN), multiplet clustering with DBSCAN, and peak
    integration with sum, Gaussian or pseudo-Voigt fitting. The package
    also includes manual editing tools, session save/load and CSV/batch
    export with shift-tolerance conflict detection.
License: GPL-3
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
Depends:
    R (>= 4.1.0)
Imports:
    shiny,
    plotly,
    dbscan,
    DT,
    data.table,
    jsonlite,
    minpack.lm,
    Matrix,
    stats,
    utils
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown,
    pkgdown,
    shinytest2,
    covr
VignetteBuilder: knitr
Config/testthat/edition: 3
URL: https://gitlab.com/yourusername/APPIN
BugReports: https://gitlab.com/yourusername/APPIN/issues
Config/roxygen2/version: 8.0.0
```

# `clean_centroids_df`: Clean imported centroids dataframe

## Description

Cleans a dataframe of imported peak centroids by converting string
representations to numeric values (handles European decimal comma).

## Usage

```r
clean_centroids_df(df)
```

## Arguments

* `df`: Data frame with columns: F2_ppm, F1_ppm, Volume

## Value

Data frame with cleaned numeric columns

# `clear_spectrum_cache`: Clear the spectrum cache

## Description

Removes all cached spectra from memory. Useful when memory is low
or when spectra files have been modified on disk.

## Usage

```r
clear_spectrum_cache()
```

## Value

NULL (invisible)

# `clip_negative_box_intensities`: For each bounding box, sum the raw spectrum intensity inside the box. If
the sum is < 0, set the matching peak's intensity column to 0 in `peaks_df`.
Boxes and peaks are kept (coordinates preserved); only the intensity is
forced to 0. This is consistent with the post-integration behaviour in
`mod_integration.R` and the export behaviour in `mod_export.R`.

## Description

For each bounding box, sum the raw spectrum intensity inside the box. If
the sum is < 0, set the matching peak's intensity column to 0 in `peaks_df`.
Boxes and peaks are kept (coordinates preserved); only the intensity is
forced to 0. This is consistent with the post-integration behaviour in
`mod_integration.R` and the export behaviour in `mod_export.R`.

## Usage

```r
clip_negative_box_intensities(peaks_df, boxes_df, spectrum_matrix)
```

## Arguments

* `peaks_df`: Data frame with `stain_id` and an intensity column
(`stain_intensity` for CNN output, `Volume` for local-max output).
* `boxes_df`: Data frame with `stain_id`, `xmin`, `xmax`, `ymin`, `ymax`.
* `spectrum_matrix`: Raw 2D NMR matrix (rownames = F1 ppm, colnames = F2 ppm).

## Value

List with `peaks_df` (clipped) and `n_clipped` (count).

# `compute_f2_territories`: Compute exclusive F2 territory bounds for each box

## Description

Two boxes compete for the same peak only if they also overlap in F1 (same
"row" of the spectrum). For each box, the F2 territory is bounded by the
midpoint to the nearest F1-overlapping neighbour on each side. Recentering
(both the search window and the final shift) is then confined to this
territory, so neighbouring integration boxes can never drift onto the same
peak or overlap after recentering — which lets `tol` be larger without
creating conflicts. Especially important along the TOCSY/COSY diagonal where
boxes are densely packed.

## Usage

```r
compute_f2_territories(boxes)
```

## Value

data.frame with columns terr_lo, terr_hi (F2 ppm bounds) per row of boxes.

# `get_box_intensity`: Calculate intensity within bounding boxes

## Description

Calculates the intensity of NMR signals within specified bounding boxes.
Supports two methods: simple sum (fast) or peak fitting (more accurate).

## Usage

```r
get_box_intensity(mat, ppm_x, ppm_y, boxes, method = "sum", model = "gaussian")
```

## Arguments

* `mat`: Numeric matrix. The spectrum data matrix
* `ppm_x`: Numeric vector. F2 (x-axis) chemical shift values
* `ppm_y`: Numeric vector. F1 (y-axis) chemical shift values
* `boxes`: Data frame. Bounding boxes with columns: xmin, xmax, ymin, ymax
* `method`: Character. Integration method: "sum" (default) or "fit"
* `model`: Character. Fitting model when method="fit": "gaussian" (default) or "voigt"

## Note

When method="fit", requires calculate_fitted_volumes() from Peak_fitting.R

## Value

Numeric vector of intensities, one per box

# `%||%`: Null-coalesce operator

## Description

Returns the first argument if it is not NULL, otherwise returns the second.
Similar to the ?? operator in C# or the // operator in Perl.

## Usage

```r
a %||% b
```

## Arguments

* `a`: First value to check
* `b`: Default value if a is NULL

## Value

a if not NULL, otherwise b

## Examples

```r
NULL %||% "default"  # returns "default"
"value" %||% "default"  # returns "value"
```

# `mod_box_editor_server`: Box Editor Module - Server

## Description

Box Editor Module - Server

## Usage

```r
mod_box_editor_server(id, rv, data_reactives, parent_input, parent_session)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values
* `data_reactives`: List. Reactive expressions
* `parent_input`: Shiny input. Parent input object for table selections
* `parent_session`: Shiny session. Parent session for plotlyProxy

## Value

NULL (side effects only)

# `mod_box_editor_ui`: Box Editor Module - UI

## Description

Box Editor Module - UI

## Usage

```r
mod_box_editor_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing box editing controls

# `mod_click_mode_server`: Click Mode Module - Server

## Description

Click Mode Module - Server

## Usage

```r
mod_click_mode_server(id, rv, data_reactives, peak_picking = NULL)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values (first_click_for_box, last_click_coords,
pending_boxes, pending_centroids, centroids_data)
* `data_reactives`: List. Reactive expressions (bounding_boxes_data, result_data)
* `peak_picking`: List. Return value from mod_peak_picking_server (for eps_value)

## Value

A list containing:

* `box_click_mode`: Reactive returning current click mode

# `mod_click_mode_ui`: Click Mode Module - UI

## Description

Click Mode Module - UI

## Usage

```r
mod_click_mode_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing click mode controls

# `mod_delete_server`: Delete Module - Server

## Description

Delete Module - Server

## Usage

```r
mod_delete_server(id, rv)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values (centroids_data, modifiable_boxes, fixed_boxes, pending_deletions)

## Value

NULL (side effects only)

# `mod_delete_ui`: Delete Module - UI

## Description

Delete Module - UI

## Usage

```r
mod_delete_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing delete controls

# `mod_export_server`: Export Module - Server

## Description

Export Module - Server

## Usage

```r
mod_export_server(id, status_msg, rv, load_data, data_reactives)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `rv`: List. Shared reactive values
* `load_data`: List. Return value from mod_load_data_server
* `data_reactives`: List. Reactive expressions
* `boxes`: data.frame with stain_id, xmin, xmax, ymin, ymax
* `tol_f2`: numeric, tolerance on the F2 (1H) axis in ppm
* `tol_f1`: numeric, tolerance on the F1 (13C) axis in ppm

## Value

A list of reactives exposing:

* `shift_tolerance_ppm()` — current tolerance slider value
* `conflict_ids()` — character vector of stain_id flagged as
potentially overlapping. Use this in the main plot to color those
boxes red (others green).

list with `pairs` (data.frame of conflicting stain_id pairs),
`conflict_ids` (vector), `n_pairs`, `n_boxes_in_conflict`

# `mod_export_ui`: Export Module - UI

## Description

Export Module - UI

## Usage

```r
mod_export_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing export controls

# `mod_fusion_server`: Fusion Module - Server

## Description

Fusion Module - Server

## Usage

```r
mod_fusion_server(id, rv)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values (centroids_data, modifiable_boxes, fixed_boxes, pending_fusions)

## Value

NULL (side effects only)

# `mod_fusion_ui`: Fusion Module - UI

## Description

Fusion Module - UI

## Usage

```r
mod_fusion_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing fusion controls

# `mod_import_server`: Import Module - Server

## Description

Import Module - Server

## Usage

```r
mod_import_server(id, rv, refresh_nmr_plot)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values
* `refresh_nmr_plot`: Function. Function to refresh the NMR plot

## Value

NULL (side effects only)

# `mod_import_ui`: Import Module - UI

## Description

Import Module - UI

## Usage

```r
mod_import_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing import controls

# `mod_integration_server`: Integration Module - Server

## Description

Server logic for the integration module. Handles method selection,
peak integration calculations, and results export.

## Usage

```r
mod_integration_server(id, status_msg, load_data, rv)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `load_data`: List. Return value from mod_load_data_server containing:
* `bruker_data`: Reactive for current spectrum data
* `rv`: List. Named list of reactive values:
* `modifiable_boxes`: ReactiveVal for editable boxes
* `fit_results_data`: ReactiveVal for fit results (will be updated)
* `last_fit_method`: ReactiveVal for last fit method used (will be updated)

## Value

A list containing:

* `effective_integration_method`: Reactive returning the selected method
* `integration_results`: ReactiveVal containing integration results
* `integration_done`: ReactiveVal indicating if integration is complete

# `mod_integration_ui`: Integration Module - UI

## Description

Creates the UI components for the integration section.
Includes method selection (Sum/Gaussian/Voigt), fitting options,
and results display.

## Usage

```r
mod_integration_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing the module's UI elements

# `mod_load_data_server`: Load Data Module - Server

## Description

Internal reactive value storing the list of loaded spectra
Internal reactive value storing the currently displayed spectrum
Parse selected directory path
Reactive that processes the shinyDirChoose selection and returns
a normalized, validated directory path.

## Usage

```r
mod_load_data_server(id, status_msg, trigger_subfolder_update = NULL)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `trigger_subfolder_update`: ReactiveVal. Trigger to notify parent
when spectra list changes (for updating subfolder selector)

## Details

Server logic for the data loading module. Handles directory selection,
spectrum detection, and loading of Bruker NMR data.

## Value

A list containing:

* `spectra_list`: Reactive containing named list of loaded spectra
* `bruker_data`: Reactive containing currently selected spectrum
* `main_directory`: Reactive containing the selected main directory path

Character. Normalized path to selected directory, or NULL if invalid
Render selected directory path
Detect Bruker spectrum subfolders
Reactive that scans the selected directory for valid Bruker NMR
spectrum folders (containing 'acqus' and either 'ser' or 'fid' files).
Character vector of paths to valid spectrum folders
Render UI for available spectra selection
Creates a dynamic UI with checkboxes for each detected spectrum,
along with "Select All" / "Deselect All" buttons and a load button.
Handle "Select All" button click
Handle "Deselect All" button click
Handle spectrum loading
Loads selected Bruker spectra when the load button is clicked.
Shows progress bar and notifications for success/failure.

## Examples

```r
# In server:
load_data <- mod_load_data_server("load_data", status_msg = status_msg)
# Access loaded spectra:
load_data$spectra_list()
```

# `mod_load_data_ui`: Load Data Module - UI

## Description

Creates the UI components for the data loading section.
This includes a directory picker, list of available spectra with
checkboxes, and load button.

## Usage

```r
mod_load_data_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing the module's UI elements

## Examples

```r
# In UI definition:
mod_load_data_ui("load_data")
```

# `mod_manual_add_server`: Manual Add Module - Server

## Description

Manual Add Module - Server

## Usage

```r
mod_manual_add_server(id, rv, data_reactives, peak_picking)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values (centroids_data, modifiable_boxes, pending_*)
* `data_reactives`: List. Reactive expressions (result_data for contour_data)
* `peak_picking`: List. Return value from mod_peak_picking_server (for eps_value)

## Value

NULL (side effects only)

# `mod_manual_add_ui`: Manual Add Module - UI

## Description

Manual Add Module - UI

## Usage

```r
mod_manual_add_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing manual add controls

# `mod_manual_editing_server`: Manual Editing Module - Server (Wrapper)

## Description

Manual Editing Module - Server (Wrapper)

## Usage

```r
mod_manual_editing_server(
  id,
  status_msg,
  load_data,
  rv,
  data_reactives,
  refresh_nmr_plot,
  peak_picking,
  parent_input,
  parent_session
)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `load_data`: List. Return value from mod_load_data_server
* `rv`: List. Named list of reactive values
* `data_reactives`: List. Named list of reactive expressions
* `refresh_nmr_plot`: Function. Function to refresh the NMR plot
* `peak_picking`: List. Return value from mod_peak_picking_server
* `parent_input`: Shiny input. Parent input object for table selections
* `parent_session`: Shiny session. Parent session for plotlyProxy

## Value

A list containing:

* `box_click_mode`: Reactive returning the current click mode

# `mod_manual_editing_ui`: Manual Editing Module - UI (Wrapper)

## Description

Manual Editing Module - UI (Wrapper)

## Usage

```r
mod_manual_editing_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing all manual editing UI components

# `mod_pending_changes_server`: Pending Changes Module - Server

## Description

Pending Changes Module - Server

## Usage

```r
mod_pending_changes_server(
  id,
  rv,
  data_reactives,
  load_data,
  refresh_nmr_plot,
  parent_input,
  parent_session
)
```

## Arguments

* `id`: Character. The module's namespace ID
* `rv`: List. Shared reactive values
* `data_reactives`: List. Reactive expressions
* `load_data`: List. Return value from mod_load_data_server (for bruker_data)
* `refresh_nmr_plot`: Function. Function to refresh the NMR plot
* `parent_input`: Shiny input. Parent input object for table selections
* `parent_session`: Shiny session. Parent session for plotlyProxy

## Value

NULL (side effects only)

# `mod_pending_changes_ui`: Pending Changes Module - UI

## Description

Pending Changes Module - UI

## Usage

```r
mod_pending_changes_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing apply/discard buttons

# `mod_reset_server`: Reset Module - Server

## Description

Reset Module - Server

## Usage

```r
mod_reset_server(id, status_msg, rv, parent_session)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `rv`: List. Shared reactive values
* `parent_session`: Shiny session. Parent session for updating inputs

## Value

A list containing:

* `reset_triggered`: ReactiveVal that increments when reset is triggered

# `mod_reset_ui`: Reset Module - UI

## Description

Reset Module - UI

## Usage

```r
mod_reset_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing reset button

# `mod_save_export_server`: Save & Export Module - Server (Wrapper)

## Description

Save & Export Module - Server (Wrapper)

## Usage

```r
mod_save_export_server(
  id,
  status_msg,
  rv,
  load_data,
  data_reactives,
  refresh_nmr_plot,
  parent_session,
  parent_input
)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `rv`: List. Named list of reactive values
* `load_data`: List. Return value from mod_load_data_server
* `data_reactives`: List. Named list of reactive expressions
* `refresh_nmr_plot`: Function. Function to refresh the NMR plot
* `parent_session`: Shiny session. Parent session for updating inputs
* `parent_input`: Shiny input. Parent input object for reading UI values

## Value

A list containing:

* `reset_triggered`: Reactive that increments when reset is triggered

# `mod_save_export_ui`: Save & Export Module - UI (Wrapper)

## Description

Save & Export Module - UI (Wrapper)

## Usage

```r
mod_save_export_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing all save/export UI components

# `mod_session_server`: Session Module - Server

## Description

Session Module - Server

## Usage

```r
mod_session_server(
  id,
  status_msg,
  rv,
  load_data,
  refresh_nmr_plot,
  parent_session,
  parent_input
)
```

## Arguments

* `id`: Character. The module's namespace ID
* `status_msg`: ReactiveVal. Shared status message reactive value
* `rv`: List. Shared reactive values
* `load_data`: List. Return value from mod_load_data_server
* `refresh_nmr_plot`: Function. Function to refresh the NMR plot
* `parent_session`: Shiny session. Parent session for updating inputs
* `parent_input`: Shiny input. Parent input object for reading UI values

## Value

NULL (side effects only)

# `mod_session_ui`: Session Module - UI

## Description

Session Module - UI

## Usage

```r
mod_session_ui(id)
```

## Arguments

* `id`: Character. The module's namespace ID

## Value

A tagList containing session save/load controls

# `parse_keep_peak_ranges`: Parse keep_peak_ranges from text input

## Description

Parses a semicolon-separated string of coordinate pairs defining
exclusion zones for peak detection (e.g., solvent regions).

## Usage

```r
parse_keep_peak_ranges(input_string)

parse_keep_peak_ranges(input_string)
```

## Arguments

* `input_string`: Character. Format: "max1,min1; max2,min2; ..."
* `text`: Character string like "0.5,-0.5; 1,0.8; 1.55,1.45;"

## Value

List of numeric vectors, each with 2 elements (min, max)
List of numeric vectors, each with 2 elements (max, min), or NULL

## Examples

```r
parse_keep_peak_ranges("0.5,-0.5; 1,0.8")
# Returns: list(c(0.5, -0.5), c(1, 0.8))
```

# `read_bruker_cached`: Read Bruker spectrum with caching

## Description

Reads a Bruker NMR spectrum from disk, caching the result for subsequent
calls with the same path. This significantly improves performance when
the same spectrum is accessed multiple times.

## Usage

```r
read_bruker_cached(path, dim = "2D")
```

## Arguments

* `path`: Character. Path to the Bruker spectrum directory
* `dim`: Character. Dimension of the spectrum, either "1D" or "2D" (default: "2D")

## Note

Requires the read_bruker() function from Read_2DNMR_spectrum.R

## Value

List containing the spectrum data (from read_bruker function)

# `recenter_box_f2_shift`: Calculate box intensities across multiple spectra

## Description

Calculates peak intensities for a set of reference bounding boxes
across multiple spectra. Handles data validation, duplicate detection,
and optional spectrum alignment.

## Usage

```r
recenter_box_f2_shift(
  mat,
  ppm_x,
  ppm_y,
  xmin,
  xmax,
  ymin,
  ymax,
  tol,
  terr_lo = -Inf,
  terr_hi = Inf
)
```

## Arguments

* `terr_lo, terr_hi`: Exclusive F2 territory bounds (ppm). Default +/-Inf
(no neighbour constraint).
* `reference_boxes`: Data frame. Reference boxes with columns:
xmin, xmax, ymin, ymax, and optionally stain_id
* `spectra_list`: Named list of spectrum data objects
* `apply_shift`: Logical. If TRUE, attempts to align spectra (default: FALSE)
* `method`: Character. Integration method: "sum" (default) or "fit"
* `model`: Character. Fitting model: "gaussian" (default) or "voigt"
* `progress`: Function. Progress callback with signature (value, detail)
* `shift_tolerance_ppm`: Numeric. Tolerance for dynamic box recentering (default: 0).
When > 0, each box is temporarily expanded by this amount in all directions,
the local maximum is found within the expanded region, and the box is
recentered on that maximum ONLY IF the new position captures more intensity
than the original position. This compensates for small chemical shift
variations between spectra (e.g., due to pH, temperature) while avoiding
false recentering on noise or artifacts.
Typical values: 0.01-0.05 ppm. Set to 0 to disable.

## Note

Requires get_box_intensity() and optionally calculate_fitted_volumes()
Compute the F2-only recentering shift for a box (Option B: bounded shift)
The search window IS widened by `tol` on F2 (so a peak that has drifted
slightly out of the reference box can still be recovered). Within that
widened window the INTENSITY-WEIGHTED F2 centroid is computed, matching the
peak-picking convention (weighted.mean(F2_ppm, Volume)) used elsewhere in
APPIN, rather than a raw pixel maximum (more stable, sub-pixel, less
sensitive to single hot pixels). The resulting displacement of the box
centre is then CLAMPED to +/- tol on F2. This keeps chemical-shift
compensation while preventing a box from snapping onto a far-away dominant
feature such as the diagonal / water ridge in homonuclear experiments
(COSY / TOCSY / UF-COSY), since such a feature lies further than `tol` away.
F1 (y) is never shifted, consistent with the F2-only UI.
Both the search window and the final box centre are additionally confined to
the exclusive F2 territory [terr_lo, terr_hi](terr_lo,%20terr_hi) (midpoints to F1-overlapping
neighbours). This stops neighbouring boxes from competing for the same peak
or overlapping after recentering, so `tol` can be larger without conflicts.

## Value

Data frame with columns:

* stain_id: Box identifier
* F2_ppm, F1_ppm: Box center coordinates
* xmin, xmax, ymin, ymax: Box boundaries
* Intensity_<spectrum_name>: One column per spectrum with intensities

Numeric F2 shift (ppm) to add to xmin/xmax. 0 if no valid signal.

