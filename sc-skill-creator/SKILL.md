---
name: sc-skill-creator
description: Generate separate-script skills for any single-cell analysis tool — Python or R. Use when user wants to create a skill for any single-cell function (Scanpy, Seurat, DESeq2, Monocle3, or any other tool). Output is a folder with skill.yaml (metadata + cognitive layer + parameter_science_guide), execution.py or execution.R (run function), critic.py or critic.R (evaluation), and README.md (auto-generated). This is a protocol converter that transforms raw function docs into Agent-executable skills with data-driven feedback loops.
---

# sc-skill-creator

A skill generator that creates scRNA-seq tool skills from **any** single-cell analysis function — Python or R.

## Scope

This tool is not limited to any specific package. It can generate skills for **any** single-cell analysis function in any language, including (but not limited to):

**Python examples**: scanpy, scvi-tools, bbknn, diffxpy, harmonypy, pySCENIC, cellrank, squidpy, decoupler, liana, …

**R examples**: Seurat, harmony, DESeq2, edgeR, MAST, Monocle3, scater, scran, Signac, muscat, miloR, …

The only requirement is that the user provides the function's documentation (docstring / `?help` output) so the skill can be generated accurately.

## Core Philosophy

This is NOT a simple code generator. It is a **protocol converter** that transforms raw function documentation into Agent-executable skills with three integrated layers:

1. **Execution Layer** — How to run the tool
2. **Cognitive Layer** — Why the tool works biologically
3. **Critic Layer** — Data-driven feedback on results

The Critic layer is the key innovation: it generates post-processing code that **extracts hard metrics** and **generates semantic feedback** rather than relying on the Agent to guess whether results are good.

---

## Four Production Robustness Requirements

Every generated Skill MUST incorporate these four patterns in BOTH Python and R variants.

### 1. Memory-First Execution

Accept both an in-memory object and a file path to avoid repeated I/O in multi-round iterations.

**Python (AnnData)**
```python
def run_skill(input_data, params_dict=None):
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)
    # ...
    if not isinstance(input_data, ad.AnnData):
        adata.write(output_path)
    return adata
```

**R (Seurat / SingleCellExperiment)**
```r
run_skill <- function(input_data, params_list = NULL) {
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }
    # ...
    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
```

### 2. Parameter Safety with Dynamic Defaults

Never use bare `{{placeholder}}`. Merge caller params over defaults so missing params fall back safely.

**Python**
```python
def run_skill(input_data, params_dict=None, default_params=None):
    default_params = default_params or {'resolution': 1.0, 'n_neighbors': 15}
    current_params = {**default_params, **(params_dict or {})}
    sc.tl.leiden(adata, resolution=current_params['resolution'])
```

**R**
```r
run_skill <- function(input_data, params_list = NULL) {
    default_params <- list(resolution = 1.0, n.neighbors = 15L)
    current_params <- modifyList(default_params, if (!is.null(params_list)) params_list else list())
    Seurat::FindClusters(obj, resolution = current_params$resolution)
```

### 3. Context-Aware Thresholds

Parameterize protocol-specific thresholds — do not hardcode.

**Python**
```python
def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')
    thresholds = {
        '10x Genomics': {'min_clusters': 3, 'max_clusters': 200, 'min_cluster_size': 5},
        'Smart-seq2':   {'min_clusters': 5, 'max_clusters': 100, 'min_cluster_size': 10},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])
```

**R**
```r
critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"
    thresholds <- list(
        "10x Genomics" = list(min_clusters = 3L, max_clusters = 200L, min_cluster_size = 5L),
        "Smart-seq2"   = list(min_clusters = 5L, max_clusters = 100L, min_cluster_size = 10L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]
```

### 4. Metadata Footprinting

Record execution to the object's unstructured metadata for downstream traceability.

**Python (adata.uns)**

⚠️ **h5ad serialization requirement**: `adata.uns` is written to HDF5 when saving. h5py does **not** accept numpy scalar types (`np.int64`, `np.float32`, `np.ndarray`, etc.) — only Python native types (`int`, `float`, `str`, `list`, `dict`, `bool`). Always pass the history entry through `_to_serializable()` before appending, otherwise `adata.write_h5ad()` will raise a `TypeError`.

```python
import numpy as np

def _to_serializable(obj):
    """Recursively convert numpy types to Python natives for h5ad compatibility."""
    if isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(i) for i in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj

# In run_skill():
if 'analysis_history' not in adata.uns:
    adata.uns['analysis_history'] = []
adata.uns['analysis_history'].append(_to_serializable({
    'skill_id': 'scanpy_leiden',
    'params': current_params,
    'metrics': {'n_clusters': n_clusters},
    'timestamp': datetime.now().isoformat()
}))
```

`_to_serializable()` MUST be included in every Python `execution.py`. It is not imported from a library — paste it directly into the file.

**R (Seurat @misc / SCE metadata)**
```r
# Seurat
if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
    skill_id  = "seurat_findclusters",
    params    = current_params,
    metrics   = list(n_clusters = n_clusters),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)))

# SingleCellExperiment
if (is.null(S4Vectors::metadata(sce)$analysis_history))
    S4Vectors::metadata(sce)$analysis_history <- list()
S4Vectors::metadata(sce)$analysis_history <- c(
    S4Vectors::metadata(sce)$analysis_history,
    list(list(skill_id = "scran_cluster", params = current_params, ...))
)
```

---

## Input Specification

User provides:
- **Tool + function name**: e.g., `scanpy.tl.leiden`, `Seurat::FindClusters`, `DESeq2::DESeq`
- **Function docstring / help output**: The official documentation
- **Language**: `python` or `R`
- **Optional context**: Parameter tuning goals, typical dataset properties

---

## Output Schema

Every generated skill is a **folder** with four files:

```
<package>_<function_name>/        ← skill_id is also the folder name
├── skill.yaml      ← structured data (metadata + cognitive layer + parameter_science_guide)
├── execution.py    ← Python run function   (or execution.R for R skills)
├── critic.py       ← Python critic         (or critic.R for R skills)
└── README.md       ← auto-generated, DO NOT edit manually
```

### skill.yaml structure

```yaml
# One-line description of what this skill does
skill_id: <package>_<function_name>    # e.g. scanpy_leiden, seurat_findclusters

language: python                        # python | R
package: scanpy                         # the package name, e.g. scanpy, seurat, deseq2, monocle3, ...

required_inputs:
  - input_data
  - params_dict                         # params_list for R skills

default_params:
  param_name: default_value             # comment explaining the param

output_objects:
  - updated_adata                       # or seurat_obj / sce_obj

# Accepted input types — Agent uses this to decide whether a conversion step is needed
# before calling this skill. List ALL types the execution script can handle.
accepted_input_types:
  - AnnData                             # Python in-memory object
  - h5ad                                # file path to .h5ad
  # R examples: Seurat, SingleCellExperiment, DESeqDataSet, cell_data_set, h5ad, rds

# output_format — what this skill writes back when input was a file path
# 'same' means write the same format that was read (default)
output_format: same                     # same | h5ad | rds

metrics_to_extract:
  - metric_name

context_params:
  - protocol

success_thresholds: "n_clusters >= t['min_clusters'] and n_clusters <= t['max_clusters']"

error_handling:                         # optional
  ErrorType: "diagnostic suggestion"

cognitive_layer:
  purpose: >
    Multi-line biological purpose of this tool.
  parameter_impact:
    param_name: "Causal relationship when you increase/decrease this parameter."

parameter_science_guide:
  param_name:
    condition_name:              # MUST match feedback keys returned by critic file
      adjust: increase           # 'increase' or 'decrease'
      delta: 0.3                 # numeric amount to change by
      causal_chain: "mechanism"
      expected_effect: "outcome"
```

### execution.py / execution.R structure

**Python**
```python
"""
<package>_<function_name> — One-line title

PURPOSE:
    Brief biological explanation.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    param=value, ...

OUTPUTS:
    AnnData with adata.obs['key'] and adata.uns['analysis_history']
"""
import scanpy as sc
import anndata as ad
import numpy as np
from datetime import datetime


def _to_serializable(obj):
    """Recursively convert numpy types to Python natives for h5ad compatibility."""
    if isinstance(obj, dict):          return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)): return [_to_serializable(i) for i in obj]
    if isinstance(obj, np.integer):    return int(obj)
    if isinstance(obj, np.floating):   return float(obj)
    if isinstance(obj, np.ndarray):    return obj.tolist()
    if isinstance(obj, np.bool_):      return bool(obj)
    return obj


def run_<function_name>(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or { ... }
    current_params = {**default_params, **(params_dict or {})}

    # ... execute sc.<module>.<function> ...

    # Metadata Footprinting — _to_serializable prevents h5ad TypeError on numpy types
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': '<package>_<function_name>',
        'params': current_params,
        'metrics': { ... },
        'timestamp': datetime.now().isoformat()
    }))
    return adata
```

**R**
```r
#' <package>_<function_name> — One-line title
#'
#' PURPOSE:
#'   Brief biological explanation.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat | SingleCellExperiment
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   param = value, ...
#'
#' OUTPUTS:
#'   Seurat/SCE object with analysis_history in @misc / metadata()

library(<package>)

run_<function_name> <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list( ... )
    current_params <- modifyList(default_params, if (!is.null(params_list)) params_list else list())

    # ... execute <Package>::<Function>() ...

    # Metadata Footprinting
    if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
    obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
        skill_id  = "<package>_<function_name>",
        params    = current_params,
        metrics   = list( ... ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
```

### critic.py / critic.R structure

**Python**
```python
"""
<package>_<function_name> critic — Evaluates result quality

FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
    condition_name → description of what triggers this feedback
"""

def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')

    thresholds = {
        '10x Genomics': { ... },
        'Smart-seq2':   { ... },
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []   # ← list of condition strings matching parameter_science_guide keys

    # ... extract metrics, check thresholds ...

    if <condition>:
        warnings.append("human-readable message")
        feedback.append('condition_name')   # ← exact key from skill.yaml guide

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,   # Agent uses this to look up parameter_science_guide
        'success': not bool(feedback),
        'context_used': protocol
    }
```

**R**
```r
#' critic — Evaluates result quality
#'
#' FEEDBACK KEYS (must match parameter_science_guide conditions in skill.yaml):
#'   condition_name → description of what triggers this feedback

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"

    thresholds <- list(
        "10x Genomics" = list( ... ),
        "Smart-seq2"   = list( ... )
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)   # ← strings matching parameter_science_guide keys

    # ... extract metrics, check thresholds ...

    if (<condition>) {
        warnings <- c(warnings, "human-readable message")
        feedback <- c(feedback, "condition_name")
    }

    list(
        metrics      = metrics,
        warnings     = warnings,
        feedback     = feedback,
        success      = length(feedback) == 0L,
        context_used = protocol
    )
}
```

**Critical contract**: `feedback` values returned by the critic file must **exactly match**
the condition keys in `skill.yaml`'s `parameter_science_guide`. The Agent does a direct
dict lookup: `guide[param][feedback_item]` → no fuzzy matching.

---

## Step-by-Step Process

### Step 1: Identify Tool and Language

Determine:
- **Language**: `python` or `R`
- **Package**: e.g., `scanpy`, `seurat`, `harmony`, `deseq2`
- **Object type**: `AnnData` (Python) | `Seurat` | `SingleCellExperiment` | `DESeqDataSet`
- **skill_id**: `<package>_<function_name>` — all lowercase, underscore-separated
  - `scanpy_leiden`, `seurat_findclusters`, `harmony_runharmony`, `deseq2_deseq`

### Step 2: Parse Function Signature

Extract from the docstring / `?FunctionName` help:
- All parameters with types and default values
- Return values and side effects (what gets added to the object)
- Prerequisite steps (e.g., must run `RunPCA` before `RunUMAP` in Seurat)

### Step 3: Generate Execution Layer

Create an execution file with:
- **Memory-First**: accept both path and in-memory object
- **Parameter Safety**: `{**default_params, **agent_params}` (Python) or `modifyList()` (R)
- **Metadata Footprinting**: append to `adata.uns['analysis_history']` or `obj@misc$analysis_history`
- Use `current_params['key']` (Python) or `current_params$key` (R) — never bare placeholders

### Step 4: Generate Critic Layer (Data-Driven)

Generate post-processing code that:
1. **Extracts hard metrics** from the result object automatically
2. **Applies context-aware thresholds** (protocol, species, expected cell count)
3. **Returns structured feedback** — a list of condition strings tied to `parameter_science_guide`

### Step 5: Generate Cognitive Layer + Parameter Science Guide

For each tunable parameter:
- **purpose**: what it controls biologically
- **parameter_impact**: causal chain (parameter change → intermediate effect → biological outcome)
- **parameter_science_guide**: `{condition: {adjust, delta, causal_chain, expected_effect}}`
  - Condition names MUST match feedback strings returned by critic

---

## Example A: Leiden Clustering (Python / Scanpy)

**`skill.yaml`**
```yaml
# Leiden clustering — community detection on kNN graph
skill_id: scanpy_leiden
language: python
package: scanpy

required_inputs:
  - input_data
  - params_dict

default_params:
  resolution: 1.0    # float  ↑→more clusters  ↓→fewer clusters
  n_neighbors: 15    # int    ↑→smoother graph  ↓→finer structure
  random_state: 42
  key_added: leiden

output_objects:
  - updated_adata

accepted_input_types:
  - AnnData
  - h5ad

output_format: same

metrics_to_extract:
  - n_clusters
  - min_cluster_size
  - max_cluster_size
  - mean_cluster_size

context_params:
  - protocol

success_thresholds: "n_clusters >= t['min_clusters'] and n_clusters <= t['max_clusters'] and min_cluster_size >= t['min_cluster_size']"

error_handling:
  MemoryError: "Reduce n_neighbors or n_pcs"
  "KeyError: neighbors": "Run sc.pp.neighbors() before sc.tl.leiden()"

cognitive_layer:
  purpose: >
    Leiden clustering partitions single cells into groups based on transcriptional
    similarity using community detection on a k-nearest neighbor graph.
    Identifies cell populations (cell types, states, subtypes) without prior labels.
  parameter_impact:
    resolution: "Controls cluster granularity. Higher values produce more, smaller clusters."
    n_neighbors: "Determines k in kNN graph. More neighbors smooth transitions between populations."

parameter_science_guide:
  resolution:
    too_few_clusters:
      adjust: increase
      delta: 0.3
      causal_chain: "Higher resolution lowers modularity threshold, allowing more communities"
      expected_effect: "More clusters, smaller sizes"
    too_many_clusters:
      adjust: decrease
      delta: 0.3
      causal_chain: "Lower resolution raises modularity threshold, causing similar cells to merge"
      expected_effect: "Fewer clusters, larger sizes"
    small_clusters:
      adjust: decrease
      delta: 0.2
      causal_chain: "Lower resolution merges small clusters with their nearest neighbors"
      expected_effect: "Small clusters absorbed into adjacent populations"
  n_neighbors:
    noisy_fragmented:
      adjust: increase
      delta: 7
      causal_chain: "More neighbors smooth the kNN graph by averaging over more cells"
      expected_effect: "Smoother clusters, fewer spurious small groups"
```

**`execution.py`**
```python
"""
scanpy_leiden — Leiden Clustering

PURPOSE:
    Partitions single cells into groups based on transcriptional similarity
    using community detection on a k-nearest neighbor graph.

INPUTS:
    input_data  : str | AnnData
    params_dict : dict (optional)

DEFAULTS:
    resolution=1.0, n_neighbors=15, random_state=42, key_added='leiden'

OUTPUTS:
    AnnData with adata.obs['leiden'] and adata.uns['analysis_history']
"""
import scanpy as sc
import anndata as ad
import numpy as np
from datetime import datetime


def _to_serializable(obj):
    """Recursively convert numpy types to Python natives for h5ad compatibility.
    adata.uns is saved to HDF5; h5py rejects numpy scalars and arrays."""
    if isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(i) for i in obj]
    if isinstance(obj, np.integer):   return int(obj)
    if isinstance(obj, np.floating):  return float(obj)
    if isinstance(obj, np.ndarray):   return obj.tolist()
    if isinstance(obj, np.bool_):     return bool(obj)
    return obj


def run_leiden(input_data, params_dict=None, default_params=None):
    # Memory-First
    adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)

    # Parameter Safety
    default_params = default_params or {
        'resolution': 1.0, 'n_neighbors': 15, 'random_state': 42, 'key_added': 'leiden'
    }
    current_params = {**default_params, **(params_dict or {})}

    sc.tl.leiden(
        adata,
        resolution=current_params['resolution'],
        n_neighbors=current_params['n_neighbors'],
        random_state=current_params['random_state'],
        key_added=current_params['key_added']
    )

    key = current_params['key_added']
    cluster_sizes = adata.obs[key].value_counts()

    # Metadata Footprinting — _to_serializable prevents h5ad write errors
    if 'analysis_history' not in adata.uns:
        adata.uns['analysis_history'] = []
    adata.uns['analysis_history'].append(_to_serializable({
        'skill_id': 'scanpy_leiden',
        'params': current_params,
        'metrics': {
            'n_clusters': adata.obs[key].nunique(),
            'n_cells': adata.n_obs,
            'min_cluster_size': int(cluster_sizes.min()) if len(cluster_sizes) else 0
        },
        'timestamp': datetime.now().isoformat()
    }))
    return adata
```

**`critic.py`**
```python
"""
scanpy_leiden critic — Evaluates clustering quality

FEEDBACK KEYS:
    too_few_clusters  → resolution too low
    too_many_clusters → resolution too high
    small_clusters    → resolution slightly too high or noise present
    noisy_fragmented  → n_neighbors too low
"""

def critic_post_process(adata, context=None):
    context = context or {}
    protocol = context.get('protocol', '10x Genomics')

    thresholds = {
        '10x Genomics': {'min_clusters': 3, 'max_clusters': 200, 'min_cluster_size': 5},
        'Smart-seq2':   {'min_clusters': 5, 'max_clusters': 100, 'min_cluster_size': 10},
        'Drop-seq':     {'min_clusters': 3, 'max_clusters': 150, 'min_cluster_size': 5},
    }
    t = thresholds.get(protocol, thresholds['10x Genomics'])

    metrics = {}
    warnings = []
    feedback = []

    if 'leiden' not in adata.obs:
        return {'metrics': {}, 'warnings': ['CRITICAL: leiden not in adata.obs'],
                'feedback': [], 'success': False, 'context_used': protocol}

    counts = adata.obs['leiden'].value_counts()
    metrics['n_clusters'] = adata.obs['leiden'].nunique()
    metrics['min_cluster_size'] = int(counts.min())
    metrics['max_cluster_size'] = int(counts.max())
    metrics['mean_cluster_size'] = float(counts.mean())

    if metrics['n_clusters'] < t['min_clusters']:
        warnings.append(f"UNDER-SEGMENTATION: {metrics['n_clusters']} clusters for {protocol}")
        feedback.append('too_few_clusters')

    if metrics['n_clusters'] > t['max_clusters']:
        warnings.append(f"OVER-SEGMENTATION: {metrics['n_clusters']} clusters for {protocol}")
        feedback.append('too_many_clusters')

    if metrics['min_cluster_size'] < t['min_cluster_size']:
        warnings.append(f"SMALL CLUSTERS: min={metrics['min_cluster_size']} cells")
        feedback.append('small_clusters')

    return {
        'metrics': metrics,
        'warnings': warnings,
        'feedback': feedback,
        'success': not bool(feedback),
        'context_used': protocol
    }
```

---

## Example B: FindClusters (R / Seurat)

**`skill.yaml`**
```yaml
# Seurat FindClusters — Louvain/Leiden community detection
skill_id: seurat_findclusters
language: R
package: seurat

required_inputs:
  - input_data
  - params_list

default_params:
  resolution: 0.8       # float  ↑→more clusters  ↓→fewer clusters
  algorithm: 1          # 1=Louvain, 2=LME, 3=SLM, 4=Leiden
  random.seed: 42

output_objects:
  - seurat_obj

metrics_to_extract:
  - n_clusters
  - min_cluster_size
  - max_cluster_size

context_params:
  - protocol

success_thresholds: "n_clusters >= t['min_clusters'] and n_clusters <= t['max_clusters']"

error_handling:
  "Error in FindClusters: graph not found": "Run FindNeighbors() before FindClusters()"

cognitive_layer:
  purpose: >
    FindClusters identifies cell populations by applying community detection on
    the shared nearest neighbor (SNN) graph built by FindNeighbors. Resolution
    controls granularity: higher values yield more, smaller clusters.
  parameter_impact:
    resolution: "Primary granularity knob. Increase to split large clusters; decrease to merge over-fragmented ones."
    algorithm: "Community detection algorithm. Leiden (4) is more accurate but slower; Louvain (1) is standard."

parameter_science_guide:
  resolution:
    too_few_clusters:
      adjust: increase
      delta: 0.2
      causal_chain: "Higher resolution lowers the modularity bar, revealing finer community structure"
      expected_effect: "More clusters, smaller median size"
    too_many_clusters:
      adjust: decrease
      delta: 0.2
      causal_chain: "Lower resolution raises the modularity bar, merging similar communities"
      expected_effect: "Fewer clusters, larger median size"
    small_clusters:
      adjust: decrease
      delta: 0.1
      causal_chain: "Marginal communities disappear when the bar rises slightly"
      expected_effect: "Small spurious clusters absorbed into neighbors"
```

**`execution.R`**
```r
#' seurat_findclusters — Seurat FindClusters
#'
#' PURPOSE:
#'   Identifies cell populations via community detection on the SNN graph.
#'
#' INPUTS:
#'   input_data  : character (RDS path) | Seurat object
#'   params_list : list (optional, overrides defaults)
#'
#' DEFAULTS:
#'   resolution = 0.8, algorithm = 1L, random.seed = 42L
#'
#' OUTPUTS:
#'   Seurat object with Idents() set and @misc$analysis_history updated

library(Seurat)

run_findclusters <- function(input_data, params_list = NULL) {
    # Memory-First
    if (is.character(input_data)) {
        obj <- readRDS(input_data)
    } else {
        obj <- input_data
    }

    # Parameter Safety
    default_params <- list(resolution = 0.8, algorithm = 1L, random.seed = 42L)
    current_params <- modifyList(default_params, if (!is.null(params_list)) params_list else list())

    obj <- Seurat::FindClusters(
        obj,
        resolution  = current_params$resolution,
        algorithm   = current_params$algorithm,
        random.seed = current_params$random.seed
    )

    cluster_sizes <- table(Idents(obj))
    n_clusters    <- length(cluster_sizes)

    # Metadata Footprinting
    if (is.null(obj@misc$analysis_history)) obj@misc$analysis_history <- list()
    obj@misc$analysis_history <- c(obj@misc$analysis_history, list(list(
        skill_id  = "seurat_findclusters",
        params    = current_params,
        metrics   = list(
            n_clusters       = n_clusters,
            n_cells          = ncol(obj),
            min_cluster_size = as.integer(min(cluster_sizes))
        ),
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    )))

    if (is.character(input_data)) saveRDS(obj, input_data)
    return(obj)
}
```

**`critic.R`**
```r
#' seurat_findclusters critic — Evaluates clustering quality
#'
#' FEEDBACK KEYS (match parameter_science_guide in skill.yaml):
#'   too_few_clusters  → resolution too low
#'   too_many_clusters → resolution too high
#'   small_clusters    → marginal communities present

critic_post_process <- function(obj, context = NULL) {
    if (is.null(context)) context <- list()
    protocol <- if (!is.null(context$protocol)) context$protocol else "10x Genomics"

    thresholds <- list(
        "10x Genomics" = list(min_clusters = 3L, max_clusters = 200L, min_cluster_size = 5L),
        "Smart-seq2"   = list(min_clusters = 5L, max_clusters = 100L, min_cluster_size = 10L),
        "Drop-seq"     = list(min_clusters = 3L, max_clusters = 150L, min_cluster_size = 5L)
    )
    t <- if (!is.null(thresholds[[protocol]])) thresholds[[protocol]] else thresholds[["10x Genomics"]]

    metrics  <- list()
    warnings <- character(0)
    feedback <- character(0)

    cluster_sizes <- table(Idents(obj))
    if (length(cluster_sizes) == 0L) {
        return(list(metrics = list(), warnings = "CRITICAL: no clusters found",
                    feedback = character(0), success = FALSE, context_used = protocol))
    }

    metrics$n_clusters       <- length(cluster_sizes)
    metrics$min_cluster_size <- as.integer(min(cluster_sizes))
    metrics$max_cluster_size <- as.integer(max(cluster_sizes))
    metrics$mean_cluster_size <- mean(cluster_sizes)

    if (metrics$n_clusters < t$min_clusters) {
        warnings <- c(warnings, sprintf("UNDER-SEGMENTATION: %d clusters for %s",
                                        metrics$n_clusters, protocol))
        feedback <- c(feedback, "too_few_clusters")
    }
    if (metrics$n_clusters > t$max_clusters) {
        warnings <- c(warnings, sprintf("OVER-SEGMENTATION: %d clusters for %s",
                                        metrics$n_clusters, protocol))
        feedback <- c(feedback, "too_many_clusters")
    }
    if (metrics$min_cluster_size < t$min_cluster_size) {
        warnings <- c(warnings, sprintf("SMALL CLUSTERS: min=%d cells", metrics$min_cluster_size))
        feedback <- c(feedback, "small_clusters")
    }

    list(
        metrics      = metrics,
        warnings     = warnings,
        feedback     = feedback,
        success      = length(feedback) == 0L,
        context_used = protocol
    )
}
```

---

## Validation Checklist

Before returning the generated skill, verify:

**skill.yaml**
- [ ] `skill_id` follows `<package>_<function_name>` format (all lowercase)
- [ ] `language` is `python` or `R`
- [ ] `package` matches the tool (scanpy, seurat, harmony, deseq2, …)
- [ ] Contains all required keys: `skill_id`, `language`, `package`, `required_inputs`, `default_params`, `output_objects`, `metrics_to_extract`, `context_params`, `success_thresholds`, `cognitive_layer`, `parameter_science_guide`
- [ ] `parameter_science_guide` condition keys match feedback strings returned by critic

**execution.py / execution.R**
- [ ] Memory-First pattern: `isinstance(input_data, ad.AnnData)` (Python) or `is.character(input_data)` (R)
- [ ] Parameter Safety: `{**default_params, **agent_params}` (Python) or `modifyList(default_params, ...)` (R)
- [ ] Metadata Footprinting: appends to `adata.uns['analysis_history']` or `obj@misc$analysis_history`
- [ ] `_to_serializable()` helper defined in every Python `execution.py` and footprint entry wrapped with it
- [ ] Module docstring (Python `"""..."""`) or roxygen header (R `#' ...`) with PURPOSE, INPUTS, DEFAULTS, OUTPUTS

**critic.py / critic.R**
- [ ] Defines `critic_post_process(adata/obj, context=None/NULL)`
- [ ] Returns `feedback` as a list/vector of strings
- [ ] Every feedback string exactly matches a condition key in `parameter_science_guide`
- [ ] Uses context-aware thresholds (`t['key']` / `t$key`)

Run `python scripts/validate_skill.py path/to/<skill_folder>/` to verify Python skills programmatically.

---

## Usage Patterns

When user says:
- "Create a skill for [any single-cell function]"
- "Generate the execution template for [tool]"
- "Build a critic for [clustering/DE/trajectory method]"
- "I need the cognitive guide for [normalization/integration method]"

Always confirm **language** (Python or R) and **package** before generating.

---

## Tooling

### Validate a skill folder (Python skills)
```bash
python scripts/validate_skill.py path/to/scanpy_leiden/
```

### Generate / refresh README.md
```bash
# Single skill
python scripts/generate_readme.py path/to/scanpy_leiden/

# All skills under a root folder
python scripts/generate_readme.py path/to/skills/ --all
```
README.md is also auto-generated when the registry first scans a new skill.

### Registry
`skill_registry.py` detects format automatically (priority: `skill.yaml` > `skill.md` > `skill.json`).
```python
registry = SkillRegistry('/path/to/skills')
registry.scan()
skill = registry.get_skill('seurat_findclusters')
```

### References
- `references/scanpy_reference.md` — biological context for Scanpy functions
- `references/r_tools_reference.md` — biological context for R packages (Seurat, Harmony, DESeq2, Monocle3, etc.)
