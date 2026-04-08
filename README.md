# sc-skill-creator

![Python](https://img.shields.io/badge/language-Python-3776AB?logo=python)
![R](https://img.shields.io/badge/language-R-276DC3?logo=r)
![License](https://img.shields.io/badge/license-MIT-green)

A **skill generator** that converts single-cell analysis function documentation into Agent-executable skills with built-in biological reasoning and data-driven feedback loops.

---

## What is this?

`sc-skill-creator` is a [Claude Code](https://claude.ai/code) skill that acts as a **protocol converter** — it takes raw function documentation (docstrings, `?help` output) from any single-cell analysis tool and outputs a structured, Agent-ready skill folder.

It is **not** a simple code generator. Every generated skill has three integrated layers:

| Layer | File | Purpose |
|-------|------|---------|
| **Execution** | `execution.py` / `execution.R` | How to run the tool |
| **Cognitive** | `skill.yaml` | Why it works biologically |
| **Critic** | `critic.py` / `critic.R` | Data-driven quality evaluation |

The **Critic layer** is the key innovation: instead of asking the Agent to guess whether results are good, it extracts hard metrics and returns structured `feedback` strings. The Agent looks these up directly in `skill.yaml`'s `parameter_science_guide` to decide how to adjust parameters — no fuzzy reasoning required.

---

## Supported Tools

`sc-skill-creator` is not limited to any specific package. It can generate skills for **any** single-cell analysis function in Python or R, including (but not limited to):

- **Python**: scanpy, scvi-tools, bbknn, harmonypy, cellrank, squidpy, diffxpy, pySCENIC, …
- **R**: Seurat, Harmony, DESeq2, edgeR, MAST, Monocle3, scater, scran, Signac, miloR, …

---

## Generated Skill Structure

```
<package>_<function>/
├── skill.yaml      # Metadata, cognitive layer, parameter science guide
├── execution.py    # Run function (or execution.R for R skills)
├── critic.py       # Quality evaluator returning feedback list (or critic.R)
└── README.md       # Auto-generated — do not edit manually
```

### Example: `scanpy_leiden`

```
scanpy_leiden/
├── skill.yaml
├── execution.py
├── critic.py
└── README.md
```

---

## Four Production Robustness Patterns

Every generated skill enforces four patterns:

### 1. Memory-First
Accepts both in-memory objects and file paths — avoids repeated I/O in multi-round Agent iterations.
```python
adata = input_data if isinstance(input_data, ad.AnnData) else sc.read(input_data)
```

### 2. Parameter Safety
Merges caller-supplied params over defaults — never crashes when the Agent omits a parameter.
```python
current_params = {**default_params, **(params_dict or {})}
```

### 3. Context-Aware Thresholds
Evaluation thresholds adapt to the sequencing protocol (10x Genomics, Smart-seq2, Drop-seq, …).
```python
t = thresholds.get(protocol, thresholds['10x Genomics'])
```

### 4. Metadata Footprinting
Every execution is recorded in `adata.uns['analysis_history']` for full traceability.
```python
adata.uns['analysis_history'].append(_to_serializable({
    'skill_id': 'scanpy_leiden',
    'params': current_params,
    'metrics': {'n_clusters': n_clusters},
    'timestamp': datetime.now().isoformat()
}))
```
> `_to_serializable()` is included in every `execution.py` to convert numpy types to Python natives, preventing `TypeError` when saving `.h5ad` files.

---

## Agent Feedback Loop

```
Agent calls run_skill(adata, params)
        ↓
execution.py runs the analysis
        ↓
Agent calls critic_post_process(adata, context)
        ↓
critic returns feedback = ['too_few_clusters', ...]
        ↓
Agent looks up skill.yaml:
    parameter_science_guide['resolution']['too_few_clusters']
    → adjust: increase, delta: 0.3
        ↓
Agent adjusts params and re-runs → until feedback is empty
```

---

## Cross-Language Data Exchange

Skills declare their accepted input types in `skill.yaml`:

```yaml
accepted_input_types:
  - AnnData       # Python in-memory
  - h5ad          # file path — universal interchange format
output_format: same
```

`.h5ad` serves as the **universal interchange format** between Python and R skills. Each R skill's `execution.R` handles conversion internally (e.g., `zellkonverter::readH5AD()` for Seurat/SCE, count matrix extraction for DESeq2/edgeR), so no separate conversion step is needed.

---

## Tooling

### Validate a generated skill
```bash
python sc-skill-creator/scripts/validate_skill.py path/to/scanpy_leiden/
```

### Generate / refresh README.md
```bash
# Single skill
python sc-skill-creator/scripts/generate_readme.py path/to/scanpy_leiden/

# All skills under a directory
python sc-skill-creator/scripts/generate_readme.py path/to/skills/ --all
```

### Load a skill in Python
```python
from sc-skill-creator.scripts.skill_registry import SkillRegistry

registry = SkillRegistry('path/to/skills/')
registry.scan()
skill = registry.get_skill('scanpy_leiden')
```

---

## References

- `sc-skill-creator/references/scanpy_reference.md` — biological context for Scanpy functions
- `sc-skill-creator/references/r_tools_reference.md` — biological context for R packages (Seurat, Harmony, DESeq2, Monocle3, …)

---

## Repository Layout

```
sc-skill-creator/
├── SKILL.md                        # Claude Code skill definition (generation spec)
├── references/
│   ├── scanpy_reference.md
│   └── r_tools_reference.md
└── scripts/
    ├── skill_registry.py           # Scan and load skills at runtime
    ├── load_skill_folder.py        # Load skill.yaml + execution + critic into dict
    ├── validate_skill.py           # Check robustness patterns
    ├── generate_readme.py          # Auto-generate README.md from skill.yaml
    ├── parse_skill_md.py           # Parse legacy skill.md format
    ├── migrate_json_to_md.py       # Migrate legacy skill.json → skill.md
    └── compile_md_to_json.py       # Compile skill.md → JSON

scanpy_leiden/                      # Example skill: Leiden clustering
├── skill.yaml
├── execution.py
├── critic.py
└── README.md
```

---

## License

MIT
