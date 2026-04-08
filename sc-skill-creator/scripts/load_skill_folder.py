#!/usr/bin/env python3
"""
load_skill_folder.py - Load a skill from the separate-scripts format.

Expected folder structure:
    scanpy_leiden/          ← Python skill
    ├── skill.yaml      ← metadata + cognitive_layer + parameter_science_guide
    ├── execution.py    ← run function (Memory-First, Parameter Safety, Footprinting)
    ├── critic.py       ← critic_post_process function (returns feedback list)
    └── README.md       ← auto-generated, not loaded at runtime

    seurat_findclusters/    ← R skill
    ├── skill.yaml
    ├── execution.R
    ├── critic.R
    └── README.md

Returns a dict with the same shape as the JSON skill schema so the registry
and validator work unchanged.
"""

from pathlib import Path
from typing import Any, Dict

try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml")


REQUIRED_YAML_KEYS = [
    'skill_id',
    'required_inputs',
    'default_params',
    'output_objects',
    'metrics_to_extract',
    'context_params',
    'success_thresholds',
    'cognitive_layer',
    'parameter_science_guide',
]


def _find_exec_critic(folder: Path, language: str):
    """Return (exec_path, critic_path) for the given language.

    Python skills: execution.py + critic.py
    R skills:      execution.R  + critic.R
    Auto-detect when language is not specified.
    """
    lang = (language or '').lower()

    if lang == 'r':
        return folder / 'execution.R', folder / 'critic.R'
    elif lang == 'python':
        return folder / 'execution.py', folder / 'critic.py'
    else:
        # Auto-detect: prefer .py, fall back to .R
        exec_py = folder / 'execution.py'
        if exec_py.exists():
            return exec_py, folder / 'critic.py'
        return folder / 'execution.R', folder / 'critic.R'


def load_skill_folder(folder_path: str) -> Dict[str, Any]:
    """Load skill.yaml + execution.py/.R + critic.py/.R into a unified skill dict.

    Parameters
    ----------
    folder_path : str
        Path to the skill folder containing skill.yaml and execution/critic scripts.

    Returns
    -------
    dict matching the standard skill schema:
        skill_id, execution_layer, cognitive_layer, critic_layer,
        parameter_science_guide

    Raises
    ------
    FileNotFoundError  : if any required file is missing
    ValueError         : if skill.yaml is missing required keys
    """
    folder = Path(folder_path)

    yaml_path = folder / 'skill.yaml'
    if not yaml_path.exists():
        raise FileNotFoundError(f"Required file missing: {yaml_path}")

    with open(yaml_path, encoding='utf-8') as f:
        meta = yaml.safe_load(f)

    if not isinstance(meta, dict):
        raise ValueError(f"skill.yaml did not parse to a dict: {yaml_path}")

    missing = [k for k in REQUIRED_YAML_KEYS if k not in meta]
    if missing:
        raise ValueError(f"skill.yaml missing required keys: {missing}")

    language = meta.get('language', '')
    exec_path, critic_path = _find_exec_critic(folder, language)

    for p in (exec_path, critic_path):
        if not p.exists():
            raise FileNotFoundError(f"Required file missing: {p}")

    execution_code = exec_path.read_text(encoding='utf-8')
    critic_code = critic_path.read_text(encoding='utf-8')

    cog = meta['cognitive_layer']

    critic_layer: Dict[str, Any] = {
        'metrics_to_extract': meta['metrics_to_extract'],
        'success_thresholds': meta['success_thresholds'],
        'context_params': meta['context_params'],
        'post_processing_code': critic_code,
    }
    if 'error_handling' in meta:
        critic_layer['error_handling'] = meta['error_handling']

    return {
        'skill_id': meta['skill_id'],
        'execution_layer': {
            'code_template': execution_code,
            'required_inputs': meta['required_inputs'],
            'default_params': meta['default_params'],
            'output_objects': meta['output_objects'],
        },
        'cognitive_layer': {
            'purpose': cog.get('purpose', ''),
            'parameter_impact': cog.get('parameter_impact', {}),
        },
        'critic_layer': critic_layer,
        'parameter_science_guide': meta['parameter_science_guide'],
    }


if __name__ == '__main__':
    import sys
    import json

    if len(sys.argv) < 2:
        print("Usage: load_skill_folder.py <skill_folder>")
        sys.exit(1)

    skill = load_skill_folder(sys.argv[1])
    print(json.dumps(skill, indent=2, ensure_ascii=False))
