#!/usr/bin/env python3
"""
generate_readme.py - Auto-generate README.md from a skill folder.

Sources (in priority order):
  1. skill.yaml  → purpose, parameter_impact, parameter_science_guide, defaults
  2. execution.py module docstring → usage notes (if present)
  3. critic.py module docstring   → feedback key descriptions (if present)

Usage:
    python generate_readme.py path/to/skill_folder
    python generate_readme.py path/to/skills_root --all
"""

import argparse
import ast
import sys
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml")


# ──────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────

def _get_module_docstring(py_path: Path) -> Optional[str]:
    """Extract the module-level docstring from a .py file without importing it."""
    try:
        tree = ast.parse(py_path.read_text(encoding='utf-8'))
        return ast.get_docstring(tree)
    except Exception:
        return None


def _truncate(text: str, max_len: int = 80) -> str:
    return text[:max_len] + '...' if len(text) > max_len else text


# ──────────────────────────────────────────
# README builder
# ──────────────────────────────────────────

def _language_badge(language: str, package: str) -> str:
    """Return a markdown badge line for language and package."""
    badges = []
    lang = (language or '').lower()
    if lang == 'r':
        badges.append('![R](https://img.shields.io/badge/language-R-276DC3?logo=r)')
    elif lang == 'python':
        badges.append('![Python](https://img.shields.io/badge/language-Python-3776AB?logo=python)')
    if package:
        pkg_display = package.capitalize()
        badges.append(f'![{pkg_display}](https://img.shields.io/badge/package-{pkg_display}-lightgrey)')
    return '  '.join(badges)


def _get_exec_docstring(folder: Path, language: str) -> Optional[str]:
    """Extract module/file docstring from execution script (Python or R)."""
    lang = (language or '').lower()
    if lang == 'r':
        r_path = folder / 'execution.R'
        if not r_path.exists():
            return None
        # R roxygen header: lines starting with #'
        lines = r_path.read_text(encoding='utf-8').splitlines()
        doc_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#'"):
                doc_lines.append(stripped[2:].lstrip())
            elif doc_lines:
                break  # first non-#' line after header ends it
        return '\n'.join(doc_lines).strip() or None
    else:
        py_path = folder / 'execution.py'
        return _get_module_docstring(py_path) if py_path.exists() else None


def _get_critic_docstring(folder: Path, language: str) -> Optional[str]:
    """Extract module/file docstring from critic script (Python or R)."""
    lang = (language or '').lower()
    if lang == 'r':
        r_path = folder / 'critic.R'
        if not r_path.exists():
            return None
        lines = r_path.read_text(encoding='utf-8').splitlines()
        doc_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#'"):
                doc_lines.append(stripped[2:].lstrip())
            elif doc_lines:
                break
        return '\n'.join(doc_lines).strip() or None
    else:
        py_path = folder / 'critic.py'
        return _get_module_docstring(py_path) if py_path.exists() else None


def build_readme(folder: Path) -> str:
    yaml_path = folder / 'skill.yaml'

    with open(yaml_path, encoding='utf-8') as f:
        meta = yaml.safe_load(f)

    skill_id = meta['skill_id']
    language: str = meta.get('language', '')
    package: str = meta.get('package', '')
    cog = meta.get('cognitive_layer', {})
    purpose = cog.get('purpose', '').strip()
    parameter_impact: Dict[str, str] = cog.get('parameter_impact', {})
    default_params: Dict[str, Any] = meta.get('default_params', {})
    guide: Dict[str, Any] = meta.get('parameter_science_guide', {})
    metrics: list = meta.get('metrics_to_extract', [])
    context_params: list = meta.get('context_params', [])
    error_handling: Dict[str, str] = meta.get('error_handling', {})

    exec_doc   = _get_exec_docstring(folder, language)
    critic_doc = _get_critic_docstring(folder, language)

    lines = []

    # ── Title ──
    lines += [f'# {skill_id}', '']

    # ── Language / package badges ──
    badge_line = _language_badge(language, package)
    if badge_line:
        lines += [badge_line, '']

    lines += [f'> Auto-generated from `skill.yaml`. Edit source files, not this README.', '']

    # ── Purpose ──
    if purpose:
        lines += [purpose, '']

    # ── Quick reference table ──
    lines += ['## Parameters', '']
    lines += ['| Parameter | Default | Effect |', '|-----------|---------|--------|']
    for param, value in default_params.items():
        impact = parameter_impact.get(param, '')
        lines.append(f'| `{param}` | `{value}` | {_truncate(impact)} |')
    lines.append('')

    # ── Parameter impact detail ──
    if parameter_impact:
        lines += ['## Parameter Impact', '']
        for param, desc in parameter_impact.items():
            lines += [f'### `{param}`', '', desc.strip(), '']

    # ── Parameter science guide ──
    if guide:
        lines += ['## Parameter Tuning Guide', '']
        lines += ['The Agent uses this table to adjust parameters when critic feedback is received.', '']
        for param, conditions in guide.items():
            if not isinstance(conditions, dict):
                continue
            lines += [f'### `{param}`', '']
            lines += ['| Condition | Adjust | Causal Chain |',
                      '|-----------|--------|--------------|']
            for condition, entry in conditions.items():
                if not isinstance(entry, dict):
                    continue
                adjust = entry.get('adjust', '')
                delta = entry.get('delta', '')
                adjust_str = f'{adjust} by {delta}' if delta else str(adjust)
                causal = _truncate(entry.get('causal_chain', ''), 70)
                lines.append(f'| `{condition}` | {adjust_str} | {causal} |')
            lines.append('')

    # ── Critic metrics ──
    if metrics:
        lines += ['## Critic Metrics', '']
        for m in metrics:
            lines.append(f'- `{m}`')
        lines.append('')

    # ── Context params ──
    if context_params:
        lines += ['## Context Parameters', '']
        lines += ['Pass via `context` dict to `critic_post_process()`.', '']
        for p in context_params:
            lines.append(f'- `{p}`')
        lines.append('')

    # ── Error handling ──
    if error_handling:
        lines += ['## Error Handling', '']
        lines += ['| Error | Solution |', '|-------|----------|']
        for error, solution in error_handling.items():
            lines.append(f'| `{error}` | {solution} |')
        lines.append('')

    # ── Execution notes from docstring ──
    if exec_doc:
        lines += ['## Execution Notes', '']
        lines += [exec_doc.strip(), '']

    # ── Critic notes from docstring ──
    if critic_doc:
        lines += ['## Critic Notes', '']
        lines += [critic_doc.strip(), '']

    return '\n'.join(lines)


# ──────────────────────────────────────────
# CLI
# ──────────────────────────────────────────

def generate_one(folder: Path, silent: bool = False) -> bool:
    """Generate README.md for a single skill folder. Returns True on success."""
    yaml_path = folder / 'skill.yaml'
    if not yaml_path.exists():
        if not silent:
            print(f"  Skipping {folder.name}: no skill.yaml")
        return False
    try:
        content = build_readme(folder)
        readme_path = folder / 'README.md'
        readme_path.write_text(content, encoding='utf-8')
        if not silent:
            print(f"  Generated: {readme_path}")
        return True
    except Exception as e:
        if not silent:
            print(f"  ERROR generating README for {folder.name}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Auto-generate README.md from skill folder'
    )
    parser.add_argument('path', help='Skill folder or skills root (with --all)')
    parser.add_argument('--all', action='store_true',
                        help='Generate for all skill folders under the given root')
    args = parser.parse_args()

    target = Path(args.path)

    if args.all:
        if not target.is_dir():
            print(f"Error: {target} is not a directory")
            sys.exit(1)
        folders = [d for d in target.iterdir()
                   if d.is_dir() and (d / 'skill.yaml').exists()]
        if not folders:
            print(f"No skill folders found under {target}")
            sys.exit(0)
        success = sum(generate_one(f) for f in folders)
        print(f"\nGenerated {success}/{len(folders)} READMEs")
    else:
        if not target.exists():
            print(f"Error: Not found: {target}")
            sys.exit(1)
        ok = generate_one(target)
        sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
