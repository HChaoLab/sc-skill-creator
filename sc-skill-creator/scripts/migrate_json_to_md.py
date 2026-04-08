#!/usr/bin/env python3
"""
migrate_json_to_md.py - Convert an existing skill.json to the Markdown-first skill.md format.

Usage:
    python migrate_json_to_md.py path/to/skill.json
    python migrate_json_to_md.py path/to/skill.json --remove-json   # delete JSON after verified round-trip
    python migrate_json_to_md.py /path/to/skills_folder --all        # migrate all skill.json in subfolders

The script performs a round-trip validation before removing JSON:
    JSON → write skill.md → parse_skill_md() → compare dicts field-by-field
If any field differs, the JSON is kept and a warning is printed.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml")

# Allow running from any directory
sys.path.insert(0, str(Path(__file__).parent))
from parse_skill_md import parse_skill_md


# ──────────────────────────────────────────────
# Assembly helpers
# ──────────────────────────────────────────────

def _yaml_block(value: Any, key: str) -> str:
    """Render a single YAML key as a block scalar for multi-line strings,
    or as a regular YAML field for simple values."""
    if isinstance(value, str) and ('\n' in value or "'" in value or '"' in value):
        # Use literal block scalar to avoid quoting issues
        dumped = yaml.dump({key: value}, default_flow_style=False,
                           allow_unicode=True, sort_keys=False)
        return dumped.rstrip()
    return yaml.dump({key: value}, default_flow_style=False,
                     allow_unicode=True, sort_keys=False).rstrip()


def build_frontmatter(skill: Dict[str, Any]) -> str:
    """Build the YAML frontmatter block from a skill dict."""
    el = skill.get('execution_layer', {})
    cl = skill.get('critic_layer', {})

    meta = {
        'skill_id': skill['skill_id'],
        'required_inputs': el.get('required_inputs', []),
        'default_params': el.get('default_params', {}),
        'output_objects': el.get('output_objects', []),
        'metrics_to_extract': cl.get('metrics_to_extract', []),
        'context_params': cl.get('context_params', []),
        'success_thresholds': cl.get('success_thresholds', ''),
    }
    # Only include error_handling if present in the original
    if 'error_handling' in cl:
        meta['error_handling'] = cl['error_handling']

    dumped = yaml.dump(meta, default_flow_style=False, allow_unicode=True,
                       sort_keys=False)
    return f"---\n{dumped.rstrip()}\n---"


def build_parameter_impact_section(parameter_impact: Dict[str, str]) -> str:
    """Build the ### Parameter Impact sub-section."""
    if not parameter_impact:
        return ''
    lines = ['\n### Parameter Impact']
    for param, description in parameter_impact.items():
        lines.append(f'\n#### {param}\n')
        lines.append(description)
    return '\n'.join(lines)


def build_parameter_science_guide(guide: Dict[str, Any]) -> str:
    """Build the ## Parameter Science Guide section."""
    if not guide:
        return ''
    lines = ['\n## Parameter Science Guide']
    for param, conditions in guide.items():
        lines.append(f'\n### {param}')
        if not isinstance(conditions, dict):
            continue
        for condition, entry in conditions.items():
            lines.append(f'\n#### {condition}\n')
            if isinstance(entry, dict):
                adjust = entry.get('adjust', '')
                causal = entry.get('causal_chain', entry.get('reason', ''))
                effect = entry.get('expected_effect', '')
                if adjust:
                    lines.append(f'- **Adjust**: {adjust}')
                if causal:
                    lines.append(f'- **Causal chain**: {causal}')
                if effect:
                    lines.append(f'- **Expected effect**: {effect}')
    return '\n'.join(lines)


def skill_to_markdown(skill: Dict[str, Any]) -> str:
    """Convert a skill dict to the Markdown skill format."""
    el = skill.get('execution_layer', {})
    cog = skill.get('cognitive_layer', {})
    cl = skill.get('critic_layer', {})
    guide = skill.get('parameter_science_guide', {})

    skill_id = skill['skill_id']
    code_template = el.get('code_template', '')
    post_processing_code = cl.get('post_processing_code', '')
    purpose = cog.get('purpose', '')
    parameter_impact = cog.get('parameter_impact', {})

    frontmatter = build_frontmatter(skill)
    param_impact_section = build_parameter_impact_section(parameter_impact)
    psg_section = build_parameter_science_guide(guide)

    lines = [
        frontmatter,
        '',
        f'# Skill: {skill_id}',
        '',
        '## Execution Layer',
        '',
        '```python',
        code_template,
        '```',
        '',
        '## Cognitive Layer',
        '',
        purpose,
        param_impact_section,
        '',
        '## Critic Layer',
        '',
        '```python',
        post_processing_code,
        '```',
        psg_section,
    ]

    return '\n'.join(lines) + '\n'


# ──────────────────────────────────────────────
# Round-trip validation
# ──────────────────────────────────────────────

def _normalize(value: Any) -> Any:
    """Normalize values for comparison: strip strings, sort list keys."""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        return [_normalize(v) for v in value]
    if isinstance(value, dict):
        return {k: _normalize(v) for k, v in value.items()}
    return value


def _diff_dicts(orig: Dict, parsed: Dict, path: str = '') -> list:
    """Return list of field differences between two dicts."""
    diffs = []
    all_keys = set(orig.keys()) | set(parsed.keys())
    for key in sorted(all_keys):
        full_path = f"{path}.{key}" if path else key
        if key not in orig:
            diffs.append(f"Extra key in parsed: {full_path}")
        elif key not in parsed:
            diffs.append(f"Missing key in parsed: {full_path}")
        else:
            o, p = _normalize(orig[key]), _normalize(parsed[key])
            if isinstance(o, dict) and isinstance(p, dict):
                diffs.extend(_diff_dicts(o, p, full_path))
            elif o != p:
                diffs.append(f"Value mismatch at {full_path}")
    return diffs


def round_trip_validate(md_path: str, original_skill: Dict) -> list:
    """Parse skill.md back to a dict and compare with the original. Returns list of diffs."""
    try:
        parsed = parse_skill_md(md_path)
    except Exception as e:
        return [f"parse_skill_md failed: {e}"]
    return _diff_dicts(original_skill, parsed)


# ──────────────────────────────────────────────
# Migration logic
# ──────────────────────────────────────────────

def migrate_one(json_path: Path, remove_json: bool = False) -> bool:
    """Migrate a single skill.json to skill.md. Returns True on success."""
    print(f"Migrating: {json_path}")

    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            skill = json.load(f)
    except Exception as e:
        print(f"  ERROR reading JSON: {e}")
        return False

    md_path = json_path.parent / 'skill.md'
    md_content = skill_to_markdown(skill)

    md_path.write_text(md_content, encoding='utf-8')
    print(f"  Written: {md_path}")

    # Round-trip validation
    diffs = round_trip_validate(str(md_path), skill)
    if diffs:
        print(f"  WARNING: Round-trip validation found {len(diffs)} difference(s):")
        for d in diffs:
            print(f"    - {d}")
        print(f"  Keeping {json_path.name} (skill.md written but may need manual review)")
        return False

    print(f"  Round-trip validation PASSED")

    if remove_json:
        json_path.unlink()
        print(f"  Removed: {json_path}")

    return True


def migrate_all(skills_root: Path, remove_json: bool = False) -> None:
    """Migrate all skill.json files found in subdirectories of skills_root."""
    found = list(skills_root.rglob('skill.json'))
    if not found:
        print(f"No skill.json files found under {skills_root}")
        return

    success = 0
    for json_path in found:
        if migrate_one(json_path, remove_json):
            success += 1

    print(f"\nMigrated {success}/{len(found)} skills successfully")


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Migrate skill.json to Markdown-first skill.md format'
    )
    parser.add_argument('path', help='Path to skill.json or skills folder (with --all)')
    parser.add_argument('--all', action='store_true',
                        help='Migrate all skill.json files under the given folder')
    parser.add_argument('--remove-json', action='store_true',
                        help='Remove skill.json after successful round-trip validation')
    args = parser.parse_args()

    target = Path(args.path)

    if args.all:
        if not target.is_dir():
            print(f"Error: {target} is not a directory")
            sys.exit(1)
        migrate_all(target, args.remove_json)
    else:
        if not target.exists():
            print(f"Error: File not found: {target}")
            sys.exit(1)
        success = migrate_one(target, args.remove_json)
        sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
