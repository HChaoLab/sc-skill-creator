#!/usr/bin/env python3
"""
compile_md_to_json.py - Compile a skill.md to skill.json for downstream tools that require JSON.

Usage:
    python compile_md_to_json.py path/to/skill.md
    python compile_md_to_json.py path/to/skill.md -o path/to/output.json
    python compile_md_to_json.py /path/to/skills_folder --all
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from parse_skill_md import parse_skill_md


def compile_one(md_path: Path, output_path: Path = None) -> Path:
    """Compile a single skill.md to JSON. Returns the output path."""
    if output_path is None:
        output_path = md_path.with_suffix('.json')

    skill = parse_skill_md(str(md_path))
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(skill, f, indent=2, ensure_ascii=False)

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Compile skill.md to skill.json'
    )
    parser.add_argument('path', help='Path to skill.md or skills folder (with --all)')
    parser.add_argument('-o', '--output', help='Output JSON path (single file mode only)')
    parser.add_argument('--all', action='store_true',
                        help='Compile all skill.md files under the given folder')
    args = parser.parse_args()

    target = Path(args.path)

    if args.all:
        if not target.is_dir():
            print(f"Error: {target} is not a directory")
            sys.exit(1)
        found = list(target.rglob('skill.md'))
        if not found:
            print(f"No skill.md files found under {target}")
            sys.exit(0)
        for md_path in found:
            try:
                out = compile_one(md_path)
                print(f"Compiled: {md_path} → {out}")
            except Exception as e:
                print(f"ERROR compiling {md_path}: {e}")
    else:
        if not target.exists():
            print(f"Error: File not found: {target}")
            sys.exit(1)
        output_path = Path(args.output) if args.output else None
        try:
            out = compile_one(target, output_path)
            print(f"Compiled: {target} → {out}")
        except Exception as e:
            print(f"ERROR: {e}")
            sys.exit(1)


if __name__ == '__main__':
    main()
