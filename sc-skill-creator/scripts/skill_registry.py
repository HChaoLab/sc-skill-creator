#!/usr/bin/env python3
"""
SkillRegistry - Scan and manage sc-skill-creator generated skills.

Each Skill lives in its own folder (folder name = tool name).
Folder must contain skill.md (preferred) or skill.json (legacy fallback).
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any

sys.path.insert(0, str(Path(__file__).parent))
from parse_skill_md import parse_skill_md
from load_skill_folder import load_skill_folder
from generate_readme import generate_one


class SkillRegistry:
    """
    Registry for scanning and accessing sc-skill-creator generated skills.
    
    Storage Rule (priority order):
        /skills/
            scanpy_leiden/          ← preferred: separate scripts
                skill.yaml
                execution.py
                critic.py
                README.md           ← auto-generated
            scanpy_filter_cells/    ← legacy: markdown
                skill.md
            scanpy_pca/             ← legacy: JSON
                skill.json
    
    Usage:
        registry = SkillRegistry('/path/to/skills')
        registry.scan()  # Initial full scan
        
        # Agent tool selection
        for skill in registry.list_skills():
            print(f"{skill['skill_id']}: {skill['purpose']}")
        
        # Get full skill for SkillRunner
        skill = registry.get_skill('scanpy_filter_cells')
        
        # Dynamic registration (when sc-skill-creator generates new skill)
        registry.register('/path/to/new_skill/skill.json')
    """
    
    def __init__(self, skills_root: str):
        """
        Initialize registry with skills root folder.
        
        Parameters:
        ----------
        skills_root : str
            Path to /skills root directory containing skill subfolders
        """
        self.skills_root = Path(skills_root)
        self._index: Dict[str, Dict[str, str]] = {}
    
    def scan(self) -> int:
        """
        Deep scan skills_root, find all valid skill folders.
        A valid folder must contain skill.json.
        
        Returns:
        --------
        int : Number of skills found and indexed
        """
        self._index.clear()
        
        if not self.skills_root.exists():
            raise FileNotFoundError(f"Skills root not found: {self.skills_root}")
        
        for subfolder in self.skills_root.iterdir():
            if not subfolder.is_dir():
                continue

            skill_yaml_path = subfolder / "skill.yaml"
            skill_md_path   = subfolder / "skill.md"
            skill_json_path = subfolder / "skill.json"

            # Priority: skill.yaml (separate scripts) > skill.md > skill.json
            if skill_yaml_path.exists():
                source_path = skill_yaml_path
                source_format = 'yaml'
            elif skill_md_path.exists():
                source_path = skill_md_path
                source_format = 'markdown'
            elif skill_json_path.exists():
                source_path = skill_json_path
                source_format = 'json'
            else:
                continue

            try:
                if source_format == 'yaml':
                    skill_data = load_skill_folder(str(subfolder))
                elif source_format == 'markdown':
                    skill_data = parse_skill_md(str(source_path))
                else:
                    with open(source_path, 'r', encoding='utf-8') as f:
                        skill_data = json.load(f)

                skill_id = skill_data.get('skill_id')
                if not skill_id:
                    print(f"Warning: {source_path} missing skill_id, skipping")
                    continue

                self._index[skill_id] = {
                    'skill_id': skill_id,
                    'purpose': skill_data.get('cognitive_layer', {}).get('purpose', 'N/A'),
                    'folder': str(subfolder),
                    'skill_source_path': str(source_path.resolve()),
                    'source_format': source_format,
                }

                # Auto-generate README.md for yaml-format skills if missing
                if source_format == 'yaml' and not (subfolder / 'README.md').exists():
                    generate_one(subfolder, silent=True)

            except Exception as e:
                print(f"Warning: Error reading {source_path}: {e}")
        
        return len(self._index)
    
    def register(self, skill_path: str) -> bool:
        """
        Dynamically register a new skill without rescanning.
        Accepts either skill.md or skill.json.

        Parameters:
        ----------
        skill_path : str
            Absolute path to skill.md or skill.json

        Returns:
        --------
        bool : True if registered successfully, False otherwise
        """
        skill_path = Path(skill_path)

        if not skill_path.exists():
            print(f"Error: File not found: {skill_path}")
            return False

        if skill_path.suffix == '.yaml':
            source_format = 'yaml'
        elif skill_path.suffix == '.md':
            source_format = 'markdown'
        else:
            source_format = 'json'

        try:
            if source_format == 'yaml':
                skill_data = load_skill_folder(str(skill_path.parent))
            elif source_format == 'markdown':
                skill_data = parse_skill_md(str(skill_path))
            else:
                with open(skill_path, 'r', encoding='utf-8') as f:
                    skill_data = json.load(f)

            skill_id = skill_data.get('skill_id')
            if not skill_id:
                print(f"Error: skill file missing skill_id")
                return False

            self._index[skill_id] = {
                'skill_id': skill_id,
                'purpose': skill_data.get('cognitive_layer', {}).get('purpose', 'N/A'),
                'folder': str(skill_path.parent),
                'skill_source_path': str(skill_path.resolve()),
                'source_format': source_format,
            }

            print(f"Registered: {skill_id} (from {source_format})")

            # Auto-generate README.md for new yaml-format skills
            if source_format == 'yaml':
                generate_one(skill_path.parent, silent=True)

            return True

        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in {skill_path}: {e}")
            return False
        except Exception as e:
            print(f"Error registering {skill_path}: {e}")
            return False
    
    def unregister(self, skill_id: str) -> bool:
        """
        Remove a skill from the registry.
        
        Parameters:
        ----------
        skill_id : str
        
        Returns:
        --------
        bool : True if removed, False if not found
        """
        if skill_id in self._index:
            del self._index[skill_id]
            return True
        return False
    
    def get_skill(self, skill_id: str) -> Optional[Dict[str, Any]]:
        """
        Get full skill object by skill_id.
        
        Parameters:
        ----------
        skill_id : str
            The skill identifier
        
        Returns:
        --------
        dict or None : Full skill JSON object, or None if not found
        """
        if not self._index:
            self.scan()
        
        info = self._index.get(skill_id)
        if not info:
            return None

        try:
            source_path = info['skill_source_path']
            fmt = info.get('source_format')
            if fmt == 'yaml':
                return load_skill_folder(info['folder'])
            elif fmt == 'markdown':
                return parse_skill_md(source_path)
            else:
                with open(source_path, 'r', encoding='utf-8') as f:
                    return json.load(f)
        except Exception as e:
            print(f"Error reading skill file: {e}")
            return None
    
    def list_skills(self) -> List[Dict[str, str]]:
        """
        Get list of all available skills for Agent tool selection.
        
        Returns:
        --------
        list : List of dicts with skill_id and purpose
        """
        if not self._index:
            self.scan()
        
        return [
            {
                'skill_id': info['skill_id'],
                'purpose': info['purpose']
            }
            for info in self._index.values()
        ]
    
    def list_skills_detailed(self) -> List[Dict[str, Any]]:
        """
        Get detailed list of all skills with full metadata.
        
        Returns:
        --------
        list : List of skill info dicts
        """
        if not self._index:
            self.scan()
        
        result = []
        for info in self._index.values():
            skill = self.get_skill(info['skill_id'])
            if skill:
                result.append({
                    'skill_id': info['skill_id'],
                    'purpose': info['purpose'],
                    'folder': info['folder'],
                    'default_params': skill.get('execution_layer', {}).get('default_params', {}),
                    'metrics': skill.get('critic_layer', {}).get('metrics_to_extract', []),
                    'context_params': skill.get('critic_layer', {}).get('context_params', [])
                })
        return result
    
    def search(self, query: str) -> List[Dict[str, str]]:
        """
        Search skills by keyword in skill_id or purpose.
        
        Parameters:
        ----------
        query : str
            Search term
        
        Returns:
        --------
        list : Matching skills
        """
        if not self._index:
            self.scan()
        
        query_lower = query.lower()
        return [
            {
                'skill_id': info['skill_id'],
                'purpose': info['purpose']
            }
            for info in self._index.values()
            if query_lower in info['skill_id'].lower() or query_lower in info['purpose'].lower()
        ]
    
    def get_toolbox_for_agent(self) -> str:
        """
        Generate a formatted toolbox清单 for Agent.
        
        Returns:
        --------
        str : Human-readable list of available tools
        """
        if not self._index:
            self.scan()
        
        lines = ["Available Tools:", ""]
        
        for info in self._index.values():
            lines.append(f"  [{info['skill_id']}]")
            lines.append(f"    Purpose: {info['purpose']}")
            lines.append("")
        
        return "\n".join(lines)
    
    @property
    def skill_ids(self) -> List[str]:
        """Get list of all registered skill IDs."""
        if not self._index:
            self.scan()
        return list(self._index.keys())
    
    def __len__(self) -> int:
        """Number of skills in registry."""
        return len(self._index)
    
    def __repr__(self) -> str:
        return f"SkillRegistry(root={self.skills_root}, skills={len(self)})"


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        root = sys.argv[1]
    else:
        root = "/home/rstudio/scAgentSkills"
    
    registry = SkillRegistry(root)
    registry.scan()
    
    print(f"=== SkillRegistry: {registry} ===\n")
    print(registry.get_toolbox_for_agent())
