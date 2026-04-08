#!/usr/bin/env python3
"""
Validate a generated sc-skill-creator skill JSON structure.
Checks for production robustness requirements:
1. Memory-First Execution pattern
2. Parameter Safety with default_params
3. Context-Aware Thresholds in critic_layer
4. Metadata Footprinting (analysis_history)
"""

import json
import sys
from pathlib import Path
from typing import Dict, Any, List

sys.path.insert(0, str(Path(__file__).parent))
from parse_skill_md import parse_skill_md, split_frontmatter, extract_code_block_under_heading
from load_skill_folder import load_skill_folder, REQUIRED_YAML_KEYS


REQUIRED_FIELDS = {
    "skill_id": str,
    "execution_layer": dict,
    "cognitive_layer": dict,
    "critic_layer": dict,
    "parameter_science_guide": dict
}

EXECUTION_LAYER_FIELDS = {
    "code_template": str,
    "required_inputs": list,
    "default_params": dict,
    "output_objects": list
}

COGNITIVE_LAYER_FIELDS = {
    "purpose": str,
    "parameter_impact": dict
}

CRITIC_LAYER_FIELDS = {
    "metrics_to_extract": list,
    "success_thresholds": str,
    "context_params": list,
    "post_processing_code": str
}
# error_handling is optional


def validate_skill(skill_json: Dict[str, Any]) -> List[str]:
    """Validate a skill JSON structure. Returns list of errors."""
    errors = []
    
    for field, expected_type in REQUIRED_FIELDS.items():
        if field not in skill_json:
            errors.append(f"Missing required field: {field}")
        elif not isinstance(skill_json[field], expected_type):
            errors.append(f"Field '{field}' should be {expected_type.__name__}, got {type(skill_json[field]).__name__}")
    
    if "execution_layer" in skill_json:
        el = skill_json["execution_layer"]
        for field in EXECUTION_LAYER_FIELDS:
            if field not in el:
                errors.append(f"Missing execution_layer field: {field}")
            elif not isinstance(el[field], EXECUTION_LAYER_FIELDS[field]):
                errors.append(f"execution_layer.{field} should be {EXECUTION_LAYER_FIELDS[field].__name__}")
        
        code = el.get("code_template", "")
        
        if "isinstance(input_data" not in code and "isinstance(input_data, ad.AnnData)" not in code:
            errors.append("code_template missing Memory-First pattern: isinstance(input_data, ad.AnnData) check")
        
        if "default_params" not in code and "{**default_params" not in code:
            errors.append("code_template missing Parameter Safety pattern: {**default_params, **agent_params} merge")
        
        if "current_params[" not in code and "current_params.get(" not in code:
            errors.append("code_template should use current_params['param'] instead of bare {{placeholder}}")
        
        if "analysis_history" not in code:
            errors.append("code_template missing Metadata Footprinting: adata.uns['analysis_history'] append")
    
    if "cognitive_layer" in skill_json:
        if "purpose" not in skill_json["cognitive_layer"]:
            errors.append("Missing cognitive_layer.purpose")
        if "parameter_impact" not in skill_json["cognitive_layer"]:
            errors.append("Missing cognitive_layer.parameter_impact")
    
    if "critic_layer" in skill_json:
        cl = skill_json["critic_layer"]
        for field in CRITIC_LAYER_FIELDS:
            if field not in cl:
                errors.append(f"Missing critic_layer field: {field}")
        
        if "post_processing_code" in cl:
            code = cl["post_processing_code"]
            if "def critic_post_process" not in code:
                errors.append("post_processing_code should define a critic_post_process function")
            if "metrics" not in code.lower() or "warnings" not in code.lower():
                errors.append("post_processing_code should extract metrics and generate warnings")
            if "context" not in code:
                errors.append("post_processing_code missing context parameter for protocol-specific thresholds")
            
            # Check context-aware thresholds: accept either in success_thresholds expression
            # OR in post_processing_code (where t['key'] is more commonly used)
            thresholds_expr = cl.get("success_thresholds", "")
            has_context_thresholds = (
                "t['" in thresholds_expr or 't["' in thresholds_expr or
                "t['" in code or 't["' in code
            )
            if not has_context_thresholds:
                errors.append("success_thresholds or post_processing_code should use context thresholds (t['min_clusters'], etc.)")
    
    if "parameter_science_guide" in skill_json:
        guide = skill_json["parameter_science_guide"]
        for param, actions in guide.items():
            if not isinstance(actions, dict):
                continue
            for condition, instruction in actions.items():
                if isinstance(instruction, dict):
                    if "adjust" not in instruction:
                        errors.append(f"parameter_science_guide.{param}.{condition} missing 'adjust' field")
                    if "causal_chain" not in instruction and "reason" not in instruction:
                        errors.append(f"parameter_science_guide.{param}.{condition} missing causal explanation")
    
    return errors


def _detect_language(folder_path: Path) -> str:
    """Detect skill language from skill.yaml, falling back to file presence."""
    yaml_path = folder_path / 'skill.yaml'
    if yaml_path.exists():
        try:
            import yaml as _yaml
            with open(yaml_path, encoding='utf-8') as f:
                meta = _yaml.safe_load(f)
            lang = (meta or {}).get('language', '').lower()
            if lang in ('python', 'r'):
                return lang
        except Exception:
            pass
    if (folder_path / 'execution.R').exists():
        return 'r'
    return 'python'


def validate_yaml_folder(folder: str) -> List[str]:
    """Check that a skill folder has valid skill.yaml, execution script, and critic script."""
    errors = []
    folder_path = Path(folder)

    yaml_path = folder_path / 'skill.yaml'
    if not yaml_path.exists():
        errors.append("Missing required file: skill.yaml")

    language = _detect_language(folder_path)
    exec_name   = 'execution.R' if language == 'r' else 'execution.py'
    critic_name = 'critic.R'    if language == 'r' else 'critic.py'

    for fname in (exec_name, critic_name):
        if not (folder_path / fname).exists():
            errors.append(f"Missing required file: {fname}")

    if yaml_path.exists():
        try:
            import yaml as _yaml
            with open(yaml_path, encoding='utf-8') as f:
                meta = _yaml.safe_load(f)
            if not isinstance(meta, dict):
                errors.append("skill.yaml did not parse to a dict")
            else:
                for k in REQUIRED_YAML_KEYS:
                    if k not in meta:
                        errors.append(f"skill.yaml missing required key: {k}")
                # Validate parameter_science_guide feedback keys match conditions
                guide = meta.get('parameter_science_guide', {})
                if not isinstance(guide, dict):
                    errors.append("parameter_science_guide must be a dict")
        except Exception as e:
            errors.append(f"skill.yaml parse error: {e}")

    exec_path = folder_path / exec_name
    if exec_path.exists():
        code = exec_path.read_text(encoding='utf-8')
        if language == 'r':
            _validate_r_execution(code, exec_name, errors)
        else:
            _validate_py_execution(code, exec_name, errors)

    critic_path = folder_path / critic_name
    if critic_path.exists():
        code = critic_path.read_text(encoding='utf-8')
        if language == 'r':
            _validate_r_critic(code, critic_name, errors)
        else:
            _validate_py_critic(code, critic_name, errors)

    return errors


def _validate_py_execution(code: str, fname: str, errors: List[str]) -> None:
    if 'isinstance(input_data' not in code:
        errors.append(f"{fname} missing Memory-First pattern: isinstance(input_data, ad.AnnData)")
    if '{**default_params' not in code:
        errors.append(f"{fname} missing Parameter Safety pattern: {{**default_params, **agent_params}}")
    if 'analysis_history' not in code:
        errors.append(f"{fname} missing Metadata Footprinting: analysis_history append")


def _validate_r_execution(code: str, fname: str, errors: List[str]) -> None:
    if 'is.character(input_data)' not in code:
        errors.append(f"{fname} missing Memory-First pattern: is.character(input_data) check")
    if 'modifyList' not in code:
        errors.append(f"{fname} missing Parameter Safety pattern: modifyList(default_params, ...)")
    if 'analysis_history' not in code:
        errors.append(f"{fname} missing Metadata Footprinting: analysis_history append")


def _validate_py_critic(code: str, fname: str, errors: List[str]) -> None:
    if 'def critic_post_process' not in code:
        errors.append(f"{fname} must define critic_post_process()")
    if 'context' not in code:
        errors.append(f"{fname} missing context parameter for protocol-specific thresholds")
    if "'feedback'" not in code and '"feedback"' not in code:
        errors.append(f"{fname} should return a 'feedback' list (maps to parameter_science_guide keys)")


def _validate_r_critic(code: str, fname: str, errors: List[str]) -> None:
    if 'critic_post_process' not in code:
        errors.append(f"{fname} must define critic_post_process()")
    if 'context' not in code:
        errors.append(f"{fname} missing context parameter for protocol-specific thresholds")
    if '"feedback"' not in code and "'feedback'" not in code and 'feedback' not in code:
        errors.append(f"{fname} should return a 'feedback' vector (maps to parameter_science_guide keys)")


def validate_markdown_structure(path: str) -> List[str]:
    """Check that a skill.md has valid frontmatter and required code sections."""
    errors = []
    try:
        text = Path(path).read_text(encoding='utf-8')
    except Exception as e:
        return [f"Cannot read file: {e}"]

    try:
        frontmatter_text, body = split_frontmatter(text)
    except ValueError as e:
        return [f"Frontmatter parse error: {e}"]

    try:
        import yaml
        meta = yaml.safe_load(frontmatter_text)
        if not isinstance(meta, dict):
            errors.append("YAML frontmatter did not parse to a dict")
    except Exception as e:
        errors.append(f"YAML parse error: {e}")
        return errors

    required_keys = ['skill_id', 'required_inputs', 'default_params', 'output_objects',
                     'metrics_to_extract', 'context_params', 'success_thresholds']
    for k in required_keys:
        if k not in meta:
            errors.append(f"YAML frontmatter missing key: {k}")

    for section in ['## Execution Layer', '## Critic Layer']:
        try:
            extract_code_block_under_heading(body, section)
        except ValueError:
            errors.append(f"Missing Python code block under '{section}'")

    return errors


def main():
    if len(sys.argv) < 2:
        print("Usage: validate_skill.py <skill_folder | skill.yaml | skill.md | skill.json>")
        sys.exit(1)

    skill_path = sys.argv[1]
    path_obj = Path(skill_path)

    if not path_obj.exists():
        print(f"Error: Not found: {skill_path}")
        sys.exit(1)

    # Detect format
    is_yaml_folder = path_obj.is_dir() or path_obj.suffix == '.yaml'
    is_markdown = path_obj.suffix == '.md'
    folder_path = path_obj if path_obj.is_dir() else path_obj.parent

    if is_yaml_folder:
        struct_errors = validate_yaml_folder(str(folder_path))
        if struct_errors:
            print("Folder structure validation FAILED:")
            for e in struct_errors:
                print(f"  - {e}")
            sys.exit(1)
        try:
            skill_json = load_skill_folder(str(folder_path))
        except Exception as e:
            print(f"Error: Failed to load skill folder: {e}")
            sys.exit(1)
        fmt = "YAML folder"

    elif is_markdown:
        # First validate markdown structure
        md_errors = validate_markdown_structure(skill_path)
        if md_errors:
            print("Markdown structure validation FAILED:")
            for e in md_errors:
                print(f"  - {e}")
            sys.exit(1)

        try:
            skill_json = parse_skill_md(skill_path)
        except Exception as e:
            print(f"Error: Failed to parse skill.md: {e}")
            sys.exit(1)
    else:
        try:
            with open(skill_path, 'r') as f:
                skill_json = json.load(f)
        except FileNotFoundError:
            print(f"Error: File not found: {skill_path}")
            sys.exit(1)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON: {e}")
            sys.exit(1)
        fmt = "JSON"

    errors = validate_skill(skill_json)

    if errors:
        print("Validation FAILED:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)
    else:
        print(f"Validation PASSED: Skill ({fmt}) is well-formed")
        print(f"  skill_id: {skill_json.get('skill_id', 'N/A')}")
        print(f"  execution_layer parameters: {len(skill_json.get('execution_layer', {}).get('default_params', {}))}")
        print(f"  cognitive_layer parameters: {len(skill_json.get('cognitive_layer', {}).get('parameter_impact', {}))}")
        print(f"  critic_layer metrics: {len(skill_json.get('critic_layer', {}).get('metrics_to_extract', []))}")
        print(f"  parameter_science_guide entries: {len(skill_json.get('parameter_science_guide', {}))}")
        print("\nProduction Robustness Checks:")
        el = skill_json.get('execution_layer', {})
        code = el.get('code_template', '')
        print(f"  [x] Memory-First pattern: {'isinstance' in code}")
        print(f"  [x] Parameter Safety: {'{**default_params' in code or '{**default_params' in code}")
        print(f"  [x] Metadata Footprinting: {'analysis_history' in code}")
        cl = skill_json.get('critic_layer', {})
        print(f"  [x] Context params: {cl.get('context_params', [])}")
        has_thresholds = 'thresholds[' in cl.get('post_processing_code', '') or 'thresholds.get' in cl.get('post_processing_code', '')
        print(f"  [x] Context thresholds in post_processing: {has_thresholds}")


if __name__ == "__main__":
    main()
