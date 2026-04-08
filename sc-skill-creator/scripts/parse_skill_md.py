#!/usr/bin/env python3
"""
parse_skill_md.py - Parse a Markdown-format skill file into a skill dict.

The Markdown skill format uses:
- YAML frontmatter for structured metadata
- ## Execution Layer section with one Python code block = code_template
- ## Cognitive Layer section with prose = purpose, H4s under ### Parameter Impact = parameter_impact
- ## Critic Layer section with one Python code block = post_processing_code
- ## Parameter Science Guide with H3 param → H4 condition → bullets

Returns a dict with the same shape as the existing JSON skill schema.
"""

import re
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

try:
    import yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml")


def split_frontmatter(text: str) -> Tuple[str, str]:
    """Split YAML frontmatter from Markdown body.

    Frontmatter is delimited by the first two '---' lines.
    Returns (frontmatter_text, body_text).
    Raises ValueError if frontmatter delimiters are not found.
    """
    lines = text.split('\n')
    if not lines or lines[0].strip() != '---':
        raise ValueError("skill.md must start with '---' YAML frontmatter delimiter")

    end_idx = None
    for i in range(1, len(lines)):
        if lines[i].strip() == '---':
            end_idx = i
            break

    if end_idx is None:
        raise ValueError("skill.md frontmatter closing '---' not found")

    frontmatter_text = '\n'.join(lines[1:end_idx])
    body_text = '\n'.join(lines[end_idx + 1:])
    return frontmatter_text, body_text


def extract_code_block_under_heading(body: str, heading: str) -> str:
    """Extract the first Python code block that appears under the given H2 heading.

    Stops scanning at the next H2 heading (##) or end of document.
    Returns the code block content (without the fence lines).
    Raises ValueError if the heading or code block is not found.
    """
    lines = body.split('\n')

    # Find the heading line
    heading_idx = None
    for i, line in enumerate(lines):
        if line.strip() == heading:
            heading_idx = i
            break

    if heading_idx is None:
        raise ValueError(f"Heading '{heading}' not found in skill.md body")

    # Scan from heading to next H2 or end
    in_code_block = False
    code_lines = []
    found_block = False

    for i in range(heading_idx + 1, len(lines)):
        line = lines[i]

        # Stop at next H2 heading
        if line.startswith('## ') and i > heading_idx:
            break

        if not in_code_block:
            if line.strip().startswith('```python') or line.strip() == '```python':
                in_code_block = True
                found_block = True
                continue
        else:
            if line.strip() == '```':
                in_code_block = False
                break  # Only collect the first code block
            code_lines.append(line)

    if not found_block:
        raise ValueError(f"No Python code block found under heading '{heading}'")

    # Strip leading/trailing blank lines but preserve internal indentation
    code = '\n'.join(code_lines)
    return code.strip('\n')


def extract_first_paragraph_under_heading(body: str, heading: str) -> str:
    """Extract the first non-empty paragraph of text under the given H2 heading.

    Stops at the first sub-heading (###) or next H2 or code block.
    Returns the paragraph text joined into a single string.
    """
    lines = body.split('\n')

    heading_idx = None
    for i, line in enumerate(lines):
        if line.strip() == heading:
            heading_idx = i
            break

    if heading_idx is None:
        return ''

    paragraph_lines = []
    collecting = False

    for i in range(heading_idx + 1, len(lines)):
        line = lines[i]

        # Stop at sub-headings or next H2
        if line.startswith('## ') or line.startswith('### ') or line.startswith('#### '):
            if collecting:
                break
            continue

        # Stop at code fence
        if line.strip().startswith('```'):
            break

        stripped = line.strip()
        if stripped:
            collecting = True
            paragraph_lines.append(stripped)
        elif collecting and paragraph_lines:
            # First blank line after content ends the paragraph
            break

    return ' '.join(paragraph_lines)


def extract_parameter_impact(body: str) -> Dict[str, str]:
    """Parse the ### Parameter Impact section within ## Cognitive Layer.

    Each H4 (####) heading is a parameter name.
    Text following it (until next H4, H3, H2, or code block) is the impact description.
    Returns {param_name: description_string}.
    """
    lines = body.split('\n')

    # Find "### Parameter Impact" heading
    pi_idx = None
    for i, line in enumerate(lines):
        if line.strip() == '### Parameter Impact':
            pi_idx = i
            break

    if pi_idx is None:
        return {}

    result = {}
    current_param = None
    current_lines = []

    for i in range(pi_idx + 1, len(lines)):
        line = lines[i]

        # Stop at next H2 or H3
        if line.startswith('## ') or line.startswith('### '):
            break

        if line.startswith('#### '):
            # Save previous param
            if current_param is not None:
                result[current_param] = ' '.join(current_lines).strip()
            current_param = line[5:].strip()
            current_lines = []
        elif current_param is not None:
            stripped = line.strip()
            if stripped:
                current_lines.append(stripped)

    # Save last param
    if current_param is not None:
        result[current_param] = ' '.join(current_lines).strip()

    return result


def extract_parameter_science_guide(body: str) -> Dict[str, Dict[str, Dict[str, str]]]:
    """Parse the ## Parameter Science Guide section.

    Structure:
        ## Parameter Science Guide
            ### param_name
                #### condition_name
                    - **Adjust**: ...
                    - **Causal chain**: ...
                    - **Expected effect**: ...

    Returns nested dict: {param: {condition: {adjust, causal_chain, expected_effect}}}
    """
    lines = body.split('\n')

    guide_idx = None
    for i, line in enumerate(lines):
        if line.strip() == '## Parameter Science Guide':
            guide_idx = i
            break

    if guide_idx is None:
        return {}

    result = {}
    current_param = None
    current_condition = None
    current_entry: Dict[str, str] = {}

    # Bullet key mappings (case-insensitive prefix match)
    BULLET_KEYS = {
        'adjust': 'adjust',
        'causal chain': 'causal_chain',
        'causal_chain': 'causal_chain',
        'expected effect': 'expected_effect',
        'expected_effect': 'expected_effect',
    }

    def save_condition():
        if current_param and current_condition and current_entry:
            if current_param not in result:
                result[current_param] = {}
            result[current_param][current_condition] = dict(current_entry)

    for i in range(guide_idx + 1, len(lines)):
        line = lines[i]

        # Stop at next H2
        if line.startswith('## ') and i > guide_idx:
            save_condition()
            break

        if line.startswith('### '):
            save_condition()
            current_param = line[4:].strip()
            current_condition = None
            current_entry = {}

        elif line.startswith('#### '):
            save_condition()
            current_condition = line[5:].strip()
            current_entry = {}

        elif line.strip().startswith('- **') and current_condition is not None:
            # Parse bullet: "- **Key**: value"
            m = re.match(r'-\s+\*\*([^*]+)\*\*\s*[:\-]\s*(.*)', line.strip())
            if m:
                raw_key = m.group(1).strip().lower()
                value = m.group(2).strip()
                mapped_key = BULLET_KEYS.get(raw_key)
                if mapped_key:
                    current_entry[mapped_key] = value

    # Save last entry
    save_condition()

    return result


def parse_skill_md(path: str) -> Dict[str, Any]:
    """Parse a Markdown skill file and return a dict matching the JSON skill schema.

    Parameters
    ----------
    path : str
        Path to the skill.md file

    Returns
    -------
    dict
        Skill dict with keys: skill_id, execution_layer, cognitive_layer,
        critic_layer, parameter_science_guide

    Raises
    ------
    ValueError
        If required sections or frontmatter fields are missing
    """
    text = Path(path).read_text(encoding='utf-8')

    frontmatter_text, body = split_frontmatter(text)
    meta = yaml.safe_load(frontmatter_text)

    if not isinstance(meta, dict):
        raise ValueError("YAML frontmatter did not parse to a dict")

    required_meta = ['skill_id', 'required_inputs', 'default_params', 'output_objects',
                     'metrics_to_extract', 'context_params', 'success_thresholds']
    missing = [k for k in required_meta if k not in meta]
    if missing:
        raise ValueError(f"YAML frontmatter missing required keys: {missing}")

    code_template = extract_code_block_under_heading(body, '## Execution Layer')
    post_processing_code = extract_code_block_under_heading(body, '## Critic Layer')
    purpose = extract_first_paragraph_under_heading(body, '## Cognitive Layer')
    parameter_impact = extract_parameter_impact(body)
    parameter_science_guide = extract_parameter_science_guide(body)

    critic_layer: Dict[str, Any] = {
        'metrics_to_extract': meta['metrics_to_extract'],
        'success_thresholds': meta['success_thresholds'],
        'context_params': meta['context_params'],
        'post_processing_code': post_processing_code,
    }
    if 'error_handling' in meta:
        critic_layer['error_handling'] = meta['error_handling']

    return {
        'skill_id': meta['skill_id'],
        'execution_layer': {
            'code_template': code_template,
            'required_inputs': meta['required_inputs'],
            'default_params': meta['default_params'],
            'output_objects': meta['output_objects'],
        },
        'cognitive_layer': {
            'purpose': purpose,
            'parameter_impact': parameter_impact,
        },
        'critic_layer': critic_layer,
        'parameter_science_guide': parameter_science_guide,
    }


if __name__ == '__main__':
    import sys
    import json

    if len(sys.argv) < 2:
        print("Usage: parse_skill_md.py <skill.md>")
        sys.exit(1)

    skill = parse_skill_md(sys.argv[1])
    print(json.dumps(skill, indent=2, ensure_ascii=False))
