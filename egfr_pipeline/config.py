"""Shared configuration loading and project root resolution."""
import json
from pathlib import Path

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


def load_config(path: str) -> dict:
    config_path = Path(path)
    if config_path.suffix in (".yaml", ".yml"):
        if not HAS_YAML:
            raise RuntimeError("pyyaml is required to load YAML config files.")
        with open(config_path, encoding="utf-8") as f:
            return yaml.safe_load(f)
    with open(config_path, encoding="utf-8") as f:
        return json.load(f)


def save_config(config: dict, path: str):
    config_path = Path(path)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    if config_path.suffix in (".yaml", ".yml"):
        if not HAS_YAML:
            raise RuntimeError("pyyaml is required to save YAML config files.")
        with open(config_path, "w", encoding="utf-8") as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    else:
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(config, f, indent=2, ensure_ascii=False)


def project_root_from_config(config: dict) -> Path:
    output_root = Path(config.get("output_root", "./output"))
    project_name = config.get("project_name")
    return output_root / project_name if project_name else output_root
