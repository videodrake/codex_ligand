"""Shared configuration loading and project root resolution."""
import copy
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
    """Legacy project root — returns workflow_a root for new layout."""
    from egfr_pipeline.paths import workflow_a_root
    return workflow_a_root(config)


RESOURCE_DEFAULTS = {
    "scheduler": "pbs",
    "detect_from_env": True,
    "vina": {
        "cpu_per_job": 8,
        "parallel_receptors": 3,
    },
    "ppi": {
        "cpu_per_job": 16,
    },
    "lightdock": {
        "cpu_per_job": 16,
    },
    "phase3": {
        "cpu_per_job": 4,
    },
    "gpu": {
        "enabled": False,
        "engine": "gnina",
        "cpu_per_job": 4,
        "devices": "auto",
    },
}


def _deep_merge_dicts(base: dict, override: dict) -> dict:
    merged = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge_dicts(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def resolve_resource_config(config: dict) -> tuple[dict, list[str]]:
    """Return canonical runtime resources plus compatibility warnings."""
    merged = _deep_merge_dicts(RESOURCE_DEFAULTS, config.get("resources") or {})
    warnings: list[str] = []

    legacy_max_workers = config.get("max_workers")
    legacy_vina = config.get("vina") if isinstance(config.get("vina"), dict) else {}

    legacy_vina_cpu = legacy_vina.get("cpu")
    if legacy_vina_cpu not in (None, ""):
        canonical = merged["vina"].get("cpu_per_job")
        if canonical != int(legacy_vina_cpu):
            warnings.append(
                "Legacy config vina.cpu overrides resources.vina.cpu_per_job; "
                "prefer resources.vina.cpu_per_job going forward."
            )
        merged["vina"]["cpu_per_job"] = int(legacy_vina_cpu)

    legacy_parallel = legacy_vina.get("parallel_receptors", legacy_max_workers)
    if legacy_parallel not in (None, ""):
        canonical = merged["vina"].get("parallel_receptors")
        if canonical != int(legacy_parallel):
            warnings.append(
                "Legacy config vina.parallel_receptors/max_workers overrides "
                "resources.vina.parallel_receptors; prefer resources.vina.parallel_receptors."
            )
        merged["vina"]["parallel_receptors"] = int(legacy_parallel)

    return merged, warnings
