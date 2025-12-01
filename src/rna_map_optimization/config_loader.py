"""Configuration file loader for optimization settings."""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional


def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """Load optimization configuration from YAML file.
    
    Args:
        config_path: Path to config file. If None, looks for config/optimization_config.yml
                    relative to package root or current directory.
        
    Returns:
        Dictionary with configuration settings
        
    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is invalid YAML
    """
    if config_path is None:
        # Try to find config file relative to package
        package_root = Path(__file__).parent.parent.parent
        config_path = package_root / "config" / "optimization_config.yml"
        
        # If not found, try current directory
        if not config_path.exists():
            config_path = Path("config/optimization_config.yml")
    
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found: {config_path}. "
            f"Please create config/optimization_config.yml or specify a custom path."
        )
    
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    return config


def get_case_config(config: Dict[str, Any], case_name: str) -> Dict[str, Any]:
    """Get configuration for a specific case, applying overrides if present.
    
    Args:
        config: Full configuration dictionary
        case_name: Name of the case to get config for
        
    Returns:
        Configuration dictionary with case-specific overrides applied
    """
    # Start with base optimization config
    case_config = config.get("optimization", {}).copy()
    
    # Apply case-specific overrides if they exist
    case_overrides = config.get("case_overrides", {})
    if case_name in case_overrides:
        case_config.update(case_overrides[case_name])
    
    return case_config


def get_parameter_ranges(config: Dict[str, Any]) -> Dict[str, Any]:
    """Extract parameter search ranges from configuration.
    
    Args:
        config: Full configuration dictionary
        
    Returns:
        Dictionary of parameter ranges for Optuna
    """
    return config.get("optimization", {}).get("parameters", {})

