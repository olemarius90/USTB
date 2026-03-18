#!/usr/bin/env python3
"""Auto-generate RST stub files for USTB MATLAB packages.

Scans each +package/ directory for .m files, determines whether each is a
classdef or a function, and writes an RST file with the appropriate Sphinx
autodoc directive.

Usage:
    python generate_api_rst.py
"""

import os
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
API_DIR = Path(__file__).resolve().parent / "api"

PACKAGES = {
    "+uff": {
        "title": "uff Data Classes",
        "description": (
            "The ``uff`` data classes represent the Ultrasound File Format "
            "(UFF) and are used to store ultrasound data, probes, scans, "
            "and related structures. All data classes can be written to and "
            "read from UFF files."
        ),
    },
    "+midprocess": {
        "title": "midprocess",
        "description": (
            "A midprocess performs the delay-and-sum operation, transforming "
            "channel data into beamformed data. The number of output pixels "
            "is determined by a scan grid object.\n\n"
            "| **Input:** ``channel_data`` |rarr| **Output:** ``beamformed_data``"
        ),
    },
    "+postprocess": {
        "title": "postprocess",
        "description": (
            "A postprocess modifies beamformed data. This includes adaptive "
            "beamforming, compounding, and image enhancement techniques.\n\n"
            "| **Input:** ``beamformed_data`` |rarr| **Output:** ``beamformed_data``"
        ),
    },
    "+preprocess": {
        "title": "preprocess",
        "description": (
            "A preprocess modifies channel data before beamforming.\n\n"
            "| **Input:** ``channel_data`` |rarr| **Output:** ``channel_data``"
        ),
    },
    "+tools": {
        "title": "tools Functions",
        "description": (
            "The ``tools`` functions provide utilities for signal "
            "processing, visualization, downloading datasets, and other "
            "helper tasks."
        ),
    },
}

CORE_CLASSES = [
    "uff", "process", "midprocess", "postprocess", "preprocess", "pipeline",
]
CORE_ENUMS = [
    "code", "dimension", "spherical_transmit_delay_model",
]
CORE_FUNCTIONS = [
    "data_path", "ustb_path",
]


def is_classdef(filepath: Path) -> bool:
    """Check if an .m file defines a class (starts with 'classdef')."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("%"):
                    continue
                return stripped.startswith("classdef")
    except OSError:
        pass
    return False


def discover_members(package_dir: Path) -> tuple[list[str], list[str]]:
    """Return sorted lists of (classes, functions) in a package directory."""
    classes = []
    functions = []
    seen = set()
    for mfile in sorted(package_dir.glob("*.m")):
        name = mfile.stem
        if name in seen:
            continue
        seen.add(name)
        if is_classdef(mfile):
            classes.append(name)
        else:
            functions.append(name)
    return classes, functions


def generate_package_rst(package_folder: str, meta: dict) -> str:
    """Generate RST content for a MATLAB package."""
    package_name = package_folder.lstrip("+")
    package_dir = REPO_ROOT / package_folder
    classes, functions = discover_members(package_dir)

    lines = [
        meta["title"],
        "=" * len(meta["title"]),
        "",
        meta["description"],
        "",
    ]

    if classes:
        lines += ["Classes", "-------", ""]
        for cls in classes:
            lines += [
                f".. autoclass:: {package_name}.{cls}",
                "   :members:",
                "",
            ]

    if functions:
        lines += ["Functions", "---------", ""]
        for func in functions:
            lines += [
                f".. autofunction:: {package_name}.{func}",
                "",
            ]

    return "\n".join(lines)


def generate_core_rst() -> str:
    """Generate RST content for root-level core classes and enumerations."""
    lines = [
        "Core Classes",
        "============",
        "",
        "Root-level classes and enumerations that form the foundation of the "
        "USTB framework.",
        "",
        "Base Classes",
        "------------",
        "",
    ]
    for cls in CORE_CLASSES:
        lines += [f".. autoclass:: {cls}", "   :members:", ""]

    lines += ["Enumerations", "------------", ""]
    for cls in CORE_ENUMS:
        lines += [f".. autoclass:: {cls}", "   :members:", ""]

    lines += ["Utilities", "---------", ""]
    for func in CORE_FUNCTIONS:
        lines += [f".. autofunction:: {func}", ""]

    return "\n".join(lines)


def main():
    API_DIR.mkdir(parents=True, exist_ok=True)

    for folder, meta in PACKAGES.items():
        package_name = folder.lstrip("+")
        rst_path = API_DIR / f"{package_name}.rst"
        content = generate_package_rst(folder, meta)
        rst_path.write_text(content, encoding="utf-8")
        print(f"Generated {rst_path.relative_to(REPO_ROOT)}")

    core_path = API_DIR / "core.rst"
    core_path.write_text(generate_core_rst(), encoding="utf-8")
    print(f"Generated {core_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
