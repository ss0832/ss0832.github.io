---
title:  【Python】 How to make pyproject.toml file
published: 2025-12-06
description: ""
tags: [toml]
category: Python
draft: false
---
Last updated: 2025-12-06

This article explains the basics of creating a `pyproject.toml` file for packaging Python projects and installing them via `pip`.

## What is pyproject.toml?

`pyproject.toml` is the standard configuration file for Python packaging (defined in PEP 518 and PEP 621). It replaces the legacy `setup.py` and `requirements.txt` approach, centralizing build dependencies and project metadata in a single, declarative file.

By configuring this file correctly, you enable users to install your package directly using `pip`, even from a Git repository.

## Basic Structure

Here is a complete example of a `pyproject.toml` file.

```toml name=pyproject.toml
[build-system]
requires = ["setuptools", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "my-awesome-package"
version = "0.1.0"
description = "A sample python package for demonstration."
readme = "README.md"
requires-python = ">=3.10"
license = {text = "MIT"}
authors = [
    {name = "Your Name", email = "your.email@example.com"}
]
dependencies = [
    "numpy>=1.24.0",
    "pandas>=2.0.0"
]

[project.scripts]
my-cli = "my_package.main:run"

[tool.setuptools.packages.find]
where = ["."]
include = ["my_package*"]
```

## Detailed Explanation of Sections

### 1. [build-system]
This section defines *how* the package is built.

*   **requires**:
    *   Specifies the build tools needed *before* installation starts.
    *   Commonly includes `"setuptools"` (the library to facilitate packaging) and `"wheel"` (to create binary wheel files).
*   **build-backend**:
    *   Tells the frontend (like pip) which Python object to use to perform the build.
    *   For setuptools, this is usually `"setuptools.build_meta"`.

### 2. [project]
This contains the core metadata of your package.

*   **name**: The name of your package on PyPI or when installed.
*   **version**: The semantic version string (e.g., "1.0.0").
*   **requires-python**: Specifies the supported Python versions.
*   **dependencies**: A list of runtime dependencies (libraries required for your code to run). `pip` will automatically install these.

### 3. [project.scripts]
This allows you to create command-line executables.

*   Format: `command_name = "package.module:function"`
*   In the example above, typing `my-cli` in the terminal will execute the `run` function inside `my_package/main.py`.

### 4. [tool.setuptools.packages.find]
This configuration helps setuptools automatically discover your source code directories.

*   **where**: The root directory to search (usually `["."]`).
*   **include**: Patterns to match package directories (e.g., `["my_package*"]`).

## Hands-on: Installing via Git

Once you have pushed your code and this `pyproject.toml` to a GitHub repository, anyone can install your package using `pip` without needing to publish it to PyPI.

### Command
```bash
pip install git+https://github.com/username/repository_name.git
```

### Specific Branch or Tag
To install a specific branch (e.g., `develop`) or tag (e.g., `v1.0`):
```bash
pip install git+https://github.com/username/repository_name.git@develop
```

This method is extremely useful for sharing private tools within a team or for Continuous Integration (CI) environments.