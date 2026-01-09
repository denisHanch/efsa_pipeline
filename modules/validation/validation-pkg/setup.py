from setuptools import setup, find_packages
from pathlib import Path

# Get the directory containing setup.py
HERE = Path(__file__).parent

# Read requirements from requirements.txt
def read_requirements():
    requirements_path = HERE / "requirements.txt"
    with open(requirements_path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]

# Read README from parent directory
readme_path = HERE.parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="validation_pkg",
    version="0.1.0",
    author="Dominika Bohuslavová",
    author_email="dominikadraesslerova@gmail.com",
    description="Simple package to manage input validation for GMO pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: EUPL-1.2 license",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest==8.4.2",
        ],
    },
)