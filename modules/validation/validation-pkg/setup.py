from setuptools import setup, find_packages
from pathlib import Path

# Read requirements from requirements.txt
def read_requirements():
    requirements_path = Path(__file__).parent / "requirements.txt"
    with open(requirements_path) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="validation_pkg",
    version="0.1.0",
    author="Dominika Bohuslavová",
    author_email="dominikadraesslerova@gmail.com",
    description="Simple package to manage input validation for GMO pipeline",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/draessld/validation-pkg",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: EUPL-1.2 license",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12",
    install_requires=read_requirements(),
)