#!/usr/bin/env python3
"""
微生物引物设计工具包安装脚本 / Microbial Primer Design Toolkit Setup Script
"""

from setuptools import setup, find_packages
from pathlib import Path

# 读取README文件 / Read README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

# 读取requirements文件 / Read requirements file
requirements = []
if (this_directory / "requirements.txt").exists():
    requirements = (this_directory / "requirements.txt").read_text().strip().split('\n')
    requirements = [req.strip() for req in requirements if req.strip() and not req.startswith('#')]

setup(
    name="microbial-primer-design",
    version="1.0.0",
    author="Microbial Primer Design Team",
    author_email="your.email@example.com",
    description="微生物特异性引物设计完整流水线工具包 / A comprehensive pipeline toolkit for microbial-specific primer design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/YOUR_USERNAME/microbial-primer-design",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "primer-design=run_primer_design:main",
            "download-genomes=utils.download_genomes:main",
            "parse-primer3=utils.parse_primer3:main",
            "rank-primers=utils.rank_primers:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yml", "*.txt", "*.md"],
    },
    keywords="bioinformatics primer design microbiology genomics",
    project_urls={
        "Bug Reports": "https://github.com/YOUR_USERNAME/microbial-primer-design/issues",
        "Source": "https://github.com/YOUR_USERNAME/microbial-primer-design",
        "Documentation": "https://github.com/YOUR_USERNAME/microbial-primer-design#readme",
    },
) 