#! /usr/bin/env python
from setuptools import find_packages, setup


def read(filename):
    with open(filename, "r") as fp:
        return fp.read()


setup(
    name="cascade",
    version="0.0.1.dev0",
    description="The CoAStal Community-lAnDscape Evolution (CASCADE) model",
    long_description=open("README.md", encoding="utf-8").read(),
    author="Katherine Anarde",
    author_email="kanarde@ncsu.edu",
    url="https://github.com/UNC-CECL/cascade",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
    ],
    keywords=["coastal", "economics"],
    install_requires=open("requirements.txt", "r").read().splitlines(),
    packages=find_packages(),
    python_requires=">=3.6,<3.9",  # because of copulas
    include_package_data=True,
)