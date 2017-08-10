#!/usr/bin/env python3

from setuptools import find_packages, setup
import sys

if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

package_data = {
    "td": ["templates/*.plt",],
}

setup(
    name="td",
    version="0.1",
    description="Parser for excited state calculations",
    url="https://github.com/eljost/td",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    package_data=package_data,
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "simplejson",
        "python-docx",
        "argcomplete",
        "jinja2",
        "python-docx",
        "argcomplete",
    ],
    entry_points={
        "console_scripts": [
            "td = td.main:run",
        ]
    },
)
