# coding: utf-8
from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

from figaro import __version__ as version

here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    description = long_description = description.read()

    name = "figaro"
    version = version

    if sys.version_info.major != 3:
        raise EnvironmentError(
            """{toolname} is a python module that requires python3, and is not compatible with python2.""".format(
                toolname=name
            )
        )

    setup(
        name=name,
        version=version,
        description="FIGARO - An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters",
        long_description=long_description,
        license="GPLv3",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Topic :: Scientific Engineering :: Bio/Informatics",
            "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 3.7",
        ],
        zip_safe=False,
        keywords="",
        packages=find_packages(exclude=["test"]),
        install_requires=list(req.strip() for req in open("requirements.txt")),
        entry_points={
            "console_scripts": ["figaro=figaro.figaro:main"],
        },
        scripts=[],
        package_data={},
        include_package_data=True,
        data_files=[],
    )
