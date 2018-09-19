#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

from prontodb import __version__


def get_requirements():
    filepath = os.path.join(
        os.path.dirname(__file__),
        "requirements.txt"
    )

    with open(filepath) as fh:
        requirements = fh.read().splitlines()

    return requirements


setup(
    name='prontodb',
    description="Refresh Pronto with the latest data "
                "from InterPro, GOA, and UniProt",
    version=__version__,
    packages=find_packages(),
    zip_safe=False,
    install_requires=get_requirements(),
    entry_points={
        'console_scripts': [
            'pronto-update = prontodb:cli',
        ]
    }
)
