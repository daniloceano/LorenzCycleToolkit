import os

from setuptools import find_packages, setup

# Utility function to read the requirements.txt file


def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    with open(requirements_path, 'r') as f:
        return f.read().splitlines()


setup(
    name='LorenzCycleToolkit',
    version='1.0.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=read_requirements(),
    entry_points={
        'console_scripts': [
            'lorenzcycletoolkit = lorenzcycletoolkit:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)
