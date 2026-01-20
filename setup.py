import os

from setuptools import find_packages, setup

# Utility function to read the requirements.txt file


def read_requirements():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    with open(requirements_path, 'r') as f:
        return f.read().splitlines()


# Read the long description from README
def read_long_description():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ''


setup(
    name='LorenzCycleToolkit',
    version='1.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=read_requirements(),
    long_description=read_long_description(),
    long_description_content_type='text/markdown',
    author='Danilo Couto de Souza',
    author_email='danilo.oceano@gmail.com',
    description='A toolkit for calculating the Lorenz Energy Cycle in atmospheric systems',
    url='https://github.com/daniloceano/LorenzCycleToolkit',
    project_urls={
        'Documentation': 'https://lorenzcycletoolkit.readthedocs.io/',
        'Source': 'https://github.com/daniloceano/LorenzCycleToolkit',
        'Tracker': 'https://github.com/daniloceano/LorenzCycleToolkit/issues',
    },
    entry_points={
        'console_scripts': [
            'lorenzcycletoolkit = lorenzcycletoolkit:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)
