from setuptools import setup, find_packages
import pkg_resources

def parse_requirements(filename):
    with open(filename, 'r') as file:
        return [str(req) for req in pkg_resources.parse_requirements(file)]

setup(
    name='LorenzCycleToolkit',
    version='0.1.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=parse_requirements('requirements.txt'),
    entry_points={
        'console_scripts': [
            'lorenzcycletoolkit=lorenzcycletoolkit.main:main',  # Adjust this as necessary
        ],
    },
    author='Danilo Couto de Souza',
    author_email='danilo.oceano@example.com',
    description='A toolkit for computing the Lorenz Energy Cycle in atmospheric regions.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/daniloceano/LorenzCycleToolkit',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
