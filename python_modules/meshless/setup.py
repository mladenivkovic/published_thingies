from setuptools import setup
#  from distutils.core import setup

setup(
        name='Meshless Effective Surfaces',
        version='0.1.0',
        author='Mladen Ivkovic',
        author_email='mladen.ivkovic@hotmail.com',
        packages=['meshless'],
        license='GLPv3',
        scripts=[
                    "examples/example_ivanova.py",
                    "examples/example_ivanova_integrand.py",
                    "examples/example_hopkins.py",
                    "examples/check_volume.py",
                    "examples/check_versions.py",
                    "examples/check_kernels.py",
                    "examples/check_directions.py"
                ],
        long_description=open('README.md').read(),
        install_requires=[
            'numpy',
            'matplotlib'
        ]
     
     
     )
