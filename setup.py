import setuptools
from setuptools import setup
from version import get_version

# This is actually all of them, but some give problems I don't
# want to deal with right now. Also I'm pretty sure a lot of
# basic requirements aren't even supposed to be included here

# install_requires=["pandas",
#         "casm",
#         "json",
#         "numpy",
#         "hashlib",
#         "copy",
#         "scipy",
#         "glob",
#         "pickle",
#         "mpl_toolkits",
#         "matplotlib",
#         "warnings",
#         "re",
#         "math",
#         "subprocess"]

install_requires=["pandas",
        "numpy",
        "scipy",
        "matplotlib"]


setup(name="thermoplotting",
        install_requires=install_requires,
        version=get_version(),
        description="Tools for visualizing data related to cluster expansions and manipulating casm output.",
        url='https://github.com/goirijo/thermoplotting',
        author='John Goiri',
        author_email='jg.goiri@gmail.com',
        license='MIT',
        packages=setuptools.find_packages(),
        zip_safe=False)

