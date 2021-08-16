from setuptools import find_packages,setup

requirements = open("requirements.txt").read().strip().split("\n")

exec(open('hotSpot/__init__.py').read())

setup(
    name='hotSpot',
    packages=find_packages(),
    entry_points={
        "console_scripts": ['hotSpot = hotSpot.run:main']
        },
    version=__version__,
    description="Script to find immunogenic hotspots",
    author="Manju Kashyap",
    license="MIT",
    install_requires=requirements,
    url='https://github.com/manju-uss'
) 
