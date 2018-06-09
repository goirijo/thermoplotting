# Thermoplotting

thermoplotting is my solution to the constant need for plotting thermodynamic
data in graduate school. It's meant primarily to interface with casm, but contains
a few modules that could be used with raw data, such as plotting ECI values or
energy-composition plots with convex hulls.

In addition to plotting data, new features are added as I need them, such as manipulating
basis functions, plotting primitive cells in reciprocal space, or integrating Monte Carlo
data to get Grand Canonical free energies.

## Things you need
There's a bunch of packages you'll need to get this module running. Most are available through
pip, and if you're interested in this package, you probably already have most of them:

```
numpy
pandas
scipy
matplotlib
```

If you're still using python2, you'll also need ```future```.

The one not available through pip yet is ```casm```, which you can get [here](https://github.com/prisms-center/CASMcode).

## Installation
The easiest way to install thermoplotting is to simply call
```pip install thermoplotting```

I try to keep it up to date, but I don't like uploading releases when I'm still testing things.
If you want to get the latest that's here you can do it by cloning the repository and following
the standard steps for installation via ```setup.py```, or just be sneaky like I am and softlink
the ```thermoplotting/thermoplotting``` into your ```PYTHONPATH```.
