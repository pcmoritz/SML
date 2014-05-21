## SML: Submodular Optimization for Julia

## Installation

The easiest way to install the package is to run

    Pkg.clone("https://github.com/pcmoritz/SML")

from the Julia REPL.

## Usage

First, import the package with `import SML`. Then, you can for example
create and minimize a modular function with

    F = SML.Modular([1:10; -1])
    SML.minimize(F)

This will give `[11]` as a result. You can also call

    SML.min_norm_point(F, [1:11], 1e-12)

where the second argument is the ordering the minimum norm point will
be initialized with and `1e-12` is the duality gap to be achieved for
the minimum norm point (which roughly corresponds to the precision of
the minimum).