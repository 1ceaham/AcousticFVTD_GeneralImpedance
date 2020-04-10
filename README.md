# AcousticFVTD_GeneralImpedance

AcousticFVTD_GeneralImpedance is an independent MATLAB implementation of the acoustic finite volume time domain (FVTD) method under general impedance boundary conditions for room acoustics as described in Bilbao et al. [1].

![Acoustic FVTD, Theater](theater4kHz.gif)

## Description

This collection of files is an initial demonstration of the cubic grid staircase approximation, which is known to be inferior for accurately generating room acoustics simulations, but is nonetheless a straightforward way to demonstrate the computational core of the FVTD algorithm using the examples given in the paper.
While this means that the scheme itself can trivially be reduced to a two-step finite difference time domain approach, it implements all of the necessary machinery to support fitted cells at boundaries or completely unstructured meshes, while hopefully preserving some degree of simplicity for understanding.
As pure MATLAB, this collection takes a non-standard approach to mesh generation and does not support common drop-in formats, though future extensions may alleviate this limitation.
Furthermore, the implementation makes tradeoffs between computational efficiency and faithfulness to the language and structure of the paper, as performance greatly influences how pleasant it is to iteratively test or simply demonstrate the most common cases.
Due to the desire for a vectorized approach, the core is more complex than a purely pedagogical representation would be, however, effort has been made to retain as much similarity to the written scheme as possible, using the same notation and avoiding excessive optimization.
Future improvements could also include a replacement with a purely loop-based or matrix-multiplication-based core to compare overall performance while improving readability.

## How to Use

To start, simply run ```template_cubic```.
Further options can be experimented with by uncommenting the first line of the script to turn it into a function or simply editing the string variables ```space``` and ```boundary```.
As always, ```help <filename>``` provides more information for any file in this collection.

## References

[1] S. Bilbao, B. Hamilton, J. Botts, and L. Savioja, “Finite Volume Time Domain Room Acoustics Simulation under General Impedance Boundary Conditions,” IEEE/ACM Transactions on Audio, Speech, and Language Processing, vol. 24, no. 1, pp. 161–173, Jan. 2016, doi:10.1109/TASLP.2015.2500018.
