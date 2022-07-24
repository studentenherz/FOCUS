G-EQDSK file format
===================

*G-EQDSK* (G-formatted EQuilibrium DiSK) is a file format used to store the data of an equilibrium from *EFIT*'s reconstruction. This format is widely used by experimentals to share the results of discharges made in tokamaks. The code implemented is based on information scattered all around the internet [1_, 2_, 3_].

*G-EQDSK* files don't contain information about the coordinate convention they are using, so people should know beforehand from the file's source. The COCOS_ convention defines an index defining the coordinate systems chosen, choice of sign and normalization.

.. _1: https://fusion.gat.com/conferences/snowmass/working/mfe/physics/p3/equilibria/g_eqdsk_s.pdf
.. _2: https://github.com/bendudson/freegs/blob/master/freegs/_geqdsk.py
.. _3: https://github.com/bendudson/pyTokamak/blob/master/tokamak/formats/geqdsk.py
.. _COCOS: https://www.sciencedirect.com/science/article/abs/pii/S0010465512002962


These files can be read with the following function

.. doxygenfunction:: read_geqdsk

