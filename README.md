# ring
Python code for Blender to extrude a polygon with a slope, and a bit of code for generating a ring.

This project is about designing a ring which could be printed on a 3D printer.
The goal is to add various drawings to the surface of the ring in a way they can be filled with a different material,
potentially the inverse of the inscription could be printed and welded in.

To achieve that goal it was necessary to extrude a polygon with a given slope. I didn't find any utility that could do that, 
so I wrote a python code which does that by triangle based mesh generation. That's carve_extrude.

hullam.py is a code for lifting a wave of the basic ring shape.

The PNG files are samples of the results.
