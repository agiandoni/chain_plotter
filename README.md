# chain_plotter
Generates the required objects to draw quasi-1D ladder chains made up of hexagon/pentagon combinations, as well as the order of their MPS.

To use: call ```generate_chain(chaintype:: int, reps:: int, closinghex:: bool)```. This returns a list of polygons and the MPS snake in order.
Currently, all 6 chains of https://arxiv.org/abs/2211.07385 are supported (chaintypes are indexed 1 to 6).

Use ```ax.add_patch(poly)``` to add each element of polygons to axis ax; use ```ax.plot(*mps_clean)``` to plot the mps.

Quantities and observables can be plotted straightforwardly on individual sites if they follow the mps order.
If multi-site quantities are to be plotted, one can use geometries such as
```geometry = ["hup", "pdn", "hdn", "pup", "hup"]```
(for chain type 1)
to generate the coordinates of the lattice bonds with
```
bonds = generate_lattice_bonds(geometry)
lat = generate_lattice_from_bonds(bonds)
```
The index tuples in lat can be used in conjuction to the MPS to plot objects on the the bonds themselves.




