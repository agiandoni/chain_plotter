# chain_plotter
Generates the required objects to quasi-1D ladder chains of certain types and the order of their MPS.

To use call generate_chain() method with the required. This returns a list of polygons and the MPS snake.
Currently, only chain I-1 of https://arxiv.org/abs/2211.07385 is supported.

Use ax.add_patch(poly) to add each element of polygons to axis ax
Use ax.plot(*mps_clean) to plot the mps
Quantities and observables can be plotted straightforwardly on individual sites if they follow the mps order.
