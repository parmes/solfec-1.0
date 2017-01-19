matnum = MATERIAL (1, 1, 0.25)

nodes = [0, 0, 0,
         1, 0, 0,
	 1, 1, 0,
	 0, 1, 0,
	 0, 0, 1,
	 1, 0, 1,
	 1, 1, 1,
	 0, 1, 1]

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

parnum = MESH (nodes, elements, matnum, colors)

CONSTRAIN (parnum, [1, 0, 0, 0, 1, 0], [1, 0, 0, 0, 1, 0, 0, 0, 1])

SPRING (parnum, (0.5, 0.5, 0.0), -1, (0.5, 0.5, -1.0), [-1, -100, 1, 100], [-1, -1, 1, 1], (0, 0, 1))

GRAVITY (0, 0, -10)

#DEM (10, 0.01)
