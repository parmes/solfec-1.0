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

color = 1

parnum = MESH (nodes, elements, matnum, color)

RESTRAIN (parnum, [1, 0, 0, 0, 1, 0], [1, 0, 0, 0, 1, 0, 0, 0, 1])

SPRING (parnum, (0.5, 0.5, 0.0), -1, (0.5, 0.5, -1.0),
        spring = [-1, -100, 1, 100], direction = (0, 0, 1))

GRAVITY (0.0, 0.0, -10.0)

#DEM (10, 0.01)
