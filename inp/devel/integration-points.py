a=1.0
b=1.0
c=1.0

print "%g x %g x %g cube" % (a, b, c)

Nodes = [0, 0, 0, a, 0, 0, a, b, 0, 0, b, 0, 0, 0, c, a, 0, c, a, b, c, 0, b, c]

Surf = [1, 2, 3, 4, 5, 6]

Mesh = HEX(Nodes, 1, 1, 1, 1, Surf)

print 'All surface integration points:'
print Mesh.surface_integration_points()
for i in Surf:
  print 'Surface %d points:' % i
  print Mesh.surface_integration_points(i)

print 'All volume integration points:'
print Mesh.volume_integration_points()
