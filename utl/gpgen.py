from math import sqrt

x0 = 1.0 / sqrt (3.0)
x1 = sqrt (3.0 / 5.0)
x2 = sqrt ((3.0 - 2.0 * sqrt (6.0 / 5.0)) / 7.0)
x3 = sqrt ((3.0 + 2.0 * sqrt (6.0 / 5.0)) / 7.0)
x4 = (1.0 / 3.0) * sqrt (5.0 - 2.0 * sqrt (10.0 / 7.0))
x5 = (1.0 / 3.0) * sqrt (5.0 + 2.0 * sqrt (10.0 / 7.0))

w2 = (18.0 + sqrt (30.0)) / 36.0
w3 = (18.0 - sqrt (30.0)) / 36.0
w4 = (322.0 + 13.0 * sqrt (70.0)) / 900.0
w5 = (322.0 - 13.0 * sqrt (70.0)) / 900.0

onedim = [(1, (0.0), (2.0)),
          (2, (-x0, x0), (1.0, 1.0)),
	  (3, (-x1, 0.0, x1), (5.0/9.0, 8.0/9.0, 5.0/9.0)),
	  (4, (-x3, -x2, x2, x3), (w3, w2, w2, w3)),
	  (5, (-x5, -x4, 0.0, x4, x5), (w5, w4, 128.0 / 225.0, w4, w5))]


order = 4

data = onedim [order-1]
zipped = zip (data [1], data [2])

X_out = []
Y_out = []
Z_out = []
W_out = []

for x1, w1 in zipped:
  for x2, w2 in zipped:
    for x3, w3 in zipped:
      X_out.append (x1)
      Y_out.append (x2)
      Z_out.append (x3)
      W_out.append (w1*w2*w3)

print 'Order', order, ' len = ', len (X_out)
print 'X = ', X_out
print 'Y = ', Y_out
print 'Z = ', Z_out
print 'W = ', W_out
