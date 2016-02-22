file = open("series.txt", "w")

for i in range (0, 10000000):
  j = -1 + 2*(i%2);
  file.write("%d %d\n"%(i, j));

file.close()
