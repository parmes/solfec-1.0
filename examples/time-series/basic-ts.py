# basic applications of TIME_SERIES (points)
import os
d0 = os.path.dirname(os.path.realpath(__file__))

tms1 = TIME_SERIES (d0+'/ts-data-0.txt');
print 'From file:'
print tms1.times
print tms1.values

list2 = [0, 10, 1, 11, 2, 12, 3, 13, 4, 14, 5, 15, 6, 16];
tms2 = TIME_SERIES (list2);
print 'From [t0,v0, t1,v1, ..] list:'
print tms2.times
print tms2.values

list3 = [[0, 10], [1, 11], [2, 12], [3, 13], [4, 14], [5, 15], [6, 16]];
tms3 = TIME_SERIES (list3);
print 'From [[t0,v0], [t1,v1], ...] list:'
print tms3.times
print tms3.values

tms4 = tms3.derivative
print 'Derivative:'
print tms4.times
print tms4.values

tms5 = tms3.integral
print 'Integral:'
print tms5.times
print tms5.values
