p2s = {}
p2s[1] = 2
p2s[4] = 3
p2s[8] = 5

nts = NEWTON_SOLVER()

hys = HYBRID_SOLVER ("parmec_file.py", 1E-5, p2s, nts)
