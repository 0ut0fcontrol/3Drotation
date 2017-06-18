import numpy as np
def ang2R(a,b,g):
    R = [[np.cos(a)*np.cos(b), -np.sin(b)*np.cos(g), -np.sin(a+g)*np.cos(b)],
         [np.cos(a)*np.sin(b),  np.cos(b)*np.cos(g), -np.sin(a+g)*np.sin(b)],
         [np.sin(a),            np.sin(g),            np.cos(a+g)],
        ]
    return np.mat(R)

print(ang2R(0.174532925199,2.79252680319,4.88692190558))
