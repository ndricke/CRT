import numpy as np


def rotate3Dxyz(vec, theta):
    R = np.zeros((3,3))
    R[1,1] = 0
    R[:1,:1] = [[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]]
    return np.dot(R, vec)

def rotate2D(vec, theta):
    print(theta)
    R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    print(R)
    print(vec)
    return np.dot(R, vec)


if __name__ == "__main__":
    #v1 = np.array([0.3,0,0]).T
    pca_n = 2
    vlen = 1.
    v1 = np.array([[0.],[vlen]])
    print(v1.shape)
    shift = np.zeros((pca_n, 2))
    #shift[0,:] = v1
    for i in range(pca_n):
        shift[i,:] = rotate2D(v1, i*(2*np.pi/pca_n)).T

    print(v1)
    print(shift)
