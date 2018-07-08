import sys
import numpy as np
import math

class quaternion(object):

    def __init__(self, px,py,pz,ps):

        c = math.sqrt(px*px + py*py + pz*pz + ps*ps)
        qx = px/c
        qy = py/c
        qz = pz/c
        qs = ps/c

        self.rot11 = qs*qs + qx*qx - qy*qy - qz*qz
        self.rot12 = 2.0 * (qx*qy - qs*qz)
        self.rot13 = 2.0 * (qx*qz + qs*qy)
        self.rot21 = 2.0 * (qx*qy + qs*qz)
        self.rot22 = qs*qs - qx*qx + qy*qy - qz*qz
        self.rot23 = 2.0 * (qy*qz - qs*qx)
        self.rot31 = 2.0 * (qx*qz - qs*qy)
        self.rot32 = 2.0 * (qy*qz + qs*qx)
        self.rot33 = qs*qs - qx*qx - qy*qy + qz*qz

        self.x = qx
        self.y = qy
        self.z = qz
        self.s = qs

    def rotate(self, x,y,z):

        nx = self.rot11 * x + self.rot12 * y + self.rot13 * z
        ny = self.rot21 * x + self.rot22 * y + self.rot23 * z
        nz = self.rot31 * x + self.rot32 * y + self.rot33 * z
    
        return nx, ny, nz

    def __str__(self):
        out_str = ''
        out_str += str(self.rot11) + ' ' + str(self.rot12) + ' ' + str(self.rot13) + '\n'
        out_str += str(self.rot21) + ' ' + str(self.rot22) + ' ' + str(self.rot23) + '\n'
        out_str += str(self.rot31) + ' ' + str(self.rot32) + ' ' + str(self.rot33) + '\n'
        return out_str


def main(argv):

    if len(argv)!=2:
        print 'Usage: ' + argv[0] + ' <num samples>'
        exit(1);

    num_samples = int(argv[1])

    # Variance in the directions of 3 coordinate axes
    sigmas = np.random.uniform(0.0, 1.0, 3)

    # Generate points in gaussian distribution with 0 mean and given variance
    x_dist = np.random.normal(0.0, sigmas[0], num_samples)
    y_dist = np.random.normal(0.0, sigmas[1], num_samples)
    z_dist = np.random.normal(0.0, sigmas[2], num_samples)

    # Generate random orientation as a quaternion
    q_param = np.random.uniform(-1.0, 1.0, 4)
    quat = quaternion(q_param[0],q_param[1],q_param[2],q_param[3])


    # Generate randam position as mean.
    #    mean = np.random.uniform(-1.0, 1.0, 3)
    mean = [0.0, 0.0, 0.0]

    print ''
    print '# Mean: (' + str(mean[0]) + ', ' + str(mean[1]) + ', ' + str(mean[2]) + ')'
    print ''
    ax = quat.rotate(1.0, 0.0, 0.0)
    print '# X:    (' + str(ax[0]) + ', ' + str(ax[1]) + ', ' + str(ax[2]) + ')'
    ay = quat.rotate(0.0, 1.0, 0.0)
    print '# Y:    (' + str(ay[0]) + ', ' + str(ay[1]) + ', ' + str(ay[2]) + ')'
    az = quat.rotate(0.0, 0.0, 1.0)
    print '# Z:    (' + str(az[0]) + ', ' + str(az[1]) + ', ' + str(az[2]) + ')'
    print ''

    for p in zip(x_dist, y_dist, z_dist):
        # Orient -> Translate to mean.
        rotated_p = quat.rotate(p[0],p[1],p[2])
        outstr = ''
        outstr += str(mean[0] + rotated_p[0]) + ', '
        outstr += str(mean[1] + rotated_p[1]) + ', '
        outstr += str(mean[2] + rotated_p[2])
        print outstr

if __name__ == "__main__":
    main(sys.argv)
