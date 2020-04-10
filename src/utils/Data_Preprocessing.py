import sys
sys.path.insert(1, '/home/sawsn/Shiernee/FileIO/src/utils')
sys.path.insert(1, '/home/sawsn/Shiernee/Utils/src/utils')

from FileIO import FileIO
from Utils import Utils
import numpy as np


class Data_Preprocessing:
    def __init__(self):
        self.v_bipolar_input = None
        self.v_unipolar_output = None
        self.w_unipolar_output = None
        return

    def assign_unipolar_sqr_grid(self, V):
        """

        :param V: array (time_pt, no_pts)
        :return: array (time_pt, sqr_grid, sqr_grid)
        """
        time_pt = V.shape[0]
        sqr_grid = int(np.sqrt(V.shape[1]))
        unipolar = np.reshape(V, [time_pt, sqr_grid, sqr_grid])
        return unipolar

    def get_bipolar_sqr_grid(self, unipolar):
        """

        :param unipolar: array (time_pt, sqr_grid, sqr_grid)
        :return:  array (time_pt, sqr_grid-2, sqr_grid-2)
        """
        bipolar = np.zeros_like(unipolar)
        bipolar[:, 1:-1:10, 1:-1:10] = unipolar[:, 2::10, 2::10] - unipolar[:, :-2:2, :-2:2]

        return bipolar


class Ellipsoid:
    def __init__(self):
        return

    @staticmethod
    def surface_area(a, b, c):
        """

        :param a: float - radius of axis a
        :param b: float - radius of axis b
        :param c: float - radius of axis c
        :return: float
        """
        p = 1.6075
        tmp = ((a * b) ** p + (a * c) ** p + (b * c) ** p) / 3
        return 4 * np.pi * np.sqrt(tmp ** (1 / p))




if __name__ == '__main__':
    forward_folder = '../../data/case2_2Dgrid_55/forward_10201pts/5%sigma_gaussian_activation/'

    fileio = FileIO()
    fileio.assign_forward_folder(forward_folder)
    i = 1
    fhn_model_instances = fileio.read_physics_model_instance(i, 'fhn')
    point_cloud_instances = fileio.read_point_cloud_instance(i)

    no_pt = point_cloud_instances['no_pt']
    t = fhn_model_instances['t']
    V = fhn_model_instances['V']
    v = fhn_model_instances['v']
    coord = point_cloud_instances['coord']
    no_pt = point_cloud_instances['no_pt']
    local_axis1 = point_cloud_instances['local_axis1']
    local_axis2 = point_cloud_instances['local_axis2']


    # create square grid mesh with 0.3 unit distance (cm) - 3mm - same as Abbott inter-distance
    number_of_time_to_map = 1 #350

    map_mesh = []
    for map in range(number_of_time_to_map):
        print('map {}'.format(map))
        idx = 34  #np.random.randint(0, no_pt, 1)]
        rand_coord = coord[idx]
        local_axis1_tmp = local_axis1[idx] * 0.3
        local_axis2_tmp = local_axis2[idx] * 0.3
        print(local_axis2[idx])
        mesh = []
        print('rand_coord:{} '.format(rand_coord))
        for i in range(4):
            tmp_x = rand_coord + local_axis1_tmp * i
            for j in range(4):
                print('i {}, j {}'.format(i, j))
                tmp = tmp_x + local_axis2_tmp * j
                print(tmp)
                mesh.append(tmp)

        map_mesh.apppend(mesh)

    quit()

    time_pt = V.shape[0]
    sqr_grid = int(np.sqrt(V.shape[1]))
    unipolar = np.reshape(V, [time_pt, sqr_grid, sqr_grid])

    bipolar = unipolar[:, 2:-1, 2:-1] - unipolar[:, 1:-2, 1:-2]


