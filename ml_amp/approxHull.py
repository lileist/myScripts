"""
DataStructure
[[R1, E1, F1],
 [R2, E2, F2],
 ...,
 [Rn, En, Fn]]
"""
import numpy as np
from scipy.spatial import ConvexHull
import copy

"""
>>> a=np.array([[0,1],[1,0],[0,0],[1,1]])
>>> hull =ConvexHull(a)
>>> hull
<scipy.spatial.qhull.ConvexHull object at 0x101bab2d0>
>>> hull.simplices
array([[2, 0],
       [2, 1],
       [3, 0],
       [3, 1]], dtype=int32)

"""
class convex_hull():
    """
    points: points that define the hull
    """
    def __init__(self, points, dataset):
        self.points = points
        self.dataset = dataset

    def get_hull(self):
        #order points to form a convex hull
        if len(self.points) < 4:
           self.hull = ConvexHull(np.array(self.dataset))
        else:   
           self.hull = ConvexHull(np.array(self.points))

    def get_facet_norm(self):
        hull_facet_norm = []
        self.get_hull()
        for facet_index in self.hull.simplices:
            cross_vect = np.cross(self.points[facet_index[1]]-self.points[facet_index[0]], \
                                  self.points[facet_index[2]]-self.points[facet_index[0]])
            cross_vect /= np.linalg.norm(cross_vect)
            hull_facet_norm.append([self.points[facet_index[0]],cross_vect])
        return np.array(hull_facet_norm)

dataset = [[1,0,0],
           [0,1,0],
           [0,0,1],
           [1,1,1],
           [0.1,0,0],
           [0.,0.1,0],
           ] 
dist_threshold = 0.2
cycle =0 
data_selected = []
data_removed = []

temp_data = copy.deepcopy(dataset)
while dataset:
  hull_points=[]
  if cycle == 0:
     min_indices = np.argmin(dataset, axis=0)
     max_indices = np.argmax(dataset, axis=0)
   
     print  min_indices
     print  max_indices
     hull_indices = list(dict.fromkeys(np.append(min_indices, max_indices)))
     for i in hull_indices:
        hull_points.append(temp_data[i])
        dataset.remove(temp_data[i])
  hull =  convex_hull(hull_points, temp_data)
  for facet_norm in hull.get_facet_norm():
     max_data = None
     distance = 0.0
     for i in range(len(dataset)):
        data = dataset[i]
        new_distance = np.dot(facet_norm[0]-data, facet_norm[1])
        #remove data point too close to convexHull
        if new_distance <= dist_threshold:
           dataset.pop(i)
           data_removed.append(data)

        if new_distance > distance:
           distance = new_distance
           max_data = data
     hull_points.append(max_data)
     dataset.remove(max_data)
  print "After", cycle, ": dataset length:", len(dataset)
  cycle += 1
  data_selected.extend(hull_points)
print "Selected:", data_selected
print "Removed:", data_removed
