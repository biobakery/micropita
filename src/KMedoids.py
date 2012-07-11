"""
k-medoids algorithm.
"""
## MLPY build 2.2.0
## This code is written by Davide Albanese, <albanese@fbk.eu>
## (C) 2009 Fondazione Bruno Kessler - Via Santa Croce 77, 38100 Trento, ITALY.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


__all__= ['Kmedoids', 'Minkowski']


import numpy as np
import mlpy


def kmedoids_core(x, med, oth, clust, cost, dist):
    """
    * for each mediod m
       * for each non-mediod data point n
         Swap m and n and compute the total cost of the configuration
    Select the configuration with the lowest cost
    """

    d = np.empty((oth.shape[0], med.shape[0]), dtype=float)
    med_n = np.empty_like(med)
    oth_n = np.empty_like(oth)
    idx = np.arange(oth.shape[0])
    
    med_cur = med.copy()
    oth_cur = oth.copy()
    clust_cur = clust.copy()
    cost_cur = cost
    
    for i, m in enumerate(med):
        for j, n in enumerate(oth[clust == i]):
           
            med_n, oth_n = med.copy(), oth.copy()

            med_n[i] = n
            tmp = oth_n[clust == i]
            tmp[j] = m
            oth_n[clust == i] = tmp
            
            for ii, nn in enumerate(oth_n):
                for jj, mm in enumerate(med_n):
                    d[ii, jj] = dist.compute(x[mm], x[nn])

            clust_n = np.argmin(d, axis=1) # clusters
            cost_n = np.sum(d[idx, clust_n]) # total cost of configuration

            if cost_n <= cost_cur:
                med_cur   = med_n.copy()
                oth_cur   = oth_n.copy()
                clust_cur = clust_n.copy()
                cost_cur  = cost_n

    return med_cur, oth_cur, clust_cur, cost_cur
            

class Kmedoids:
    """k-medoids algorithm.
    """

    def __init__(self, k, dist, maxloops=100, rs=0):
        """Initialize Kmedoids.
        
        :Parameters:
   
          k : int
              Number of clusters/medoids
          dist : class
                 class with a .compute(x, y) method which
                 returns a distance
          maxloops : int
                     maximum number of loops
          rs : int
               random seed

        Example:

        >>> import numpy as np
        >>> import mlpy
        >>> x = np.array([[ 1. ,  1.5],
        ...               [ 1.1,  1.8],
        ...               [ 2. ,  2.8],
        ...               [ 3.2,  3.1],
        ...               [ 3.4,  3.2]])
        >>> dtw = mlpy.Dtw(onlydist=True)
        >>> km = mlpy.Kmedoids(k=3, dist=dtw)
        >>> km.compute(x)
        (array([4, 0, 2]), array([3, 1]), array([0, 1]), 0.072499999999999981)

        Samples 4, 0, 2 are medoids and represent cluster 0, 1, 2 respectively.

         * cluster 0: samples 4 (medoid) and 3
         * cluster 1: samples 0 (medoid)  and 1
         * cluster 2: sample 2 (medoid)
        """
        
        self.__k = k
        self.__maxloops = maxloops
        self.__rs = rs
        self.__dist = dist

        np.random.seed(self.__rs)


    def compute(self, x):
        """Compute Kmedoids.
        
        :Parameters:
           x : ndarray
               An 2-dimensional vector (sample x features).
   
        :Returns:
           m : ndarray (1-dimensional vector)
               medoids indexes
           n : ndarray (1-dimensional vector)
               non-medoids indexes
           cl : ndarray 1-dimensional vector)
                cluster membership for non-medoids.
                Groups are in 0, ..., k-1
           co : double
                total cost of configuration
        """

        # randomly select k of the n data points as the mediods
        idx = np.arange(x.shape[0])
        np.random.shuffle(idx)
        med = idx[0:self.__k]
        oth = idx[self.__k::]

        # compute distances
        d = np.empty((oth.shape[0], med.shape[0]), dtype=float)
        for i, n in enumerate(oth):
            for j, m in enumerate(med):
                d[i, j] = self.__dist.compute(x[m], x[n])

        # associate each data point to the closest medoid
        clust = np.argmin(d, axis=1)

        # total cost of configuration
        cost = np.sum(d[np.arange(d.shape[0]), clust])

        # repeat kmedoids_core until there is no change in the medoid
        for l in range(self.__maxloops):
          
            med_n, oth_n, clust_n, cost_n = kmedoids_core(x, med, oth, clust, cost, self.__dist)
                      
            if (cost_n < cost):
                med, oth, clust, cost = med_n, oth_n, clust_n, cost_n
            else:
                break
            
        return med, oth, clust, cost


class Minkowski:
    """
    Computes the Minkowski distance between two vectors ``x`` and ``y``.

    .. math::

      {||x-y||}_p = (\sum{|x_i - y_i|^p})^{1/p}.
    """

    def __init__(self, p):
        """
        Initialize Minkowski class.
        
        :Parameters:
          p : float
              The norm of the difference :math:`{||x-y||}_p`
        """

        self.__p = p
                

    def compute(self, x, y):
        """
        Compute Minkowski distance
        
        :Parameters:
           x : ndarray
               An 1-dimensional vector.
           y : ndarray
               An 1-dimensional vector.

        :Returns:
           d : float
               The Minkowski distance between vectors ``x`` and ``y``           
        """
        
        return (abs(x - y)**self.__p).sum() ** (1.0 / self.__p)

    
