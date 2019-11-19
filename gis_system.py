import copy

import pandas as pd
import geopandas as gp

from shapely.ops import cascaded_union


class GISSystem:

    def __init__(self, gis_file_path):
        """
        Create instance using shape file of flood affected areas (gis_file_path)
        """
        self.geodf = gp.read_file(gis_file_path)
        self.neighbours_geometry_map = {}

        # create intermediate variables
        self._create_intermediate_variables()

    def _get_category(self, flow):
        """
        Categorize how the cells are affected.
        GREEN: Equal to or more than 30% coming in -- safe areas.
        RED: Equal to or more than 30% going out -- danger areas.
        YELLOW: Less than 30% either coming in or going out -- unaffected areas.
        """
        if flow >= 30:
            return 'GREEN'
        elif flow <= -30:
            return 'RED'
        else:
            return 'YELLOW'

    def _create_intermediate_variables(self):
        # calculate centroid of the geometries to make calculations/processing easier
        self.geodf['centroid'] = self.geodf['geometry'].centroid

        # break down centroid into longitude and latitude
        self.geodf['longitude'] = self.geodf['centroid'].apply(lambda point: point.x)
        self.geodf['latitude'] = self.geodf['centroid'].apply(lambda point: point.y)

        # reset the indices there and put default indices provided by pandas
        self.geodf = self.geodf.reset_index(drop=True)

        # seeing if column "VALUE" represents percentage
        print("Seeing if column 'VALUE' represents percentage...")
        print(self.geodf['VALUE'].max(), self.geodf['VALUE'].min(), self.geodf['VALUE'].mean())

        # categorize the cells
        self.geodf['celltype'] = self.geodf['VALUE'].apply(self._get_category)

        # create sorted lists of longitude and latitude
        self.long = sorted(set(self.geodf['longitude']))
        self.lat = sorted(set(self.geodf['latitude']))

        # map longitudes and latitudes to indices in the sorted lists
        self.long_map = dict((j, i) for i, j in enumerate(self.long))
        self.lat_map = dict((j, i) for i, j in enumerate(self.lat))

        # create new column centroid_point based on longitudes and latitudes
        self.geodf['centroid_point'] = list(zip(self.geodf['latitude'], self.geodf['longitude']))

        # create map from centroid_point to their indices in the dataframe and vice versa
        self.cp2ind = dict(zip(self.geodf['centroid_point'], self.geodf.index))
        self.ind2cp = dict((j, i) for i, j in self.cp2ind.items())

        # create dataframe with sorted latitudes and longitudes as indices and columns
        # respectively, and geodf indices as the value
        self.df = pd.DataFrame([[self.cp2ind.get((j, i)) for i in self.long] for j in self.lat])

    def _get_neighbours_indices(self, j, n, statistic_type):
        """
        Gets the neighbours indices for each row in geodataframe, geodf.
        24 surrounding cells at maximum are considered to generate neighbouring indices.
        """
        try:
            neighbours = set(
                self.df.iloc[
                    self.lat_map[j[0]] + k, self.long_map[j[1]] + l
                ]
                for k in range(-n, n+1) for l in range(-n, n+1) if ((l != 0) or (k != 0))
            )
        except:
            neighbours = set()

            for k in range(-n, n+1):
                for l in range(-n, n+1):

                    if k == l == 0:
                        continue

                    try:
                        if statistic_type == 1:
                            neighbours.add(self.df.iloc[self.lat_map[j[0]] + k, self.long_map[j[1]] + l])
                        else:
                            if self.geodf.loc[
                                self.df.iloc[self.lat_map[j[0]] + k, self.long_map[j[1]] + l], 'celltype'
                            ] == 'RED':
                                neighbours.add(self.df.iloc[self.lat_map[j[0]] + k, self.long_map[j[1]] + l])
                    except:
                        pass

        return neighbours

    def _get_neighbour_geometry(self, neighbour_set):

        if neighbour_set == set():
            return None

        neighbour_set = copy.deepcopy(neighbour_set)

        if frozenset(neighbour_set) in self.neighbours_geometry_map:
            return self.neighbours_geometry_map[frozenset(neighbour_set)]

        if len(neighbour_set) == 1:
            return self.geodf.iloc[int(neighbour_set.pop()), 1]

        popped_index = int(neighbour_set.pop())

        unioned_geometry = cascaded_union(
            [self._get_neighbour_geometry(neighbour_set), self.geodf.iloc[popped_index, 1]]
        )

        neighbour_set.add(popped_index)
        self.neighbours_geometry_map[frozenset(neighbour_set)] = unioned_geometry

        return unioned_geometry

    def append_neighbour_geometry(self, no_of_neighbours, statistic_type):

        # find the depth of neighbourhood
        n = int(((no_of_neighbours + 1) ** 0.5 - 1) / 2)

        # add neighbour_indices column for neighbouring indices
        self.geodf['neighbour_indices'] = self.geodf['centroid_point'].apply(self._get_neighbours_indices, args=(n, statistic_type))

        # clean the neighbour_indices column by getting rid of nans
        self.geodf['neighbour_indices'] = self.geodf['neighbour_indices'].apply(
            lambda x: set(i for i in x if str(i) != 'nan')
        )

        # union all the neighbours and produce another column containing geometry of unioned neighbours
        # use memoization to speed up the process
        neighbour_geometry = self.geodf['neighbour_indices'].apply(self._get_neighbour_geometry)

        self.geodf['neighbour_geometry'] = neighbour_geometry

        # remove the null neighbour_geometry values from the GeoDataFrame
        self.geodf = self.geodf[~self.geodf['neighbour_geometry'].isnull()]

    def touch_and_get_dataframe(self, no_of_neighbours, statistic_type):
        """
        Returns the dataframe with new column 'neighbour_geometry' appended. The new column 'neighbour_geometry' has
        unioned geometry of all the neighbour cells.
        """

        if int(((no_of_neighbours + 1) ** 0.5 - 1) / 2) != ((no_of_neighbours + 1) ** 0.5 - 1) / 2:
            raise Exception("Invalid number of neighbours entered.")

        self.append_neighbour_geometry(no_of_neighbours, statistic_type)

        return self.geodf
