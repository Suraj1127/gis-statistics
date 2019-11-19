import os

from rasterstats import zonal_stats

from gis_system import GISSystem


def generate_statistic_shapefile(geodf, statistics_type, no_of_neighbours, output_folder):

    if statistics_type == 1:
        column_name = 'surr_alt'
    elif statistics_type == 2:
        column_name = 'red_alt'
    else:
        raise Exception("Invalid statistics type entered. Only statistics 1 and 2 are accepted.")

    # save the neighbours geometry file as intermediate shape file; this will be used later for zonalstats
    geodf.drop(['geometry', 'centroid', 'centroid_point', 'neighbour_indices'], axis=1).rename(
        columns={'neighbour_geometry': 'geometry'}).to_file(
        '/tmp/geometry.shp')

    # calculate mean altitude for surrounding cells and save in another column 'surrounding_alt'
    stats_temp = zonal_stats('/tmp/geometry.shp',
                             '/home/regmi/suraj/World Bank/FB Flood Data Advanced/Data/ASTER_DEM_30M.tif')
    geodf[column_name] = [i['mean'] for i in stats_temp]

    # calculate mean altitude for the cells and save in column 'cell_alt'
    stats_temp = zonal_stats(geodf['geometry'], '/home/regmi/suraj/World Bank/FB Flood Data Advanced/Data/ASTER_DEM_30M.tif')
    geodf['cell_alt'] = [i['mean'] for i in stats_temp]

    # calculate relative height of the cells with respect to the surrounding and save as statistic_1
    geodf['stat_1'] = geodf['cell_alt'] - geodf[column_name]

    # drop redundant columns and put only the relevant columns and save as stats shapefile
    geodf.drop(
        ['centroid', 'longitude', 'latitude', 'centroid_point', 'neighbour_indices', 'neighbour_geometry'], axis=1
    ).to_file(os.path.join(output_folder, 'stat_{}_{}_cells.shp'.format(statistics_type, no_of_neighbours)))


if __name__ == "__main__":

    # create instance of GISSystem
    gis_system = GISSystem('/home/regmi/suraj/World Bank/FB Flood Data Advanced/Data/raw_pct_change_11to15_4326.shp')

    # define the variables
    no_of_neighbours = 8
    statistic_type = 2
    output_path = ''

    # by passing the maximum number of neighbour cells, return the dataframe with the neighbour_geometry column included
    geodf = gis_system.touch_and_get_dataframe(no_of_neighbours, statistic_type)

    # generate statistic shapefile and export to output_path
    generate_statistic_shapefile(geodf, statistic_type, no_of_neighbours, output_path)
    print(
        "Saved the statistic {} in the path: {}"
            .format(
            statistic_type,
            os.path.abspath('') + output_path + '/stat_' + str(statistic_type) + '_' + str(no_of_neighbours) + '_cells.shp')
    )
