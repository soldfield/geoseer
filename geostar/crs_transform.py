import gdal


input_epsg = 4326
output_epsg = 25832
# https://eng.sdfe.dk/product-and-services/professional-users/

input_crs = gdal.osr.SpatialReference()
input_crs.ImportFromEPSG(input_epsg)

output_crs = gdal.osr.SpatialReference()
output_crs.ImportFromEPSG(output_epsg)

coord_transform = gdal.osr.CoordinateTransformation(input_crs, output_crs)

x_lon = 9.997917
y_lat = 57.050412

mapx, mapy, z = coord_transform.TransformPoint(x_lon, y_lat)

print(mapx, mapy, z)


#https://twcc.fr/en/#
#X 560538.097
#Y 6323439.984