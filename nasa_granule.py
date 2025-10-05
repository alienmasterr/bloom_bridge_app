from qgis.core import QgsRasterLayer, QgsProject
import numpy as np
from osgeo import gdal
import os

# Функція для обчислення геотрансформації для однієї координати
def calculate_geotransform(lat, lon, width, height, pixel_size=0.00054):
    ulx = lon - (width * pixel_size / 2)
    uly = lat + (height * pixel_size / 2)
    return (ulx, pixel_size, 0, uly, 0, -pixel_size)

# Функція для витягнення унікальних координат з CSV шару
def get_unique_coords(csv_layer_name):
    csv_layer = None
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == csv_layer_name:
            csv_layer = layer
            break
    if not csv_layer:
        print(f"Error: CSV layer {csv_layer_name} not found!")
        return []
    
    unique_coords = set()
    for feature in csv_layer.getFeatures():
        lat = feature['decimalLatitude']
        lon = feature['decimalLongitude']
        if lat is not None and lon is not None:
            unique_coords.add((lat, lon))
    
    if not unique_coords:
        print("Error: No valid coordinates in CSV!")
    return list(unique_coords)

# Функція для геореференсування TIFF файлу для кожної координати
def georeference_tiff(tiff_path, output_dir, geotransform, coord_idx):
    if not os.path.exists(tiff_path):
        print(f"Error: TIFF file {tiff_path} not found!")
        return
    
    gdal_ds = gdal.Open(tiff_path)
    if gdal_ds is None:
        print(f"Error: Failed to open TIFF file {tiff_path}")
        return
    
    output_path = os.path.join(output_dir, f"georeferenced_tiff_coord_{coord_idx}.tif")
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.CreateCopy(output_path, gdal_ds, 0)
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
    out_ds = None
    gdal_ds = None
    
    layer_name = f"Georeferenced_TIFF_coord_{coord_idx}"
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == layer_name:
            QgsProject.instance().removeMapLayer(layer.id())
    
    new_layer = QgsRasterLayer(output_path, layer_name)
    if new_layer.isValid():
        QgsProject.instance().addMapLayer(new_layer)
        print(f"{layer_name} created/updated!")
    else:
        print(f"Error creating {layer_name}")

# Основна логіка
def georeference_tiff_from_csv():
    csv_layer_name = "Zinnia_acerosa_DC_A_Gray_OpenFlowers"  # Заміни на назву твоєї CSV шару
    tiff_path = "/Users/alina/Downloads/EMIT_L2A_RFL_001_20230224T181000_2305512_020.tiff"
    output_dir = "/Users/alina/Documents/Alina/education/hakatons/georeferenced_tiff"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    unique_coords = get_unique_coords(csv_layer_name)
    if not unique_coords:
        return
    
    gdal_ds = gdal.Open(tiff_path)
    if gdal_ds is None:
        print(f"Error: Failed to open TIFF file {tiff_path}")
        return
    height = gdal_ds.RasterYSize
    width = gdal_ds.RasterXSize
    gdal_ds = None
    
    for idx, (lat, lon) in enumerate(unique_coords):
        geotransform = calculate_geotransform(lat, lon, width, height)
        georeference_tiff(tiff_path, output_dir, geotransform, idx)
    
    print("Georeferencing completed!")

georeference_tiff_from_csv()