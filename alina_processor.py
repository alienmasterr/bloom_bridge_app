from qgis.core import QgsRasterLayer, QgsProject, QgsVectorLayer, QgsSingleBandPseudoColorRenderer, QgsColorRampShader, QgsRasterShader, QgsRasterBandStats
from qgis.PyQt.QtCore import QDate
import netCDF4 as nc
import numpy as np
from osgeo import gdal
import os
import time
from qgis.PyQt.QtGui import QColor

# Налаштування
project_dir = "/Users/alina/Downloads"
output_base_dir = "/Users/alina/Documents/Alina/education/hakatons"

# Словник пар NetCDF–CSV
nc_csv_pairs = {
    "/Users/alina/Downloads/EMIT_L2A_RFL_001_20230107T021500_2300701_002.nc": "/Users/alina/Downloads/Zinnia_acerosa_DC_A_Gray_OpenFlowers.csv",
}

# TIFF файл для геореференсування
tiff_path = "/Users/alina/Downloads/EMIT_L2A_RFL_001_20230224T181000_2305512_020.tiff"

# Словник індексів
all_indices = {
    'NDVI': {'wl': [850, 650], 'formula': '(nir - red) / (nir + red)'},
    'PRI': {'wl': [531, 570], 'formula': '(r531 - r570) / (r531 + r570)'},
    'ARI': {'wl': [550, 700], 'formula': '1 / r550 - 1 / r700'},
    'CRI': {'wl': [510, 550], 'formula': '1 / r510 - 1 / r550'},
    'EVI': {'wl': [850, 650, 450], 'formula': '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)'},
    'FLI': {'wl': [530, 670, 700, 550], 'formula': '((r530 - r670) / (r530 + r670)) * (r700 / r550)'},
    'REIP': {'wl': [680, 750], 'formula': 'spectral_derivative_min(680-750)'},
    'mNDVI': {'wl': [750, 705], 'formula': '(r750 - r705) / (r750 + r705)'},
    'mSR': {'wl': [750, 705], 'formula': 'r750 / r705'},
    'Flower_Fraction': {'wl': [500, 700], 'formula': '(r500 - r_background) / (r_leaf - r_background)'}
}

# Функція для пошуку найближчої band
def find_band(wavelength, wavelengths):
    idx = np.argmin(np.abs(wavelengths - wavelength))
    if idx >= len(wavelengths):
        print(f"Error: Wavelength {wavelength} not found in dataset!")
        return None
    return idx + 1

# Функція для обчислення геотрансформації
def calculate_geotransform(lat, lon, width, height):
    pixel_size = 0.00054
    ulx = lon - (width * pixel_size / 2)
    uly = lat + (height * pixel_size / 2)
    print(f"Calculated center: Lat {lat:.6f}, Lon {lon:.6f}, Width: {width}, Height: {height}")
    return (ulx, pixel_size, 0, uly, 0, -pixel_size)

# Функція для вибору Delimited Text Layer
def get_delimited_text_layer(csv_path):
    layer_name = os.path.splitext(os.path.basename(csv_path))[0]
    for layer in QgsProject.instance().mapLayers().values():
        if isinstance(layer, QgsVectorLayer) and layer.name() == layer_name:
            return layer
    uri = f"file://{csv_path}?delimiter=,&yField=decimalLatitude&xField=decimalLongitude&crs=EPSG:4326"
    layer = QgsVectorLayer(uri, layer_name, "delimitedtext")
    if layer.isValid():
        QgsProject.instance().addMapLayer(layer)
        print(f"Loaded CSV layer: {layer_name}")
        return layer
    print(f"Error: Failed to load CSV layer from {csv_path}")
    return None

# Функція для витягнення дат і унікальних локацій
def extract_dates_locations(layer, nc_date):
    dates = []
    locations = set()
    for feature in layer.getFeatures():
        date_field = 'date_x'
        lat_field, lon_field = 'decimalLatitude', 'decimalLongitude'
        if date_field in feature.fields().names():
            date_value = feature[date_field]
            if date_value:  # Перевіряємо, що дата не None
                dates.append(date_value)
        if lat_field in feature.fields().names() and lon_field in feature.fields().names():
            lat, lon = feature[lat_field], feature[lon_field]
            if lat is not None and lon is not None:
                locations.add((lat, lon))
                print(f"Added location: Lat {lat:.6f}, Lon {lon:.6f}")
    date_from = min(dates) if dates else None
    date_to = max(dates) if dates else None
    print(f"Date range: {date_from} to {date_to}, Found {len(locations)} locations")
    return date_from, date_to, list(locations)

# Функція для стилізації растрового шару
def style_raster_layer(layer, index_name):
    if not layer.isValid():
        print(f"Error: Layer {layer.name()} is not valid")
        return
    shader = QgsRasterShader()
    color_ramp = QgsColorRampShader()
    
    # Словник унікальних кольорів для кожного індексу
    color_schemes = {
        'NDVI': (QColor(0, 255, 0), QColor(0, 100, 0)),  # Зелений → темно-зелений
        'PRI': (QColor(0, 0, 255), QColor(0, 255, 255)),  # Синій → блакитний
        'ARI': (QColor(255, 165, 0), QColor(255, 69, 0)),  # Помаранчевий → темно-помаранчевий
        'CRI': (QColor(128, 0, 128), QColor(255, 0, 255)),  # Фіолетовий → світло-фіолетовий
        'EVI': (QColor(0, 128, 0), QColor(0, 255, 128)),  # Темно-зелений → світло-зелений
        'FLI': (QColor(255, 255, 0), QColor(255, 165, 0)),  # Жовтий → помаранчевий
        'REIP': (QColor(255, 0, 0), QColor(139, 0, 0)),  # Червоний → темно-червоний
        'mNDVI': (QColor(0, 255, 0), QColor(0, 100, 0)),  # Зелений → темно-зелений (як NDVI)
        'mSR': (QColor(0, 255, 255), QColor(0, 128, 128)),  # Блакитний → темно-блакитний
        'Flower_Fraction': (QColor(255, 192, 203), QColor(255, 20, 147)),  # Рожевий → темно-рожевий
        'Unmixing': (QColor(255, 255, 0), QColor(255, 0, 0)),  # Жовтий → червоний
        'TIFF': (QColor(128, 128, 128), QColor(255, 255, 255))  # Сірий → білий
    }
    
    # Обчислення реального діапазону значень
    stats = layer.dataProvider().bandStatistics(1, QgsRasterBandStats.All)
    min_value = stats.minimumValue
    max_value = stats.maximumValue
    if index_name in ['NDVI', 'mNDVI']:
        min_value = max(min_value, -1)  # Обмежуємо NDVI/mNDVI до -1
        max_value = min(max_value, 1)   # Обмежуємо NDVI/mNDVI до 1
    
    # Вибір кольорів для індексу
    color_min, color_max = color_schemes.get(index_name, color_schemes['Unmixing'])
    items = [
        QgsColorRampShader.ColorRampItem(min_value, color_min),
        QgsColorRampShader.ColorRampItem(max_value, color_max)
    ]
    color_ramp.setColorRampType(QgsColorRampShader.Interpolated)
    color_ramp.setColorRampItemList(items)
    shader.setRasterShaderFunction(color_ramp)
    renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)
    layer.setRenderer(renderer)
    layer.setOpacity(0.7)  # Додаємо прозорість 70%
    layer.triggerRepaint()
    print(f"Applied {index_name} color ramp to {layer.name()}: min {color_min.name()}, max {color_max.name()}")

# Функція для геореференсування TIFF
def georeference_tiff(tiff_path, output_dir, geotransform, coord_idx, csv_name):
    start_time = time.time()
    print(f"=== Georeferencing TIFF for coord {coord_idx} ===")
    output_path = os.path.join(output_dir, f"{csv_name}_Georeferenced_TIFF_coord_{coord_idx}.tif")
    ds = gdal.Open(tiff_path)
    if ds is None:
        print(f"Error: Failed to open TIFF file {tiff_path}")
        return
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.CreateCopy(output_path, ds, 0)
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
    out_ds = None
    ds = None
    layer_name = f"{csv_name}_Georeferenced_TIFF_coord_{coord_idx}"
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == layer_name:
            QgsProject.instance().removeMapLayer(layer.id())
    new_layer = QgsRasterLayer(output_path, layer_name)
    if new_layer.isValid():
        QgsProject.instance().addMapLayer(new_layer)
        style_raster_layer(new_layer, "TIFF")
        print(f"Georeferenced TIFF layer coord {coord_idx} created/updated! Took {time.time() - start_time:.2f} seconds")
    else:
        print(f"Error creating georeferenced TIFF layer coord {coord_idx}")

# Функція для створення RGB-шару
def add_netcdf_rgb_layer(reflectance_data, wavelengths, output_dir, geotransform, coord_idx, csv_name):
    start_time = time.time()
    print(f"=== Adding NetCDF RGB Layer for coord {coord_idx} ===")
    rgb_bands = [650, 550, 450]
    band_indices = [find_band(wl, wavelengths) for wl in rgb_bands]
    if None in band_indices:
        print(f"Error: One or more RGB bands ({rgb_bands}) not found in dataset!")
        return
    rgb_data = reflectance_data[:, :, [idx-1 for idx in band_indices]]
    rgb_data = np.rot90(rgb_data, k=1)
    if np.any(np.isnan(rgb_data)) or np.any(np.isinf(rgb_data)):
        print("Error: RGB data contains NaN or infinite values")
        return
    rgb_data = np.clip(rgb_data, 0, 1)
    rgb_data = (rgb_data * 255).astype(np.uint8)
    output_path = os.path.join(output_dir, f"{csv_name}_NetCDF_RGB_coord_{coord_idx}.tif")
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_path, rgb_data.shape[1], rgb_data.shape[0], 3, gdal.GDT_Byte)
    for i in range(3):
        out_ds.GetRasterBand(i + 1).WriteArray(rgb_data[:, :, i])
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
    out_ds = None
    layer_name = f"{csv_name}_NetCDF_RGB_Layer_coord_{coord_idx}"
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == layer_name:
            QgsProject.instance().removeMapLayer(layer.id())
    new_layer = QgsRasterLayer(output_path, layer_name)
    if new_layer.isValid():
        QgsProject.instance().addMapLayer(new_layer)
        print(f"NetCDF RGB layer coord {coord_idx} created/updated! Took {time.time() - start_time:.2f} seconds")
    else:
        print(f"Error creating NetCDF RGB layer coord {coord_idx}")

# Функція для розрахунку індексу
def calculate_index(index_name, wl_list, formula, reflectance_data, wavelengths, output_dir, geotransform, coord_idx, csv_name):
    start_time = time.time()
    bands = [find_band(wl, wavelengths) for wl in wl_list]
    if None in bands:
        print(f"Error: Invalid bands for {index_name}")
        return
    data = np.rot90(reflectance_data, k=1)
    output_path = os.path.join(output_dir, f"{csv_name}_{index_name}_coord_{coord_idx}.tif")
    if index_name == 'REIP':
        red_edge_data = data[..., bands[0]-1:bands[1]]
        deriv = np.gradient(red_edge_data, axis=2)
        reip_values = wavelengths[np.argmax(deriv, axis=2)]
    else:
        var_map = {
            'nir': data[..., bands[0]-1] if 'nir' in formula else None,
            'red': data[..., bands[1]-1] if 'red' in formula else None,
            'blue': data[..., bands[2]-1] if len(bands) > 2 and 'blue' in formula else None,
            'r531': data[..., bands[0]-1] if 'r531' in formula else None,
            'r570': data[..., bands[1]-1] if 'r570' in formula else None,
            'r550': data[..., bands[0]-1] if 'r550' in formula else (data[..., bands[3]-1] if len(bands) > 3 else None),
            'r700': data[..., bands[1]-1] if 'r700' in formula else (data[..., bands[2]-1] if len(bands) > 2 else None),
            'r510': data[..., bands[0]-1] if 'r510' in formula else None,
            'r530': data[..., bands[0]-1] if 'r530' in formula else None,
            'r670': data[..., bands[1]-1] if 'r670' in formula else None,
            'r750': data[..., bands[0]-1] if 'r750' in formula else None,
            'r705': data[..., bands[1]-1] if 'r705' in formula else None,
            'r500': data[..., bands[0]-1] if 'r500' in formula else None,
            'r_background': 0.05,
            'r_leaf': 0.1
        }
        calc_formula = formula
        for var in var_map:
            if var_map[var] is None and var in calc_formula:
                print(f"Error: Variable {var} not defined for {index_name}")
                return
        locals_dict = {k: v for k, v in var_map.items() if v is not None}
        try:
            result = eval(calc_formula, {}, locals_dict)
        except Exception as e:
            print(f"Error calculating {index_name}: {e}")
            return
        reip_values = result
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_path, reip_values.shape[1], reip_values.shape[0], 1, gdal.GDT_Float32)
    out_ds.GetRasterBand(1).WriteArray(reip_values)
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
    out_ds = None
    layer_name = f"{csv_name}_{index_name}_Layer_coord_{coord_idx}"
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == layer_name:
            QgsProject.instance().removeMapLayer(layer.id())
    new_layer = QgsRasterLayer(output_path, layer_name)
    if new_layer.isValid():
        QgsProject.instance().addMapLayer(new_layer)
        style_raster_layer(new_layer, index_name)
        print(f"{index_name} layer coord {coord_idx} created/updated! Took {time.time() - start_time:.2f} seconds")
    else:
        print(f"Error creating {index_name} layer coord {coord_idx}")

# Функція для unmixing
def calculate_unmixing(reflectance_data, wavelengths, output_dir, geotransform, coord_idx, csv_name):
    start_time = time.time()
    reflectance_data = np.rot90(reflectance_data, k=1)
    endmembers = np.array([
        [0.05, 0.1],  # Flower: 500 nm, 700 nm
        [0.1, 0.2],   # Leaf: 500 nm, 700 nm
        [0.3, 0.4]    # Soil: 500 nm, 700 nm
    ]).T
    band_500 = find_band(500, wavelengths)
    band_700 = find_band(700, wavelengths)
    if band_500 is None or band_700 is None:
        print(f"Error: Invalid band indices for unmixing. 500 nm: {band_500}, 700 nm: {band_700}")
        return
    try:
        observed = reflectance_data[:, :, [band_500-1, band_700-1]]
        print(f"Observed shape: {observed.shape}, Endmembers shape: {endmembers.shape}")
    except IndexError as e:
        print(f"Error accessing bands: {e}")
        return
    if observed.shape[2] != endmembers.shape[0]:
        print(f"Error: Dimension mismatch in unmixing. Observed bands: {observed.shape[2]}, Endmembers bands: {endmembers.shape[0]}")
        return
    if np.any(np.isnan(observed)) or np.any(np.isinf(observed)):
        print("Error: Observed data contains NaN or infinite values")
        return
    try:
        observed_reshaped = observed.reshape(-1, observed.shape[2]).T
        print(f"Observed reshaped (transposed): {observed_reshaped.shape}")
        fractions, _, _, _ = np.linalg.lstsq(endmembers, observed_reshaped, rcond=None)
        flower_frac = fractions[0, :].reshape(observed.shape[:2])
    except np.linalg.LinAlgError as e:
        print(f"Error in unmixing: {e}")
        return
    output_path = os.path.join(output_dir, f"{csv_name}_Flower_Fraction_Unmixing_coord_{coord_idx}.tif")
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_path, flower_frac.shape[1], flower_frac.shape[0], 1, gdal.GDT_Float32)
    out_ds.GetRasterBand(1).WriteArray(flower_frac)
    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]')
    out_ds = None
    layer_name = f"{csv_name}_Unmixing_Flower_Fraction_coord_{coord_idx}"
    for layer in QgsProject.instance().mapLayers().values():
        if layer.name() == layer_name:
            QgsProject.instance().removeMapLayer(layer.id())
    unmix_layer = QgsRasterLayer(output_path, layer_name)
    if unmix_layer.isValid():
        QgsProject.instance().addMapLayer(unmix_layer)
        style_raster_layer(unmix_layer, "Unmixing")
        print(f"Unmixing layer coord {coord_idx} created/updated! Took {time.time() - start_time:.2f} seconds")
    else:
        print(f"Error creating unmixing layer coord {coord_idx}")

# Функція для забезпечення CSV-шару зверху
def ensure_csv_layer_on_top(csv_layer_name):
    root = QgsProject.instance().layerTreeRoot()
    csv_layer = None
    for layer in root.children():
        if layer.name() == csv_layer_name:
            csv_layer = layer
            break
    if csv_layer:
        cloned_layer = csv_layer.clone()
        root.insertChildNode(0, cloned_layer)
        root.removeChildNode(csv_layer)
        print(f"Moved CSV layer {csv_layer_name} to top")

# Основна логіка
def main():
    start_time = time.time()
    print("=== STARTING MAIN ===")
    
    for nc_path, csv_path in nc_csv_pairs.items():
        print(f"\nProcessing NetCDF: {nc_path}, CSV: {csv_path}")
        # Зчитування NetCDF
        if not os.path.exists(nc_path):
            print(f"Error: NetCDF file {nc_path} not found! Skipping.")
            continue
        try:
            ds = nc.Dataset(nc_path)
            wavelengths = ds['sensor_band_parameters/wavelengths'][:]
            reflectance_data = ds['reflectance'][:, :, :]
            print(f"Reflectance data shape: {reflectance_data.shape}, Wavelengths shape: {wavelengths.shape}")
        except Exception as e:
            print(f"Error reading NetCDF {nc_path}: {e}")
            continue
        
        # Завантаження CSV шару
        csv_name = os.path.splitext(os.path.basename(csv_path))[0]
        csv_layer = get_delimited_text_layer(csv_path)
        if not csv_layer:
            print(f"Processing aborted for {csv_path}: No valid CSV layer!")
            continue
        
        # Витягнення дати NetCDF із назви файлу (для дебагу)
        nc_date = os.path.basename(nc_path).split('_')[3][:8]  # Наприклад, '20230107'
        print(f"NetCDF date: {nc_date}")
        
        # Витягнення дат і координат
        date_from, date_to, unique_locations = extract_dates_locations(csv_layer, nc_date)
        print(f"Selected CSV: {csv_name}, Date range: {date_from} to {date_to}, Unique Locations: {len(unique_locations)}")
        if not unique_locations:
            print(f"Processing aborted for {csv_path}: No valid coordinates found!")
            continue
        
        # Створення вихідної директорії
        output_dir = os.path.join(output_base_dir, csv_name)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Обробка кожної координати
        for idx, (lat, lon) in enumerate(unique_locations):
            coord_start = time.time()
            print(f"\nProcessing unique coordinate {idx}: Lat {lat:.6f}, Lon {lon:.6f}")
            geotransform = calculate_geotransform(lat, lon, reflectance_data.shape[1], reflectance_data.shape[0])
            layer_names = [f"{csv_name}_{name}_Layer_coord_{idx}" for name in all_indices.keys()] + \
                          [f"{csv_name}_Unmixing_Flower_Fraction_coord_{idx}", 
                           f"{csv_name}_NetCDF_RGB_Layer_coord_{idx}",
                           f"{csv_name}_Georeferenced_TIFF_coord_{idx}"]
            for layer in QgsProject.instance().mapLayers().values():
                if layer.name() in layer_names:
                    QgsProject.instance().removeMapLayer(layer.id())
            
            # Обчислення індексів
            for name, data in all_indices.items():
                calculate_index(name, data['wl'], data['formula'], reflectance_data, wavelengths, output_dir, geotransform, idx, csv_name)
            
            # Unmixing
            calculate_unmixing(reflectance_data, wavelengths, output_dir, geotransform, idx, csv_name)
            
            # RGB-шар
            add_netcdf_rgb_layer(reflectance_data, wavelengths, output_dir, geotransform, idx, csv_name)
            
            # Геореференсування TIFF
            gdal_ds = gdal.Open(tiff_path)
            if gdal_ds is None:
                print(f"Error: Failed to open TIFF file {tiff_path}")
                tiff_width, tiff_height = reflectance_data.shape[1], reflectance_data.shape[0]
            else:
                tiff_width, tiff_height = gdal_ds.RasterXSize, gdal_ds.RasterYSize
                gdal_ds = None
            tiff_geotransform = calculate_geotransform(lat, lon, tiff_width, tiff_height)
            georeference_tiff(tiff_path, output_dir, tiff_geotransform, idx, csv_name)
            
            print(f"Coordinate {idx} took {time.time() - coord_start:.2f} seconds")
        
        # Переміщення CSV-шару на верх
        ensure_csv_layer_on_top(csv_name)
        
        ds.close()
        print(f"Finished processing pair: {nc_path}, {csv_path}")
    
    print(f"Processing completed! Total time: {time.time() - start_time:.2f} seconds")
    QgsProject.instance().layerTreeRoot().findLayerIds()  # Оновлення дерева шарів

if __name__ == "__main__":
    main()