// CHIRPS v2 / v3-sat / v3-rnl daily -> one multiband GeoTIFF per year (1 band = 1 day)

var bbox = ee.FeatureCollection('projects/ee-ryonaba/assets/bbox_utm').geometry();

var START_YEAR = 2011;
var END_YEAR = 2020;

var products = [
  {name: 'chirps_v2',     id: 'UCSB-CHG/CHIRPS/DAILY',        band: 'precipitation'},
  {name: 'chirps_v3_sat', id: 'UCSB-CHC/CHIRPS/V3/DAILY_SAT', band: 'precipitation'},
  {name: 'chirps_v3_rnl', id: 'UCSB-CHC/CHIRPS/V3/DAILY_RNL', band: 'precipitation'}
];

Map.centerObject(bbox, 8);
Map.addLayer(bbox, {color: 'red'}, 'bbox');

products.forEach(function (p) {
  var col = ee.ImageCollection(p.id).select(p.band);

  // Native grid (CRS + affine transform) of the product
  var proj = col.first().projection();
  var crs = proj.crs().getInfo();
  var transform = proj.getInfo().transform;

  for (var year = START_YEAR; year <= END_YEAR; year++) {
    var start = ee.Date.fromYMD(year, 1, 1);
    var yearCol = col.filterDate(start, start.advance(1, 'year'));

    // Skip years with no data (e.g. CHIRPS v3 "sat" does not cover the early years)
    if (yearCol.size().getInfo() === 0) {
      print(p.name + ': no images in ' + year + ', skipped');
      continue;
    }

    // Rename each band after its date (YYYYMMDD)
    var dates = yearCol.aggregate_array('system:time_start')
      .map(function (t) { return ee.Date(t).format('YYYYMMdd'); });

    var img = yearCol.toBands()
      .rename(dates)
      .toFloat()
      .clip(bbox);

    Export.image.toDrive({
      image: img,
      description: p.name + '_' + year,
      folder: p.name,
      fileNamePrefix: p.name + '_' + year,
      region: bbox,
      crs: crs,
      crsTransform: transform,   // native resolution (0.05 deg), no resampling
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF'
    });
  }
});