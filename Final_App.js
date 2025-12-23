function getSeasonDates(year, season) {
  var dates = {
    'Full Year': [ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year, 12, 31)],
    'Kharif (Jun-Oct)': [ee.Date.fromYMD(year, 6, 1), ee.Date.fromYMD(year, 10, 31)],
    'Rabi (Nov-Mar)': [ee.Date.fromYMD(year, 11, 1), ee.Date.fromYMD(year + 1, 3, 31)],
    'Pre-monsoon (Apr-May)': [ee.Date.fromYMD(year, 4, 1), ee.Date.fromYMD(year, 5, 31)]
  };
  return dates[season] || dates['Full Year'];
}

// Enhanced cloud masking for Sentinel-2
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var scl = image.select('SCL');
  
  // Bits 10 and 11 are clouds and cirrus
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  
  // QA60 cloud mask
  var qaMask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  
  // SCL cloud mask (if available)
  var sclMask = image.select('SCL').neq(3)  // Cloud shadows
    .and(image.select('SCL').neq(8))   // Cloud medium probability
    .and(image.select('SCL').neq(9))   // Cloud high probability
    .and(image.select('SCL').neq(10)); // Cirrus
  
  // Combine masks
  var finalMask = qaMask.and(sclMask);
  
  return image.updateMask(finalMask)
    .divide(10000)
    .select(image.bandNames().removeAll(['QA60', 'SCL']))
    .copyProperties(image, ['system:time_start']);
}

// Calculate vegetation indices
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
    {
      'NIR': image.select('B8'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
    }
  ).rename('EVI');
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5',
    {
      'NIR': image.select('B8'),
      'RED': image.select('B4')
    }
  ).rename('SAVI');
  
  return image.addBands([ndvi, evi, ndwi, savi]);
}

// Process Sentinel-2 collection
function processS2Collection(startDate, endDate) {
  var collection = ee.ImageCollection(CONFIG.s2Collection)
    .filterBounds(selectedAOI)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CONFIG.cloudThreshold))
    .map(maskS2clouds)
    .map(addIndices)
    .select(['NDVI', 'EVI', 'NDWI', 'SAVI']);
  
  return collection;
}
var YEARS = [2018,2019,2020]; // UI and exports use these
var DEFAULT_SCALE = 500;
var DEFAULT_DRIVE_FOLDER = 'GEE_Exports';

// ----------------- 1) AOI Setup -----------------
// ----------------- 1) AOI Setup -----------------
var AOIs = {
  // existing
  "Faridpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
                .filter(ee.Filter.eq('ADM2_NAME','Faridpur')),
  "Barisal": ee.FeatureCollection("FAO/GAUL/2015/level1")
                .filter(ee.Filter.eq('ADM1_NAME','Barisal')),
  "Khulna": ee.FeatureCollection("FAO/GAUL/2015/level1")
                .filter(ee.Filter.eq('ADM1_NAME','Khulna')),

  // new AOIs — mostly ADM2 (districts)
  "Dhaka": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Dhaka')),
  "Gazipur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Gazipur')),
  "Gopalganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Gopalganj')),
  "Jamalpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Jamalpur')),
  "Kishoreganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Kishoreganj')),
  "Madaripur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Madaripur')),
  "Manikganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Manikganj')),
  "Munshiganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Munshiganj')),
  "Mymensingh": ee.FeatureCollection("FAO/GAUL/2015/level1")
              .filter(ee.Filter.eq('ADM1_NAME','Mymensingh')),
  "Narayanganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Narayanganj')),
  "Narsingdi": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Narsingdi')),
  "Netrokona": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Netrokona')),
  "Rajbari": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Rajbari')),
  "Shariatpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Shariatpur')),
  "Sherpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Sherpur')),
  "Tangail": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Tangail')),
  "Bogra": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Bogra')),
  "Joypurhat": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Joypurhat')),
  "Naogaon": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Naogaon')),
  "Natore": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Natore')),
  "Nawabganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Nawabganj')),
  "Pabna": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Pabna')),
  "Rajshahi": ee.FeatureCollection("FAO/GAUL/2015/level1")
              .filter(ee.Filter.eq('ADM1_NAME','Rajshahi')),
  "Sirajgonj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Sirajganj')),
  "Dinajpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Dinajpur')),
  "Gaibandha": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Gaibandha')),
  "Kurigram": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Kurigram')),
  "Lalmonirhat": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Lalmonirhat')),
  "Nilphamari": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Nilphamari')),
  "Panchagarh": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Panchagarh')),
  "Rangpur": ee.FeatureCollection("FAO/GAUL/2015/level1")
              .filter(ee.Filter.eq('ADM1_NAME','Rangpur')),
  "Thakurgaon": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Thakurgaon')),
  "Barguna": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Barguna')),
  "Bhola": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Bhola')),
  "Jhalokati": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Jhalokati')),
  "Patuakhali": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Patuakhali')),
  "Pirojpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Pirojpur')),
  "Bandarban": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Bandarban')),
  "Brahmanbaria": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Brahmanbaria')),
  "Chandpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Chandpur')),
  "Chittagong": ee.FeatureCollection("FAO/GAUL/2015/level1")
              .filter(ee.Filter.eq('ADM1_NAME','Chittagong')),
  "Comilla": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Comilla')),
  "CoxsBazar": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Feni')),
  "Khagrachari": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Khagrachhari')),
  "Lakshmipur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Lakshmipur')),
  "Noakhali": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Noakhali')),
  "Rangamati": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Rangamati')),
  "Habiganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Habiganj')),
  "Maulvibazar": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Maulvibazar')),
  "Sunamganj": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Sunamganj')),
  "Sylhet": ee.FeatureCollection("FAO/GAUL/2015/level1")
              .filter(ee.Filter.eq('ADM1_NAME','Sylhet')),
  "Bagerhat": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Bagerhat')),
  "Chuadanga": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Chuadanga')),
  "Jessore": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Jessore')),
  "Jhenaidah": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Jhenaidah')),
  "Kushtia": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Kushtia')),
  "Magura": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Magura')),
  "Meherpur": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Meherpur')),
  "Narail": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Narail')),
  "Satkhira": ee.FeatureCollection("FAO/GAUL/2015/level2")
              .filter(ee.Filter.eq('ADM2_NAME','Satkhira'))
};
var selectedAOI = AOIs["Faridpur"];
Map.setCenter(89.15, 23.81, 8);
Map.addLayer(selectedAOI, {color:'red'}, 'AOI');

// ----------------- 2) UI Panel -----------------
var root = ui.Panel({style:{width:'380px'}});
ui.root.insert(0, root);
root.add(ui.Label('Verdent_view(2018-2020)', {fontWeight:'bold', fontSize:'16px'}));

// AOI selector
var aoiSelect = ui.Select({
  items: Object.keys(AOIs),
  value: "Faridpur",
  onChange: function(name) {
    selectedAOI = AOIs[name];
    Map.clear();
    Map.addLayer(selectedAOI, {color:'red'}, 'AOI');
  }
});
root.add(ui.Label('Select AOI'));
root.add(aoiSelect);

// Year selector (2015-2020)
var yearItems = YEARS.map(function(y){return String(y);});
var yearSelect = ui.Select({items: yearItems, value: String(YEARS[YEARS.length-1])});
root.add(ui.Label('Select Year (2015–2020)'));
root.add(yearSelect);

// Category selector
var categories = [
  'Vegetation Trend (VCF - percent tree cover from Sentinel-2 NDVI >0.3)',
  'Vegetation Indices (MODIS NDVI - annual mean)',
  'Crops & Agri Health (Sentinel-2 NDVI median)',
  'Rainfall Proxy (MODIS NDVI mean)',
  'Mineral/Nutrients (MODIS MCD43A3 Albedo_BSA_shortwave mean)'
];
var categorySelect = ui.Select({items: categories, value: categories[0]});
root.add(ui.Label('Select Category'));
root.add(categorySelect);

// Dataset description panel (we'll update by clearing/adding)
var descPanel = ui.Panel({layout: ui.Panel.Layout.flow('vertical')});
descPanel.add(ui.Label('Dataset description:'));
var datasetDescLabel = ui.Label('', {fontSize:'12px', color: '#0b5fff'});
descPanel.add(datasetDescLabel);
root.add(descPanel);

// Scale input and Export folder name
var scaleBox = ui.Textbox({value: String(DEFAULT_SCALE)});
var folderBox = ui.Textbox({value: DEFAULT_DRIVE_FOLDER});
root.add(ui.Panel([ui.Label('Scale (m):'), scaleBox], ui.Panel.Layout.flow('horizontal')));
root.add(ui.Panel([ui.Label('Drive Folder:'), folderBox], ui.Panel.Layout.flow('horizontal')));

// Apply button
var applyBtn = ui.Button({label: 'Show Layer & Export', style: {stretch: 'horizontal'}});
root.add(applyBtn);

// Bulk export note (no automatic bulk to avoid spam — user can loop themselves)
root.add(ui.Label('Tip: Use Show Layer & Export per AOI/year. Bulk exporting many tasks may exceed your Task limits.', {fontSize:'11px', color:'#555'}));

// ----------------- 3) Palettes & Legend helper -----------------
var palettes = {
  Class: {palette: ['#ff0000','#ffff00','#00ff00'], min: 0, max: 2},
  VCF: {palette: ['#ffffcc','#a1d99b','#006400'], min: 0, max: 100},
  NDVI: {palette: ['#ffffe5','#a1d99b','#006d2c'], min: 0, max: 1},
  ETProxy: {palette: ['#f7fcf5','#a1d99b','#00441b'], min: 0, max: 1000},
  GPP: {palette: ['#ffffcc','#41ab5d','#006400'], min: 0, max: 2000},
  Phenology: {palette: ['#ffffcc','#fd8d3c','#e31a1c'], min: 250, max: 320},
  LC_Type1: {palette: ['#006400','#ffff00','#ff0000','#7f7f7f'], min: 0, max: 5},
  Mineral: {palette: ['#f7fcf0','#c7e9c0','#00441b'], min: 0, max: 0.3},
  RainfallProxy: {palette: ['#ffffe5','#a1d99b','#006d2c'], min: -1, max: 1}
};
var legendPanel = null;
function updateLegend(mainBandName) {
  if (legendPanel) {
    try { ui.root.remove(legendPanel); } catch(e) {}
    legendPanel = null;
  }
  legendPanel = ui.Panel({style:{position:'bottom-left', padding:'8px 12px', backgroundColor:'white'}});
  legendPanel.add(ui.Label('Legend: ' + mainBandName, {fontWeight:'bold'}));
  var palObj = palettes[mainBandName] || palettes.Class;
  var pal = palObj.palette;
  // small color gradient
  legendPanel.add(ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {bbox: [0,0,1,0.1], dimensions: '200x10', format: 'png', palette: pal},
    style: {stretch: 'horizontal', margin: '4px 0 4px 0'}
  }));
  // class explanation
  var cls = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
  cls.add(ui.Label('🔴 Low '));
  cls.add(ui.Label('🟡 Moderate ', {margin: '0px'}));
  cls.add(ui.Label('🟢 Good'));
  legendPanel.add(cls);
  ui.root.add(legendPanel);
}

// ----------------- 4) Utility: safe image builder -----------------
// Returns an ee.Image with a band named outName, either col.mean()*scaleFactor (if col has data)
// or a constant image -9999 (nodata) with that band name. Also returns obs_count band if asked.
function safeMeanWithObs(col, bandName, scaleFactor, outName) {
  scaleFactor = (scaleFactor === undefined) ? 1 : scaleFactor;
  var has = col.size();
  var img = ee.Image(ee.Algorithms.If(has.gt(0),
    col.select(bandName).mean().multiply(scaleFactor).rename(outName),
    ee.Image.constant(-9999).rename(outName)
  ));
  // obs count when band exists
  var obs = ee.Image(ee.Algorithms.If(has.gt(0),
    col.select(bandName).count().rename('obs_count'),
    ee.Image.constant(0).rename('obs_count')
  ));
  return img.addBands(obs);
}

// ----------------- 5) Dataset functions (hardened) -----------------

// S2 cloud mask (SCL-based for S2_SR)
function maskS2SCL(img){
  var scl = img.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return img.updateMask(mask).copyProperties(img, ['system:time_start']);
}

// 1) Vegetation Trend — percent of pixels NDVI>0.3 (Sentinel-2), returns VCF (%)
function getVegetationTrend(year) {
  datasetDescLabel.setValue('Sentinel-2: percent of pixels with NDVI > 0.3 (vegetation fraction). 2015+ only.');
  var start = ee.Date.fromYMD(year,1,1);
  var end = ee.Date.fromYMD(year,12,31);
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterDate(start, end).filterBounds(selectedAOI).map(maskS2SCL)
            .map(function(i){ return i.addBands(i.normalizedDifference(['B8','B4']).rename('NDVI')); });
  var has = col.size();
  var out = ee.Image(ee.Algorithms.If(has.gt(0), (function(){
    var ndviCol = col.select('NDVI');
    var totalCount = ndviCol.count().rename('obs_count');
    // sum of binary ndvi>0.3 over collection:
    var binSum = ndviCol.map(function(i){ return i.gt(0.3).rename('v'); }).sum();
    var percent = binSum.divide(totalCount.where(totalCount.eq(0),1)).multiply(100).rename('VCF');
    percent = percent.where(totalCount.eq(0), -9999); // nodata where no obs
    var classImg = percent.expression("b('VCF')<10?0:b('VCF')<50?1:2").rename('Class');
    return percent.addBands(totalCount).addBands(classImg);
  })(), ee.Image.constant(-9999).rename('VCF').addBands(ee.Image.constant(0).rename('obs_count')).addBands(ee.Image.constant(0).rename('Class'))));
  return ee.Image(out).clip(selectedAOI);
}

// 2) Vegetation Indices (MODIS NDVI MCD13Q1) — NDVI annual mean
function getVegetationIndices(year) {
  datasetDescLabel.setValue('MODIS MCD13Q1: annual mean NDVI (scale 0.0001).');
  var start = ee.Date.fromYMD(year,1,1);
  var end = ee.Date.fromYMD(year,12,31);
  var col = ee.ImageCollection('MODIS/061/MOD13Q1').filterDate(start,end).filterBounds(selectedAOI);
  var safe = safeMeanWithObs(col, 'NDVI', 0.0001, 'NDVI');
  // class
  var cls = safe.select('NDVI').expression("b('NDVI')<-998?0:b('NDVI')<0.2?0:b('NDVI')<0.5?1:2").rename('Class');
  return safe.addBands(cls).clip(selectedAOI);
}

// 3) Crops & Ag Health (Sentinel-2 NDVI median)
function getCrops(year) {
  datasetDescLabel.setValue('Sentinel-2 L2A/SR: NDVI median (cloud-masked). 2015+ only.');
  var start = ee.Date.fromYMD(year,1,1);
  var end = ee.Date.fromYMD(year,12,31);
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterDate(start,end).filterBounds(selectedAOI).map(maskS2SCL)
            .map(function(i){ return i.addBands(i.normalizedDifference(['B8','B4']).rename('NDVI')); });
  var has = col.size();
  var out = ee.Image(ee.Algorithms.If(has.gt(0),
    col.select('NDVI').median().rename('NDVI').addBands(col.select('NDVI').count().rename('obs_count')),
    ee.Image.constant(-9999).rename('NDVI').addBands(ee.Image.constant(0).rename('obs_count'))
  ));
  var cls = out.select('NDVI').expression("b('NDVI')<-998?0:b('NDVI')<0.2?0:b('NDVI')<0.5?1:2").rename('Class');
  return out.addBands(cls).clip(selectedAOI);
}

// 4) Rainfall Proxy (MODIS NDVI mean)
function getRainfallProxy(year) {
  datasetDescLabel.setValue('Rainfall proxy: MODIS NDVI annual mean (vegetation response).');
  var start = ee.Date.fromYMD(year,1,1);
  var end = ee.Date.fromYMD(year,12,31);
  var col = ee.ImageCollection('MODIS/061/MOD13Q1').filterDate(start,end).filterBounds(selectedAOI);
  var safe = safeMeanWithObs(col, 'NDVI', 0.0001, 'RainfallProxy');
  var cls = safe.select('RainfallProxy').expression("b('RainfallProxy')<-998?0:b('RainfallProxy')<0.2?0:b('RainfallProxy')<0.5?1:2").rename('Class');
  return safe.addBands(cls).clip(selectedAOI);
}

// 10) Mineral / Nutrients (MCD43A3 Albedo_BSA_shortwave)
function getMinerals(year) {
  datasetDescLabel.setValue('MCD43A3 Albedo: Albedo_BSA_shortwave mean as proxy.');
  var col = ee.ImageCollection('MODIS/061/MCD43A3').filterDate(ee.Date.fromYMD(year,1,1), ee.Date.fromYMD(year,12,31)).filterBounds(selectedAOI).select('Albedo_BSA_shortwave');
  var has = col.size();
  var safe = ee.Image(ee.Algorithms.If(has.gt(0),
    col.mean().rename('Mineral').addBands(col.select('Albedo_BSA_shortwave').count().rename('obs_count')),
    ee.Image.constant(-9999).rename('Mineral').addBands(ee.Image.constant(0).rename('obs_count'))
  ));
  var cls = safe.select('Mineral').expression("b('Mineral')<-998?0:b('Mineral')<0.08?0:b('Mineral')<0.18?1:2").rename('Class');
  return safe.addBands(cls).clip(selectedAOI);
}

// ----------------- 6) Export helper -----------------
// mainBandName: string name of main band in img (e.g., 'NDVI', 'VCF', 'ETProxy', ...)
function exportAllFormats(img, mainBandName, catName, year, aoiName) {
  var folder = folderBox.getValue() || DEFAULT_DRIVE_FOLDER;
  var scale = parseInt(scaleBox.getValue(), 10) || DEFAULT_SCALE;
  var base = catName.replace(/\s+/g,'') + '' + year + '' + aoiName.replace(/\s+/g,'');

  // Ensure mainBand exists (img may be constant with -9999). We assume functions always created the band.
  // 1) Export main band (GeoTIFF)
  Export.image.toDrive({
    image: img.select(mainBandName).toFloat(),
    description: base + '_main',
    folder: folder,
    fileNamePrefix: base + '_main',
    region: selectedAOI.geometry(),
    scale: scale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });

  // 2) Export Class band (GeoTIFF)
  Export.image.toDrive({
    image: img.select('Class').toInt(),
    description: base + '_class',
    folder: folder,
    fileNamePrefix: base + '_class',
    region: selectedAOI.geometry(),
    scale: scale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });

  // 3) Export obs_count (if present)
  if (img.bandNames().contains('obs_count')) {
    Export.image.toDrive({
      image: img.select('obs_count').toInt(),
      description: base + '_obs_count',
      folder: folder,
      fileNamePrefix: base + '_obs_count',
      region: selectedAOI.geometry(),
      scale: scale,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });
  }

  // 4) Export CSV stats (min,max,mean)
  var stats = img.select(mainBandName).reduceRegion({
    reducer: ee.Reducer.min().combine(ee.Reducer.max(), '', true).combine(ee.Reducer.mean(), '', true),
    geometry: selectedAOI.geometry(),
    scale: scale,
    maxPixels: 1e13
  });
  var statsFeat = ee.Feature(null, stats).set({ 'Category': catName, 'Year': year, 'AOI': aoiName });
  Export.table.toDrive({
    collection: ee.FeatureCollection([statsFeat]),
    description: base + '_stats',
    folder: folder,
    fileNamePrefix: base + '_stats',
    fileFormat: 'CSV'
  });

  // 5) Export thumbnail: use visualize with palette if available
  var vis = palettes[mainBandName] || palettes.NDVI;
  var visParams = {min: vis.min, max: vis.max, palette: vis.palette};
  var thumb = img.select(mainBandName).visualize(visParams);
  Export.image.toDrive({
    image: thumb,
    description: base + '_thumbnail',
    folder: folder,
    fileNamePrefix: base + '_thumbnail',
    region: selectedAOI.geometry(),
    scale: Math.max(30, scale * 4), // downsampled quicklook
    maxPixels: 1e13
  });

  print('\u2705 Export scheduled:');
}

// ----------------- 7) Pixel inspector (map click) -----------------
var lastDisplayedImage = null;
Map.onClick(function(coords){
  if (!lastDisplayedImage) {
    print('No layer displayed yet. Use Show Layer & Export first.');
    return;
  }
  var pt = ee.Geometry.Point(coords.lon, coords.lat);
  var scale = parseInt(scaleBox.getValue(), 10) || DEFAULT_SCALE;
  var sample = lastDisplayedImage.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: pt,
    scale: scale,
    maxPixels: 1e9
  });
  print('Pixel at');
});

// ----------------- 8) Apply button behavior -----------------
applyBtn.onClick(function(){
  var year = parseInt(yearSelect.getValue(), 10);
  var cat = categorySelect.getValue();
  var img, mainBandName;

  // map category to function
  if (cat.indexOf('Vegetation Trend') === 0) { img = getVegetationTrend(year); mainBandName = 'VCF'; }
  else if (cat.indexOf('Vegetation Indices') === 0) { img = getVegetationIndices(year); mainBandName = 'NDVI'; }
  else if (cat.indexOf('Crops') === 0) { img = getCrops(year); mainBandName = 'NDVI'; }
  else if (cat.indexOf('Rainfall Proxy') === 0) { img = getRainfallProxy(year); mainBandName = 'RainfallProxy'; }
  else if (cat.indexOf('Surface Groundwater') === 0) { img = getGroundwaterProxy(year); mainBandName = 'ETProxy'; }
  else if (cat.indexOf('Photosynthesis') === 0) { img = getPhotosynthesis(year); mainBandName = 'GPP'; }
  else if (cat.indexOf('Phenology') === 0) { img = getPhenology(year); mainBandName = 'Phenology'; }
  else if (cat.indexOf('Deforestation') === 0) { img = getDeforestation(year); mainBandName = 'LC_Type1'; }
  else if (cat.indexOf('Water Evaporation') === 0) { img = getWaterEvap(year); mainBandName = 'ETProxy'; }
  else if (cat.indexOf('Mineral') === 0) { img = getMinerals(year); mainBandName = 'Mineral'; }
  else { ui.alert('Category not implemented'); return; }

  // Ensure returned is ee.Image
  img = ee.Image(img);

  // Display
  Map.clear();
  Map.addLayer(selectedAOI, {color:'red'}, 'AOI');
  var visObj = palettes[mainBandName] || palettes.NDVI;
  Map.addLayer(img.select(mainBandName), visObj, cat + ' ' + year, true, 0.85);
  Map.addLayer(img.select('Class'), palettes.Class, cat + ' Class', true, 0.6);

  // Update legend and lastDisplayedImage
  updateLegend(mainBandName);
  lastDisplayedImage = img;

  // Schedule exports to Drive for Colab ingestion
  var aoiName = aoiSelect.getValue();
  exportAllFormats(img, mainBandName, cat, year, aoiName);

  // Print basic stats to console
  var scale = parseInt(scaleBox.getValue(), 10) || DEFAULT_SCALE;
  var stats = img.select(mainBandName).reduceRegion({
    reducer: ee.Reducer.min().combine(ee.Reducer.max(), '', true).combine(ee.Reducer.mean(), '', true),
    geometry: selectedAOI.geometry(),
    scale: scale,
    maxPixels: 1e13
  });
  print(cat + ' ' + year + ' stats (min, max, mean):');
});

// add apply button at bottom (again)
root.add(ui.Label('')); // spacing


// initial legend
updateLegend('VCF');
print('Ready — choose AOI → year → category → click "Show Layer & Export".');
print('Exports will appear in the Tasks tab and then in your Google Drive folder.');

/* End of script */
