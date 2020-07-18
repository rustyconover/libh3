use std::mem::MaybeUninit;

/// Convert degrees to radians
///
/// ```
/// use libh3::degs_to_rads;
/// assert_eq!(2.413790355508158, degs_to_rads(138.3));
/// ```
pub fn degs_to_rads(degrees: f64) -> f64 {
    unsafe { libh3_sys::degsToRads(degrees) }
}

/// Convert radians to degrees
///
/// ```
/// use libh3::rads_to_degs;
/// assert_eq!(138.3, rads_to_degs(2.413790355508158));
/// ```
pub fn rads_to_degs(radians: f64) -> f64 {
    unsafe { libh3_sys::radsToDegs(radians) }
}

/// Represent a coordinate
#[derive(Debug)]
pub struct GeoCoord {
    // The latitute of the coordinate, typcially this should be specified using
    // radians but it is easy to convert using [degs_to_rads](degs_to_rads)
    pub lat: f64,
    // The longitude of the coordinate, typcially this should be specified using
    // radians but it is easy to convert using [degs_to_rads](degs_to_rads)
    pub lon: f64,
}

impl GeoCoord {
    /// Create a new GeoCoord representing a coordinate
    ///
    /// # Arguments
    ///
    /// * `lat` - The latitude of the coordinate
    /// * `long` - The longitude of the coordinate
    ///
    pub fn new(lat: f64, lon: f64) -> GeoCoord {
        GeoCoord { lat, lon }
    }
}

impl From<&GeoCoord> for libh3_sys::GeoCoord {
    fn from(coord: &GeoCoord) -> Self {
        libh3_sys::GeoCoord {
            lat: coord.lat,
            lon: coord.lon,
        }
    }
}

/// A H3 index value a unique address of a hexagon or more unlikely
/// a pentagon.
pub type H3Index = libh3_sys::H3Index;

/// A resolution that ranges from 0 to 15.
///
/// See the [resolution table](https://h3geo.org/docs/core-library/restable)
/// for the sizes of resolution.
pub type Resolution = u8;

/// Return the edge length of a hexagon at a particular resolution in kilometers.
///
/// ```
/// use libh3::edge_length_km;
/// assert_eq!(edge_length_km(5), 8.544408276);
/// ```
pub fn edge_length_km(resolution: Resolution) -> f64 {
    unsafe { libh3_sys::edgeLengthKm(resolution as i32) }
}

/// Return the number of hexagons at a particular resolution.
///
/// ```
/// use libh3::num_hexagons;
/// assert_eq!(num_hexagons(5), 2016842);
/// ```
pub fn num_hexagons(resolution: Resolution) -> u64 {
    unsafe { libh3_sys::numHexagons(resolution as i32) as u64 }
}

/// Return the edge length of a hexagon at a particular resolution in meters.
///
/// ```
/// use libh3::edge_length_m;
/// assert_eq!(edge_length_m(5), 8544.408276);
/// ```
pub fn edge_length_m(resolution: Resolution) -> f64 {
    unsafe { libh3_sys::edgeLengthM(resolution as i32) }
}

/// Convert a GeoCoord to a H3 index.
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3};
/// let coords = GeoCoord {
///   lat: degs_to_rads(40.689167),
///   lon: degs_to_rads(-74.044444),
/// };
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(v.unwrap(), 0x8a2a1072b59ffff);
/// ```
pub fn geo_to_h3(coord: &GeoCoord, resolution: Resolution) -> Result<H3Index, ()> {
    unsafe {
        match libh3_sys::geoToH3(&libh3_sys::GeoCoord::from(coord), resolution as i32) {
            0 => Err(()),
            x => Ok(x),
        }
    }
}

/// Convert a H3 index value to a GeoCoord
///
/// ```
/// use libh3::h3_to_geo;
/// let r = h3_to_geo(0x8a2a1072b59ffff);
/// assert_eq!(r.lat, 0.7101643819054542);
/// assert_eq!(r.lon, -1.2923191206954798);
/// ```
pub fn h3_to_geo(h3: H3Index) -> GeoCoord {
    let result = unsafe {
        let mut result: MaybeUninit<libh3_sys::GeoCoord> = MaybeUninit::uninit().assume_init();
        libh3_sys::h3ToGeo(h3, result.as_mut_ptr());
        result.assume_init()
    };
    GeoCoord::new(result.lat, result.lon)
}

type GeoBoundary = Vec<GeoCoord>;

/// Convert a H3 index value to a GeoBoundary which are a
/// vector of points that describe a H3 Index's boundary
///
/// ```
/// use libh3::h3_to_geo_boundary;
/// let foo = h3_to_geo_boundary(0x8a2a1072b59ffff);
/// assert_eq!(foo.len(), 6);
/// ```
pub fn h3_to_geo_boundary(h3: H3Index) -> GeoBoundary {
    unsafe {
        let mut boundary_result: MaybeUninit<libh3_sys::GeoBoundary> =
            MaybeUninit::uninit().assume_init();
        libh3_sys::h3ToGeoBoundary(h3, boundary_result.as_mut_ptr());
        let boundary = boundary_result.assume_init();
        let mut result = Vec::with_capacity(boundary.numVerts as usize);
        for i in 0..boundary.numVerts as usize {
            result.push(GeoCoord::new(boundary.verts[i].lat, boundary.verts[i].lon));
        }
        result
    }
}

/// Return the resolution of a H3 index
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3, h3_get_resolution };
/// let coords = GeoCoord::new(
///   degs_to_rads(40.689167),
///   degs_to_rads(-74.044444),
/// );
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(h3_get_resolution(v.unwrap()), 10);
/// ````
pub fn h3_get_resolution(h3: H3Index) -> Resolution {
    unsafe { libh3_sys::h3GetResolution(h3) as Resolution }
}

/// Determine if H3 index is valid
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3, h3_is_valid};
/// let coords = GeoCoord {
///   lat: degs_to_rads(40.689167),
///   lon: degs_to_rads(-74.044444),
/// };
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(h3_is_valid(v.unwrap()),true);
/// ```
pub fn h3_is_valid(h3: H3Index) -> bool {
    unsafe {
        match libh3_sys::h3IsValid(h3) {
            0 => false,
            _ => true,
        }
    }
}

/// Determine if two H3 indexes are neighbors
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3, h3_indexes_are_neighbors};
/// let coords = GeoCoord {
///   lat: degs_to_rads(40.689167),
///   lon: degs_to_rads(-74.044444),
/// };
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(h3_indexes_are_neighbors(v.unwrap(), v.unwrap()), false);
/// ```
pub fn h3_indexes_are_neighbors(origin: H3Index, destination: H3Index) -> bool {
    unsafe {
        match libh3_sys::h3IndexesAreNeighbors(origin, destination) {
            0 => false,
            _ => true,
        }
    }
}

/// Determine the area of a hexagon at a particular resolution.
///
/// ```
/// use libh3::{ hex_area_km_2};
/// assert_eq!(hex_area_km_2(10), 0.0150475);
/// ```
pub fn hex_area_km_2(resolution: i32) -> f64 {
    unsafe { libh3_sys::hexAreaKm2(resolution) }
}

/// Determine if the specified H3 index is a pentagon.
/// ```
/// assert_eq!(libh3::h3_is_pentagon(0x8a2a1072b59ffff), false);
/// ```
pub fn h3_is_pentagon(h3: H3Index) -> bool {
    unsafe {
        match libh3_sys::h3IsPentagon(h3) {
            0 => false,
            _ => true,
        }
    }
}

/// Get the number of the base cell for a given H3 index
///
/// ```
/// use libh3;
/// assert_eq!(libh3::h3_get_base_cell(0x8a2a1072b59ffff), 21);
/// ```
pub fn h3_get_base_cell(h3: H3Index) -> i32 {
    unsafe { libh3_sys::h3GetBaseCell(h3) }
}

/// Get all hexagons in a k-ring around a given center. The order of the hexagons is undefined.
///
/// # Arguments
///
/// * `origin` - The center of the ring.
/// * `radius` - The radis of the ring in hexagons, which is the same resolution as the origin.
///
/// ```
/// let expected_kring = vec![
///   0x8a2a1072b59ffff,
///   0x8a2a1072b597fff,
///   0x8a2a1070c96ffff,
///   0x8a2a1072b4b7fff,
///   0x8a2a1072b4a7fff,
///   0x8a2a1072b58ffff,
///   0x8a2a1072b587fff,
/// ];
/// let r = libh3::k_ring(0x8a2a1072b59ffff, 1);
/// assert_eq!(r, expected_kring);
/// ```
pub fn k_ring(origin: H3Index, radius: i32) -> Vec<H3Index> {
    unsafe {
        let max = libh3_sys::maxKringSize(radius);
        let mut r = Vec::<H3Index>::with_capacity(max as usize);
        libh3_sys::kRing(origin, radius, r.as_mut_ptr());
        r.set_len(max as usize);
        r = r.into_iter().filter(|v| *v != 0).collect();
        r
    }
}

/// Get all hexagons in a k-ring around a given center, in an array of arrays
/// ordered by distance from the origin. The order of the hexagons within each ring is undefined.
///
/// # Arguments
///
/// * `origin` - The center of the ring.
/// * `radius` - The radis of the ring in hexagons, which is the same resolution as the origin.
///
/// ```
/// let expected_kring_distances = vec![
///   (0x8a2a1072b59ffff, 0),
///   (0x8a2a1072b597fff, 1),
///   (0x8a2a1070c96ffff, 1),
///   (0x8a2a1072b4b7fff, 1),
///   (0x8a2a1072b4a7fff, 1),
///   (0x8a2a1072b58ffff, 1),
///   (0x8a2a1072b587fff, 1),
/// ];
/// let r = libh3::k_ring_distances(0x8a2a1072b59ffff, 1);
/// assert_eq!(r, expected_kring_distances);
/// ```
pub fn k_ring_distances(origin: H3Index, radius: i32) -> Vec<(H3Index, i32)> {
    unsafe {
        let max = libh3_sys::maxKringSize(radius);
        let mut indexes = Vec::<H3Index>::with_capacity(max as usize);
        let mut distances = Vec::<i32>::with_capacity(max as usize);

        libh3_sys::kRingDistances(origin, radius, indexes.as_mut_ptr(), distances.as_mut_ptr());
        indexes.set_len(max as usize);
        distances.set_len(max as usize);

        indexes
            .into_iter()
            .zip(distances.into_iter())
            .filter(|v| v.0 != 0)
            .collect::<Vec<(H3Index, i32)>>()
    }
}

pub fn hex_range(origin: H3Index, k: i32) -> (bool, Vec<H3Index>) {
    unsafe {
        let max = libh3_sys::maxKringSize(k);
        let mut r = Vec::<H3Index>::with_capacity(max as usize);
        let distortion = libh3_sys::hexRange(origin, k, r.as_mut_ptr());
        r.set_len(max as usize);
        r = r.into_iter().filter(|v| *v != 0).collect();
        (distortion == 0, r)
    }
}

pub fn hex_range_distances(origin: H3Index, k: i32) -> (bool, Vec<(H3Index, i32)>) {
    unsafe {
        let max = libh3_sys::maxKringSize(k);
        let mut indexes = Vec::<H3Index>::with_capacity(max as usize);
        let mut distances = Vec::<i32>::with_capacity(max as usize);

        let distortion =
            libh3_sys::hexRangeDistances(origin, k, indexes.as_mut_ptr(), distances.as_mut_ptr());
        indexes.set_len(max as usize);
        distances.set_len(max as usize);

        (
            distortion == 0,
            indexes
                .into_iter()
                .zip(distances.into_iter())
                .filter(|v| v.0 != 0)
                .collect::<Vec<(H3Index, i32)>>(),
        )
    }
}

/// Get the grid distance between two hex indexes. This function may fail
/// to find the distance between two indexes if they are very far apart or
/// on opposite sides of a pentagon.
/// # Arguments
///
/// * `origin` - The starting H3 index
/// * `end` - The ending H3 index
///
/// ```
/// use libh3::h3_distance;
/// assert_eq!(h3_distance(0x8a2a1072b4a7fff, 0x8a2a1072b58ffff), Ok(1));
/// ```
pub fn h3_distance(origin: H3Index, end: H3Index) -> Result<i32, ()> {
    unsafe {
        let r = libh3_sys::h3Distance(origin, end);
        if r < 0 {
            Err(())
        } else {
            Ok(r)
        }
    }
}

/// Get all hexagons with centers contained in a given polygon. The polygon
/// is specified with GeoJson semantics as an array of loops. The first loop
/// is the perimeter of the polygon, and subsequent loops are
/// expected to be holes.
///
/// # Arguments
///
/// * `polygon` - The vector of polygons.
/// * `resolution` - The resolution of the generated hexagons
///
/// ```
/// use libh3::{polyfill, GeoCoord, degs_to_rads};
/// /// Some vertexes around San Francisco
/// let sf_verts = vec![
///   (0.659966917655, -2.1364398519396),
///   (0.6595011102219, -2.1359434279405),
///   (0.6583348114025, -2.1354884206045),
///   (0.6581220034068, -2.1382437718946),
///   (0.6594479998527, -2.1384597563896),
///   (0.6599990002976, -2.1376771158464),
/// ]
/// .iter()
/// .map(|v| GeoCoord::new(v.0, v.1))
/// .collect();
///
/// let h = polyfill(&vec![sf_verts], 9);
/// assert_eq!(h.len(), 1253);
///
///
/// /// Fill a polygon around Wellington, NZ
/// /// Coordinates are in GeoJSON lon, lat order.
///
/// let wellington_verts: Vec<GeoCoord> = vec![
///    (174.800937866947, -41.22501356278325),
///    (174.8079721211159, -41.226732341115365),
///    (174.82262997231396, -41.231639803277986),
///    (174.83561105648377, -41.23873115217201),
///    (174.84634815587896, -41.24769587535717),
///    (174.85069634833735, -41.252194466801384),
///    (174.8587207192276, -41.26264015566857),
///    (174.8636809159909, -41.27410982273725),
///    (174.86536017144866, -41.28610182569948),
///    (174.8653611411562, -41.29165993179072),
///    (174.8636858034186, -41.30364998147521),
///    (174.85872883206798, -41.31511423541933),
///    (174.8507068877321, -41.32555201329627),
///    (174.846359294586, -41.33004662498431),
///    (174.8356222320399, -41.33900220579377),
///    (174.82263983163466, -41.346084796073406),
///    (174.807979604528, -41.35098543989819),
///    (174.80094378133927, -41.35270175860709),
///    (174.78524618901284, -41.35520670405109),
///    (174.76919781098724, -41.35520670405109),
///    (174.75350021866086, -41.35270175860712),
///    (174.7464643954721, -41.35098543989822),
///    (174.7318041683653, -41.346084796073406),
///    (174.71882176795995, -41.33900220579369),
///    (174.7080847054138, -41.330046624984135),
///    (174.70373711226773, -41.32555201329609),
///    (174.69571516793187, -41.31511423541913),
///    (174.69075819658127, -41.303649981474955),
///    (174.68908285884382, -41.29165993179046),
///    (174.68908382855136, -41.2861018256992),
///    (174.69076308400918, -41.274109822737074),
///    (174.69572328077246, -41.26264015566849),
///    (174.70374765166264, -41.252194466801384),
///    (174.70809584412103, -41.24769587535717),
///    (174.71883294351622, -41.23873115217201),
///    (174.731814027686, -41.231639803278014),
///    (174.746471878884, -41.22673234111539),
///    (174.75350613305287, -41.22501356278328),
///    (174.7691998725514, -41.222504896122565),
///    (174.78524412744844, -41.222504896122565),
///    (174.800937866947, -41.22501356278325),
/// ].iter().map(|v| GeoCoord::new(degs_to_rads(v.1), degs_to_rads(v.0))).collect();
///
/// let mut h = polyfill(&vec![wellington_verts], 6);
/// assert_eq!(h.len(), 5);
/// h.sort_unstable();
/// assert_eq!(h, vec![606774924341673983, 606774925281198079, 606774925549633535, 606774929307729919, 606774929441947647]);
/// ```
pub fn polyfill(polygon: &[Vec<GeoCoord>], resolution: Resolution) -> Vec<H3Index> {
    let real_polygon = polygon
        .iter()
        .map(|p| p.iter().map(libh3_sys::GeoCoord::from).collect())
        .collect::<Vec<Vec<libh3_sys::GeoCoord>>>();

    unsafe {
        let fence = libh3_sys::Geofence {
            numVerts: real_polygon[0].len() as i32,
            verts: real_polygon[0].as_ptr(),
        };

        let holes = real_polygon
            .iter()
            .skip(1)
            .map(|p| libh3_sys::Geofence {
                numVerts: p.len() as i32,
                verts: p.as_ptr(),
            })
            .collect::<Vec<libh3_sys::Geofence>>();

        let p = libh3_sys::GeoPolygon {
            geofence: fence,
            numHoles: (real_polygon.len() - 1) as i32,
            holes: holes.as_ptr(),
        };

        let max = libh3_sys::maxPolyfillSize(&p, resolution as i32);
        let mut r = Vec::<H3Index>::with_capacity(max as usize);
        libh3_sys::polyfill(&p, resolution as i32, r.as_mut_ptr());
        r.set_len(max as usize);
        r.retain(|&v| v != 0);
        r
    }
}

/// Returns the size of the array needed by h3ToChildren for these inputs.
///
/// # Arguments
///
/// * `h` - The index of the parent resolution.
/// * `resolution` - The resolution of the desired level.
///
/// ```
/// use libh3::max_h3_to_children_size;
/// assert_eq!(max_h3_to_children_size(0x852a1073fffffff, 6), 7);
/// ```
pub fn max_h3_to_children_size(h: H3Index, resolution: Resolution) -> i32 {
    unsafe { libh3_sys::maxH3ToChildrenSize(h, resolution as i32) }
}

/// Returns the parent (coarser) index containing h3.
///
/// # Arguments
///
/// * `h` - The index of the child resolution.
/// * `resolution` - The resolution of the desired level.
///
/// ```
/// use libh3::h3_to_parent;
/// assert_eq!(h3_to_parent(0x8a2a1072b4a7fff, 5), 0x852a1073fffffff);
/// ```
pub fn h3_to_parent(h: H3Index, resolution: Resolution) -> H3Index {
    unsafe { libh3_sys::h3ToParent(h, resolution as i32) }
}

/// Returns children indexes contained by the given index at the given resolution.
///
/// # Arguments
///
/// * `h` - The index of the child resolution.
/// * `resolution` - The resolution of the desired level.
///
/// ```
/// use libh3::h3_to_children;
/// assert_eq!(
///     h3_to_children(0x852a1073fffffff, 6),
///     vec![
///         0x862a10707ffffff,
///         0x862a1070fffffff,
///         0x862a10717ffffff,
///         0x862a1071fffffff,
///         0x862a10727ffffff,
///         0x862a1072fffffff,
///         0x862a10737ffffff
///     ]
/// );
/// ```
pub fn h3_to_children(h: H3Index, resolution: Resolution) -> Vec<H3Index> {
    let max = max_h3_to_children_size(h, resolution) as usize;
    let mut result = Vec::<H3Index>::with_capacity(max as usize);
    unsafe {
        libh3_sys::h3ToChildren(h, resolution as i32, result.as_mut_ptr());
        result.set_len(max as usize);
        result
    }
}

/// Number of resolution 0 H3 indexes.
///
/// ```
/// use libh3::res_0_index_count;
/// assert_eq!(res_0_index_count(), 122);
/// ```
pub fn res_0_index_count() -> i32 {
    unsafe { libh3_sys::res0IndexCount() }
}

/// All the resolution 0 H3 indexes.
///
/// ```
/// use libh3::get_res_0_indexes;
/// assert_eq!(
///     get_res_0_indexes().len(),
///     122
/// );
/// ```
pub fn get_res_0_indexes() -> Vec<H3Index> {
    let max = res_0_index_count() as usize;
    let mut result = Vec::<H3Index>::with_capacity(max as usize);
    unsafe {
        libh3_sys::getRes0Indexes(result.as_mut_ptr());
        result.set_len(max);
        result
    }
}