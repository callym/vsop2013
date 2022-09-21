fn main() {
  let manifest_dir = std::env!("CARGO_MANIFEST_DIR");
  let path = format!("{}/data/VSOP2013.p2000", manifest_dir);
  let file = std::fs::read_to_string(path).unwrap();

  let (_, mut file) = vsop2013::File::parse(&file).unwrap();

  // let start = time::macros::datetime!(1900-01-01 0:00 UTC).to_julian_day() as f64;
  // let end = time::macros::datetime!(2100-01-01 0:00 UTC).to_julian_day() as f64;

  // file.header.start = start;
  // file.header.end = end;

  // file.tables = file
  //   .tables
  //   .into_iter()
  //   .filter(|table| table.start > start && table.end < end)
  //   .collect();

  let bytes = file.into_bytes().unwrap();
  let path = format!("{}/data/VSOP2013.p2000.postcard", manifest_dir);
  std::fs::write(path, bytes).unwrap();
}
