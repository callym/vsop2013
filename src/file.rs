#[cfg(feature = "nom")]
use nom::{
  bytes::complete::tag,
  character::complete::{self, multispace0, space0},
  multi::count,
  number::complete::double,
  sequence::preceded,
  IResult,
};

use crate::Planet;

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct File {
  pub header: Header,
  pub tables: Vec<Table>,
}

#[derive(Debug)]
pub struct Calculation {
  pub x: f64,
  pub y: f64,
  pub z: f64,
  pub x_vel: f64,
  pub y_vel: f64,
  pub z_vel: f64,
}

impl File {
  pub fn from_bytes(bytes: &[u8]) -> Result<Self, postcard::Error> {
    postcard::from_bytes(bytes)
  }

  pub fn as_bytes(self) -> Result<Vec<u8>, postcard::Error> {
    postcard::to_allocvec(&self)
  }

  fn get_table(&self, jde: f64) -> Option<&Table> {
    if jde < self.header.start || jde >= self.header.end {
      return None;
    }

    let offset = jde - self.header.start;
    let offset = offset / self.header.interval;
    let offset = offset as usize;

    let table = self.tables.get(offset);

    #[cfg(debug_assertions)]
    if let Some(table) = table {
      assert!(jde >= table.start);
      assert!(jde < table.end);
    }

    table
  }

  pub fn calculate(&self, planet: Planet, jde: f64) -> Option<Calculation> {
    let table = self.get_table(jde)?;

    let planet_rank = self.header.planet_rank[planet as usize];
    // make 0-indexed
    let planet_rank = planet_rank - 1;

    let coefficient_number = self.header.planet_coefficients[planet as usize];
    let sub_interval_number = self.header.planet_sub_intervals[planet as usize];

    let interval = self.header.interval as u64;
    let delta = interval / sub_interval_number;

    let ik = jde - table.start;
    let ik = ik as u64;
    let ik = ik / delta;

    let ik = if ik == sub_interval_number {
      ik - 1
    } else {
      ik
    };

    let iloc = planet_rank + (6 * coefficient_number * ik);
    let dj0 = table.start + (ik as f64 * delta as f64);

    let x = (2.0 * (jde - dj0) / delta as f64) - 1.0;
    let mut tn = vec![0.0; coefficient_number as usize];
    tn[0] = 1.0;
    tn[1] = x;

    for i in 2..coefficient_number {
      let i = i as usize;
      tn[i] = (2.0 * x * tn[i - 1]) - tn[i - 2];
    }

    let mut r = vec![0.0; 6];

    for i in 0..6 {
      for j in 0..coefficient_number {
        let jp = coefficient_number - j - 1;
        let jt = iloc + (coefficient_number * i) + jp;

        let i = i as usize;
        let jp = jp as usize;
        let jt = jt as usize;

        r[i] = r[i] + (tn[jp] * table.coords[jt as usize]);
      }
    }

    let calc = Calculation {
      x: r[0],
      y: r[1],
      z: r[2],
      x_vel: r[3],
      y_vel: r[4],
      z_vel: r[5],
    };

    Some(calc)
  }
}

#[cfg(feature = "nom")]
impl File {
  pub fn parse(i: &str) -> IResult<&str, Self> {
    let (i, header) = Header::parse(i)?;
    let (i, _) = multispace0(i)?;

    let (i, tables) = count(
      preceded(multispace0, |i| Table::parse(&header, i)),
      header.tables as usize,
    )(i)?;

    Ok((i, Self { header, tables }))
  }
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Header {
  pub start: f64,
  pub end: f64,
  pub interval: f64,
  pub tables: u64,
  pub coefficients: u64,
  pub planet_rank: [u64; 9],
  pub planet_coefficients: [u64; 9],
  pub planet_sub_intervals: [u64; 9],
}

#[cfg(feature = "nom")]
impl Header {
  fn parse(i: &str) -> IResult<&str, Header> {
    let (i, _) = preceded(multispace0, tag("2013"))(i)?;
    let (i, start) = preceded(multispace0, double)(i)?;
    let (i, end) = preceded(multispace0, double)(i)?;

    let (i, interval) = preceded(multispace0, double)(i)?;
    let (i, tables) = preceded(multispace0, complete::u64)(i)?;
    let (i, coefficients) = preceded(multispace0, complete::u64)(i)?;

    let (i, planet_rank) = count(preceded(multispace0, complete::u64), 9)(i)?;
    let planet_rank = planet_rank.try_into().unwrap();

    let (i, planet_coefficients) = count(preceded(multispace0, complete::u64), 9)(i)?;
    let planet_coefficients = planet_coefficients.try_into().unwrap();

    let (i, planet_sub_intervals) = count(preceded(multispace0, complete::u64), 9)(i)?;
    let planet_sub_intervals = planet_sub_intervals.try_into().unwrap();

    Ok((
      i,
      Self {
        start,
        end,
        interval,
        tables,
        coefficients,
        planet_rank,
        planet_coefficients,
        planet_sub_intervals,
      },
    ))
  }
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Table {
  pub start: f64,
  pub end: f64,
  pub coords: Vec<f64>,
}

#[cfg(feature = "nom")]
impl Table {
  pub fn parse<'a>(header: &'_ Header, i: &'a str) -> IResult<&'a str, Self> {
    let (i, start) = preceded(space0, double)(i)?;
    let (i, end) = preceded(space0, double)(i)?;

    let (i, _) = multispace0(i)?;

    let lines = header.coefficients / 6;

    let (i, coords) = count(
      preceded(multispace0, count(preceded(space0, parse_coord), 6)),
      lines as usize,
    )(i)?;
    let coords = coords.into_iter().flatten().collect::<Vec<_>>();

    Ok((i, Self { start, end, coords }))
  }
}

#[cfg(feature = "nom")]
fn parse_coord(i: &str) -> IResult<&str, f64> {
  let (i, float) = preceded(space0, double)(i)?;
  let (i, int) = preceded(space0, complete::i8)(i)?;

  let e = 10_f64.powi(int as i32);

  Ok((i, float * e))
}

#[cfg(test)]
mod get_tables {
  use super::*;

  fn get_file() -> File {
    let manifest_dir = std::env!("CARGO_MANIFEST_DIR");
    let path = format!("{}/data/VSOP2013.p2000.postcard", manifest_dir);
    let bytes = std::fs::read(path).unwrap();

    File::from_bytes(&bytes).unwrap()
  }

  #[test]
  fn out_of_range_start() {
    let file = get_file();

    assert!(file.get_table(file.header.start - 1.0).is_none());
  }

  #[test]
  fn out_of_range_end() {
    let file = get_file();

    assert!(file.get_table(file.header.end + 1.0).is_none());
  }

  #[test]
  fn header_start() {
    let file = get_file();

    assert!(file.get_table(file.header.start).is_some());
  }

  #[test]
  fn header_end() {
    let file = get_file();

    assert!(file.get_table(file.header.end - 0.1).is_some());
  }

  #[test]
  fn middle() {
    let file = get_file();

    assert!(file.get_table(2624600.0).is_some());
  }
}

#[cfg(test)]
mod compute {
  use super::*;

  fn get_file() -> File {
    let manifest_dir = std::env!("CARGO_MANIFEST_DIR");
    let path = format!("{}/data/VSOP2013.p2000.postcard", manifest_dir);
    let bytes = std::fs::read(path).unwrap();

    File::from_bytes(&bytes).unwrap()
  }

  fn eq(a: f64, b: f64) {
    let a = format!("{a:>17.12}");
    // let a = format!("{:.15}", a);
    let b = format!("{:>17.12}", b);

    assert_eq!(a, b);
  }

  #[test]
  fn mercury_1() {
    let file = get_file();
    let res = file.calculate(Planet::Mercury, 2268932.5).unwrap();

    eq(res.x, 0.126270312608);
    eq(res.y, -0.430594516679);
    eq(res.z, -0.046642617731);

    eq(res.x_vel, 0.021376105515);
    eq(res.y_vel, 0.009431691610);
    eq(res.z_vel, -0.001225191385);
  }

  #[test]
  fn mercury_2() {
    let file = get_file();
    let res = file.calculate(Planet::Mercury, 2405730.5).unwrap();

    eq(res.x, 0.242020971329);
    eq(res.y, -0.352713705683);
    eq(res.z, -0.051047411323);

    eq(res.x_vel, 0.017592067303);
    eq(res.y_vel, 0.017315449357);
    eq(res.z_vel, -0.000208715093);
  }

  #[test]
  fn venus_1() {
    let file = get_file();
    let res = file.calculate(Planet::Venus, 2268932.5).unwrap();

    eq(res.x, 0.285941515135);
    eq(res.y, -0.669221964935);
    eq(res.z, -0.024828670826);

    eq(res.x_vel, 0.018458097789);
    eq(res.y_vel, 0.007874211087);
    eq(res.z_vel, -0.000975617536);
  }
}
