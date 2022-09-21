# Required files

Download the `VSOP2013.p2000` file from https://ftp.imcce.fr/pub/ephem/planets/vsop2013/ephemerides/, place it in the `./data` directory

Then run `cargo run --example convert --release --features=nom`

This'll generate the `VSOP2013.p2000.postcard` file which can then be loaded and used to get data
