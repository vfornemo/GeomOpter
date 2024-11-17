# Geometry Optimizer

Simple geometry optimizer program designed for saturated hydrocarbons.

## Prerequisites

### Rust

The program is written in Rust, so the Rust compiler must be installed. The installation instructions can be found [here](https://www.rust-lang.org/tools/install).

### OpenBLAS

OpenBLAS is an optimized BLAS library, required for the linear algebra operations in the program. The installation instructions can be found [here](https://www.openblas.net/).

## Configuration

The configuration files are located in `./config` directory, as follows:

### ctrl.toml

The control file is a `.toml` file that contains input arguments to the program, containing the following fields:

- `name` : Name of the job.
- `path` : Path to the input file.
- `calc_type` : Calculation type, which can be one of the following: `energy`, `car`(cartesian optimization), `intl`(internal optimization).
- `max_cycle` : Maximum number of cycles for geometry optimization. (Default: 200)
- `rms_tol` : RMS tolerance for the gradient in geometry optimization. (Default: 0.001)
- `car_conver` : Convergence criteria for the Cartesian coordinates in internal coordinate optimization. (Default: 0.00001)
- `log_level` : The level of the log, which can be one of the following: `off`, `error`, `warn`, `info`, `debug`, `trace`. (Default: `info`)
- `output` : Location of the output file. (Default: `./output.log`)

### log4rs.yaml

This file is no longer used, as the logging configuration is now done in the `ctrl.toml` file.

## Usage

After configuring the `ctrl.toml` file, the program can now be executed. For the first time, the program must be compiled with the following command:

```$ cargo build --release```

If there is an error message saying "linking with 'cc' failed: exit status: 1" and "ld: library not found for -lopenblas", then one can try the following command:

`RUSTFLAGS='-L/[path to openblas lib]' cargo build --release`

then one can run with the following command:

```$ ./target/release/GeomOpter ./config/ctrl.toml```

## Benchmark

| Molecule        | $N_\text{atoms}$|$t_\text{energy}\text{(ms)}$|$t_\text{car}\text{(ms)}$|$t_\text{intl}\text{(ms)}$|
|-----------------|-----------------|----------------------------|------------------------------|-------------------------------|
|methane          |      5          |            < 1             |             1                |          1                    |
|ethane           |      8          |            < 1             |             4                |          16                   |
|nbutane          |      14         |            < 1             |             21               |          45                   |
|isobutane        |      14         |            < 1             |             16               |          49                   |
|methylcyclohexane|      21         |            < 1             |             56               |          245                  |
|pinane           |      25         |            < 1             |             65               |          635                  |
|cholestane       |      75         |            < 1             |             2238             |          12322                |
