# Geometry Optimizer

Simple geometry optimizer program designed for saturated hydrocarbons.

## Configuration

The configuration files are located in `./config` directory, as follows:

### ctrl.toml

The control file is a `.toml` file that contains input arguments to the program, containing the following fields:

- `path` : Path to the input file.
- `calc_type` : Calculation type, which can be one of the following: `energy`, `car`(cartesian optimization), `intl`(internal optimization).
- `max_cycle` : Maximum number of cycles for geometry optimization.
- `rms_tol` : RMS tolerance for the gradient in geometry optimization.
- `car_conver` : Convergence criteria for the Cartesian coordinates in internal coordinate optimization.

### log4rs.yaml

The logging configuration file is a `.yaml` file, which controls the output log. There are several keywords that can be configured in `root` field:

`level` : The level of the log, which can be one of the following: `off`, `error`, `warn`, `info`, `debug`, `trace`.
`appenders` : The output of the log, which can be the following:

- `stdout` : Standard output in terminal.
- `file` : Output to a file. The output file name can be adjusted in `appenders/file/path`

## Usage

After configuring the `ctrl.toml` file and `log4rs.yaml`, the program can be run with the following command:

```$ ./target/release/GeomOpter [ctrl.toml]```
