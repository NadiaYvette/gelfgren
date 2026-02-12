# Gelfgren SWI-Prolog Bindings

SWI-Prolog bindings for Gelfgren using foreign predicate interface.

## Requirements

- SWI-Prolog 8.0+
- Gelfgren C library

## Installation

```prolog
?- pack_install(gelfgren).
```

## Usage

```prolog
:- use_module(library(gelfgren)).

?- bernstein_create([1.0, 2.0, 6.0], 0.0, 1.0, Poly),
   bernstein_evaluate(Poly, 0.5, Y),
   bernstein_free(Poly).
```

## Running Examples

```bash
export LD_LIBRARY_PATH=../../target/release:$LD_LIBRARY_PATH
swipl examples/basic_usage.pl
```

## License

MIT or Apache-2.0
