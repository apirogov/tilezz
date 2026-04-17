# Agent Instructions

## Development Workflow

* This is a Rust project, ALWAYS use `cargo fmt`, `cargo clippy` and `cargo test` and fix issues.
* When the code is formatted, linted and successfully tested, ALWAYS commit your changes.

## Agent Behavior

* ALWAYS ask if there is any ambiguity in the request or multiple viable solutions.
* NEVER guess, think through issues, read docs and existing code, discuss with the user.
* NEVER assume anything if you are not confident or need more details.
* If documentation is missing, try installing it or ask the user for help.

## Code Quality

* Avoid data duplication: NEVER hardcode a numeric value or string in two places, define constants.
* Avoid code duplication: NEVER write the same logic twice, factor out a common generic function.

* Before writing any new non-trivial code, ALWAYS check whether similar functionality already exists.
* Before writing any new non-trivial code, ALWAYS consider using existing popular third-party crates.
* When actually writing new code, ALWAYS try to write it in a reusable and extensible way from the start.

* ALWAYS strive towards generic, modular, maintainable code.
* ALWAYS prioritize good architecture and clarity over performance unless told otherwise..
* ALWAYS write testable code with easily mockable inputs.
* ALWAYS write small functions that ideally fit on a single screen.
* ALWAYS avoid deep nesting, use more helper functions instead.
* ALWAYS ask for permission to let a function be untested.
* ALWAYS plan out suitable tests for the intended design BEFORE implementation.
* ALWAYS strive to test all code paths, covering edge cases and boundary input values.
