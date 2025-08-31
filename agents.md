# Rust Development Best Practices

This document outlines basic best practices and essential commands for Rust development in this project.

## Code Formatting

Maintain consistent code formatting using Rust's built-in formatter:

```bash
cargo fmt
```

Run this command to format your code according to Rust's standard formatting guidelines. Always run `cargo fmt` before committing your changes to ensure consistent code style.

## Testing

Rust has built-in support for unit and integration tests. Run all tests with:

```bash
cargo test
```

This command will:
- Compile your code and dependencies
- Run all unit tests
- Run all integration tests (located in `tests/` directory)
- Provide detailed output including test results and coverage if configured

For more verbose output, use:
```bash
cargo test -- --nocapture
```

## Building

To build your project:

```bash
cargo build
```

For optimized release builds (recommended for production):
```bash
cargo build --release
```

## Running the Application

To run your Rust application in development mode:

```bash
cargo run
```

For release mode:
```bash
cargo run --release
```

## Dependencies

To add a new dependency to your `Cargo.toml`:

1. Edit the `Cargo.toml` file and add the dependency under `[dependencies]`
2. Run `cargo build` or `cargo update` to download and compile the new dependency

Example:
```toml
[dependencies]
serde = { version = "1.0", features = ["derive"] }
```

## Git Workflow

Follow these steps before committing your changes:

1. Format your code:
   ```bash
   cargo fmt
   ```

2. Run tests to ensure everything works:
   ```bash
   cargo test
   ```

3. Build the project:
   ```bash
   cargo build
   ```

4. Fix any compiler warnings or errors

5. Commit your changes:
   ```bash
   git add .
   git commit -m "Your descriptive commit message"
   ```

## Clippy (Linting)

Use Clippy to catch common mistakes and improve your code:

```bash
cargo clippy
```

Clippy is Rust's collection of lints to help ensure your code follows best practices. Run it regularly and fix any warnings it reports.

## Documentation

To generate documentation for your project:

```bash
cargo doc --open
```

This will generate HTML documentation and open it in your browser. Write documentation comments using `///` for functions and structs to ensure good API documentation.

## Additional Tools

Consider installing these useful tools:

- `rust-analyzer` - Language server for IDE support
- `cargo-watch` - Automatically rebuild and rerun tests on file changes:
  ```bash
  cargo install cargo-watch
  cargo watch -x test
  ```

## Troubleshooting

- If you encounter compilation errors, check the error messages carefully - Rust's compiler often provides helpful suggestions
- Ensure you're using an up-to-date version of Rust: `rustup update`
- For dependency issues, try `cargo update` to update your Cargo.lock file