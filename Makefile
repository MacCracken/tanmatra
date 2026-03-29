.PHONY: check fmt clippy test audit deny bench coverage build doc clean examples

# Run all CI checks locally
check: fmt clippy test audit deny

# Format check
fmt:
	cargo fmt --all -- --check

# Lint (zero warnings)
clippy:
	cargo clippy --all-features --all-targets -- -D warnings

# Run test suite
test:
	cargo test --all-features
	cargo test --no-default-features

# Security audit
audit:
	cargo audit

# Supply-chain checks (cargo-deny)
deny:
	cargo deny check

# Run benchmarks
bench:
	cargo bench

# Generate coverage report
coverage:
	cargo tarpaulin --all-features --skip-clean

# Build release
build:
	cargo build --release --all-features

# Generate documentation
doc:
	RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features

# Run examples
examples:
	cargo run --example basic
	cargo run --example nuclear
	cargo run --example spectral
	cargo run --example relativity

# Clean build artifacts
clean:
	cargo clean
	rm -rf coverage/
