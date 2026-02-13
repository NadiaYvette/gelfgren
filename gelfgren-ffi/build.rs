//! Build script to generate C headers using cbindgen.

use std::env;
use std::path::PathBuf;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let package_name = env::var("CARGO_PKG_NAME").unwrap();
    let output_file = PathBuf::from(&crate_dir)
        .join("../include")
        .join("gelfgren.h");

    // Ensure output directory exists
    if let Some(parent) = output_file.parent() {
        std::fs::create_dir_all(parent).expect("Failed to create include directory");
    }

    println!("cargo:rerun-if-changed=src/");
    println!("cargo:rerun-if-changed=cbindgen.toml");

    match cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_config(cbindgen::Config::from_file("cbindgen.toml").unwrap())
        .generate()
    {
        Ok(bindings) => {
            bindings.write_to_file(&output_file);
            println!("cargo:info=Generated C header at {:?}", output_file);
        }
        Err(e) => {
            println!("cargo:warning=Unable to generate bindings: {}", e);
            println!("cargo:warning=Continuing without header generation");
        }
    }
}
