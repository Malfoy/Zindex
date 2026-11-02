#[path = "../../scripts/build_version.rs"]
mod build_version;

fn main() {
    build_version::emit(&std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../.."));
}
