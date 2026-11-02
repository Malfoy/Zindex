#[path = "../utils/scripts/build_version.rs"]
mod build_version;

fn main() {
    println!("cargo:rerun-if-changed=.cargo/build.rs");
    build_version::emit(std::path::Path::new(env!("CARGO_MANIFEST_DIR")));
    let features = std::env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    let has_simd = features.split(',').any(|f| f == "avx2" || f == "neon");
    if !has_simd {
        println!(
            "cargo:warning=Neither AVX2 nor NEON is enabled for this target; \
             building with scalar minimizer fallback (slower)."
        );
    }
}
