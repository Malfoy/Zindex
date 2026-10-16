fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let features = std::env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    let has_simd = features.split(',').any(|f| f == "avx2" || f == "neon");
    if !has_simd {
        println!(
            "cargo:warning=Neither AVX2 nor NEON is enabled for this target; \
             building with scalar minimizer fallback (slower)."
        );
    }
}
