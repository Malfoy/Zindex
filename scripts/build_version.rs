use std::path::Path;
use std::process::Command;

fn git(root: &Path, args: &[&str]) -> Option<String> {
    let result = Command::new("git")
        .current_dir(root)
        .args(args)
        .output()
        .ok()?;
    result
        .status
        .success()
        .then(|| String::from_utf8_lossy(&result.stdout).trim().to_owned())
}

pub fn emit(root: &Path) {
    println!(
        "cargo:rerun-if-changed={}",
        root.join("scripts/build_version.rs").display()
    );
    let package = std::env::var("CARGO_PKG_VERSION").unwrap();
    let version = if git(root, &["rev-parse", "--show-toplevel"])
        .and_then(|p| std::fs::canonicalize(p).ok())
        == std::fs::canonicalize(root).ok()
        && root.join(".git").exists()
    {
        // Editing tracked sources after a clean build must refresh the dirty
        // indicator even when HEAD has not moved.
        if let Some(files) = git(root, &["ls-files"]) {
            for path in files.lines() {
                println!("cargo:rerun-if-changed={}", root.join(path).display());
            }
        }
        for name in ["HEAD", "index", "packed-refs"] {
            if let Some(path) = git(root, &["rev-parse", "--git-path", name]) {
                let path = root.join(path);
                if path.exists() {
                    println!("cargo:rerun-if-changed={}", path.display());
                }
            }
        }
        if let Some(reference) = git(root, &["symbolic-ref", "-q", "HEAD"]) {
            if let Some(path) = git(root, &["rev-parse", "--git-path", &reference]) {
                println!("cargo:rerun-if-changed={}", root.join(path).display());
            }
        }
        let sha =
            git(root, &["rev-parse", "--short=12", "HEAD"]).unwrap_or_else(|| "unknown".into());
        let dirty = match git(root, &["status", "--porcelain", "--untracked-files=normal"]) {
            Some(status) if status.is_empty() => "",
            Some(_) => ", dirty",
            None => ", state unknown",
        };
        format!("{package} (git {sha}{dirty})")
    } else {
        format!("{package} (source archive; Git revision unavailable)")
    };
    println!("cargo:rustc-env=ZINDEX_BUILD_VERSION={version}");
}
