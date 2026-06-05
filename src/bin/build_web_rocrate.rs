//! Build the top-level RO-Crate for a hosted web/ directory.
//!
//! Scans `<web-dir>/data/*/ro-crate-metadata.json`, hoists each
//! per-dataset crate's structured properties into a single
//! `<web-dir>/ro-crate-metadata.json`, and writes the page-level
//! WebApplication / Person / SoftwareSourceCode identity entities
//! into the same graph. Used by the Pages workflow after dataset
//! generation so the deployed bundle has a fresh top-level crate
//! reflecting whatever datasets are present.
//!
//! Usage:
//!
//! ```sh
//! cargo run --release --bin build_web_rocrate -- \
//!     --web-dir web/ratdb/ \
//!     --page-url https://ratdb.app.pirogov.de/
//! ```

use clap::Parser;
use tilezz::stringmatch::dafsa::write_collection_ro_crate;

#[derive(Parser, Debug)]
#[command(about = "Build the top-level RO-Crate for a hosted web/ directory.")]
struct Cli {
    /// Path to the web root (must already exist; will contain
    /// `data/<asset>/ro-crate-metadata.json` files produced by
    /// `rat_enum --mode dafsa-blocks`).
    #[arg(long)]
    web_dir: std::path::PathBuf,

    /// Canonical hosted URL of the page, surfaced into the
    /// WebApplication entity's `url` property.
    #[arg(long, default_value = "https://apirogov.github.io/tilezz/")]
    page_url: String,
}

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    if !cli.web_dir.is_dir() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("--web-dir {:?} is not a directory", cli.web_dir),
        ));
    }
    write_collection_ro_crate(&cli.web_dir, &cli.page_url)?;
    println!(
        "wrote {}",
        cli.web_dir.join("ro-crate-metadata.json").display()
    );
    Ok(())
}
