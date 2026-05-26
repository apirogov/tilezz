//! `ProfileGuard`: tiny RAII helper for binary flamegraph profiling.
//!
//! Each CLI binary that wants a `--pprof <PATH>` option follows the same
//! pattern: start a `pprof::ProfilerGuard` at the top of `main`, write a
//! flamegraph SVG to the requested path before exit, and silently no-op
//! when the `pprof` cargo feature is disabled. This helper centralizes
//! that pattern so the binaries stay free of `#[cfg(feature = "pprof")]`
//! ladders.
//!
//! Usage:
//! ```ignore
//! let guard = ProfileGuard::start(cli.pprof.as_deref());
//! // ... work ...
//! guard.finish();
//! ```

#[cfg(feature = "pprof")]
type Inner = Option<(pprof::ProfilerGuard<'static>, String)>;

#[cfg(not(feature = "pprof"))]
type Inner = ();

pub struct ProfileGuard {
    #[allow(dead_code)] // unused when `pprof` feature is off
    inner: Inner,
}

impl ProfileGuard {
    /// Start a profiler. If `path` is `None`, the guard is a no-op. If
    /// the `pprof` feature is disabled, a warning is printed (when a
    /// path was requested) and the guard is a no-op.
    pub fn start(path: Option<&str>) -> Self {
        #[cfg(feature = "pprof")]
        {
            let inner = path.map(|p| {
                let g = pprof::ProfilerGuardBuilder::default()
                    .frequency(1000)
                    .blocklist(&["libc", "libgcc", "pthread", "vdso"])
                    .build()
                    .expect("profiler start");
                (g, p.to_string())
            });
            Self { inner }
        }
        #[cfg(not(feature = "pprof"))]
        {
            if path.is_some() {
                eprintln!("warning: --pprof requires building with `--features pprof`");
            }
            Self { inner: () }
        }
    }

    /// Write the flamegraph (if one was requested and the feature is
    /// enabled) and consume the guard. Any I/O error is reported to
    /// stderr but does not panic; profiling failures should not break
    /// the binary's primary output.
    pub fn finish(self) {
        #[cfg(feature = "pprof")]
        if let Some((g, path)) = self.inner {
            match g.report().build() {
                Ok(report) => match std::fs::File::create(&path) {
                    Ok(file) => {
                        if let Err(e) = report.flamegraph(file) {
                            eprintln!("flamegraph write failed: {e}");
                        } else {
                            eprintln!("wrote flamegraph: {path}");
                        }
                    }
                    Err(e) => eprintln!("flamegraph create failed: {e}"),
                },
                Err(e) => eprintln!("profiler report failed: {e}"),
            }
        }
    }
}
